#!/usr/bin/env python

from __future__ import print_function   # for python 2.7 compatibility

import os
import time
import argparse
import numpy as np
import tables

import pyzfp

import twixtools
import twixtools.mdh_def as mdh_def
import twixtools.hdr_def as hdr_def

import cProfile
import pstats

# helper functions:

def to_freqdomain(data, x_in_timedomain):
    if not x_in_timedomain:
        return data, False
    else:
        return np.fft.ifft(data), False


def to_timedomain(data, x_in_timedomain):
    if x_in_timedomain:
        return data, True
    else:
        return np.fft.fft(data), True


def reduce_data(data, mdh, remove_os=False, cc_mode=False, mtx=None):

    if data.dtype == np.dtype("S1"):
        # nothing to do in case of bytearray
        return data, False, False 

    x_in_timedomain = True

    rm_os_active = remove_os
    if mdh_def.is_flag_set(mdh, 'NOISEADJSCAN'):
        rm_os_active = False

    cc_active = False
    if cc_mode and mtx is not None and data.shape[0]==mtx.shape[-1]:
        cc_active = True

    reflect_data = False
    if rm_os_active or (cc_active and cc_mode=='gcc'):
        reflect_data = bool(mdh['aulEvalInfoMask'][0] & (1 << 24))
        if reflect_data:
            data = data[:,::-1]
    
    if rm_os_active:
        nx = data.shape[-1]
        data, x_in_timedomain = to_freqdomain(data, x_in_timedomain)
        data = np.delete(data, slice(nx//4, nx*3//4), -1)

    if cc_active:
        nc, nx = data.shape
        ncc = mtx.shape[1]
        if cc_mode == 'scc':
            data = mtx[0] @ data
        elif cc_mode == 'gcc':
            if nx!=mtx.shape[0]:
                # nx mismatch; deactivate cc mode
                cc_active = False
            else:
                data, x_in_timedomain = to_freqdomain(data, x_in_timedomain)
                for x in range(nx):
                    data[:ncc,x] = mtx[x] @ data[:,x]
                data = data[:ncc,:]

    data, x_in_timedomain = to_timedomain(data, x_in_timedomain)
    
    if reflect_data:
        data = data[:,::-1]

    return np.complex64(data), rm_os_active, cc_active


def expand_data(data, mdh, remove_os=False, cc_mode=False, inv_mtx=None):
      
    if data.dtype == np.dtype("S1"):
        return data # nothing to do in case of bytearray
    
    if inv_mtx is None:
        inv_mtx = False

    nc, nx = data.shape

    x_in_timedomain = True
    reflect_data = False
    if remove_os or cc_mode=='gcc':
        # for performance reasons, x dim was stored in freq. domain
        reflect_data = bool(mdh['aulEvalInfoMask'][0] & (1 << 24))
        if reflect_data:
            data = data[:,::-1]
    
    if cc_mode:
        nc = inv_mtx.shape[1]
        ncc = inv_mtx.shape[-1]
        if cc_mode=='scc':
            try:
                data = inv_mtx[0] @ data
            except:
                print('error during inv_mtx @ data')
                print('mdh flags:', mdh_def.get_active_flags(mdh))
                print('data shape: ', data.shape)
                print('inv_mtx shape: ', inv_mtx.shape)
        else: # 'gcc'
            if nx!=inv_mtx.shape[0]:
                # nx mismatch; deactivate cc mode
                cc_active = False
            else:
                data, x_in_timedomain = to_freqdomain(data, x_in_timedomain)
                # pad missing channels in data with zeros
                data = np.pad(data, [(0, nc-ncc), (0, 0)])
                for x in range(nx):
                    data[:,x] = inv_mtx[x,:,:] @ data[:ncc,x]

    if remove_os:
        data, x_in_timedomain = to_freqdomain(data, x_in_timedomain)
        data = np.insert(data, nx//2, np.zeros((nx, 1), dtype=data.dtype), -1)
        
    if reflect_data:
        data = data[:,::-1]

    data, x_in_timedomain = to_timedomain(data, x_in_timedomain)

    return np.complex64(data)


def get_restrictions(mdh_flags):
    flags = {key for key, item in mdh_flags.items() if item}

    restrictions = 'NO_RESTRICTIONS'
    if flags.intersection({'ACQEND', 'SYNCDATA'}):
        restrictions = 'BYTEARRAY'
    elif flags.intersection({'NOISEADJSCAN', 'noname60'}):
        restrictions = 'NO_COILCOMP'

    return restrictions


def calculate_prewhitening(noise, scale_factor=1.0):
    '''Calculates the noise prewhitening matrix

    :param noise: Input noise data (array or matrix), ``[coil, nsamples]``
    :scale_factor: Applied on the noise covariance matrix. Used to
                   adjust for effective noise bandwith and difference in
                   sampling rate between noise calibration and actual measurement:
                   scale_factor = (T_acq_dwell/T_noise_dwell)*NoiseReceiverBandwidthRatio

    :returns w: Prewhitening matrix, ``[coil, coil]``, w*data is prewhitened
    '''

    noise_int = noise.reshape((noise.shape[0], -1))
    M = float(noise_int.shape[1])
    dmtx = (1/(M-1))*(noise_int.dot(np.conj(noise_int).T))
    mtx = np.linalg.cholesky(dmtx)/np.sqrt(2*scale_factor)
    dmtx = np.linalg.inv(mtx)
    return dmtx, mtx


def scc_calibrate_mtx(data):
    [ny, nc, nx] = np.shape(data)
    data = np.moveaxis(data,1,-1)
    data = data.flatten().reshape((-1, nc))
    U, s, V = np.linalg.svd(data, full_matrices=False)
    V = np.conj(V)
    np.save('scc_mtx.npy', s)
    return V[np.newaxis, :, :], s


def gcc_calibrate_mtx(data):
    [ny, nc, nx] = np.shape(data)
    data = np.moveaxis(data,1,-1)
    im = np.fft.ifft(data, axis=1)
    mtx = np.zeros((nx, nc, nc),dtype = 'complex64')
    s = np.zeros((nx, nc), dtype='float32')
    for x in range(nx):
        U, s[x], V = np.linalg.svd(im[:,x,:], full_matrices=False)
        mtx[x,:,:] = np.conj(V) 
    np.save('gcc_mtx.npy', s)
    return mtx, s.mean(axis=0)


def calibrate_mtx(data, cc_mode):
    if cc_mode == 'scc':
        mtx, s = scc_calibrate_mtx(data)
    elif cc_mode == 'gcc':
        mtx, s = gcc_calibrate_mtx(data)
    else:
        print('unknown cc_mode "%s"'%(cc_mode))
        raise ValueError
    return mtx, s


def get_cal_data(meas, cc_mode, remove_os):
    cal_list, cal_isima, cal_nx, cal_nc = list(), list(), list(), list()
    for mdb_key, mdb in enumerate(meas['mdb']):
        restrictions =  get_restrictions(mdb.get_flags())
        if restrictions=='NO_RESTRICTIONS': # and not mdb.is_flag_set('RAWDATACORRECTION'): # wip!?
            cal_list.append(int(mdb_key))
            cal_isima.append(mdh_def.is_image_scan(mdb.mdh))
            cal_nx.append(int(mdb.data.shape[-1]))
            cal_nc.append(int(mdb.data.shape[0]))
    
    mtx = None
    max_calib_samples = 4096 # wip, higher is better
    if len(cal_list)>0:
        # make sure that all blocks have same read size and same number of coils
        import scipy.stats
        cal_list = np.asarray(cal_list)
        cal_isima = np.asarray(cal_isima)
        cal_nx = np.asarray(cal_nx)
        cal_nc = np.asarray(cal_nc)
        if cal_isima.sum()>=16:
            mode_x = scipy.stats.mode(cal_nx[cal_isima]).mode[0]
            mode_c = scipy.stats.mode(cal_nc[cal_isima]).mode[0]
        else:
            mode_x = scipy.stats.mode(cal_nx).mode[0]
            mode_c = scipy.stats.mode(cal_nc).mode[0]
        mask = (cal_nx==mode_x) & (cal_nc==mode_c)
        cal_list = cal_list[mask]
        cal_isima = cal_isima[mask]
        # first pick image scans
        if cal_list.size<=max_calib_samples:
            # pick all
            pick_mask = np.ones(cal_list.shape, dtype=bool)
        else:
            # pick all image scans
            pick_mask = cal_isima.copy()
            missing_picks = max_calib_samples - pick_mask.sum()
            print('missing picks: ', missing_picks)
            if missing_picks>max_calib_samples//4:
                # supplement with non-image scans
                tmp = np.full(pick_mask.size()-pick_mask.sum(), False)
                tmp[:missing_picks] = True
                np.random.shuffle(tmp)
                pick_mask[~pick_mask] = tmp
            elif missing_picks<0:
                # pick random subset of image scans
                tmp = np.full(pick_mask.sum(), False)
                tmp[:max_calib_samples] = True
                np.random.shuffle(tmp)
                pick_mask[pick_mask] = tmp

        cal_list = cal_list[pick_mask].tolist()
        if remove_os:
            cal_data = np.zeros((len(cal_list), mode_c, mode_x//2), dtype=np.complex64)
        else:
            cal_data = np.zeros((len(cal_list), mode_c, mode_x), dtype=np.complex64)
        for cnt, mdb_idx in enumerate(cal_list):
            data = meas['mdb'][mdb_idx].data
            if meas['mdb'][mdb_idx].is_flag_set('REFLECT'):
                data = data[:,::-1]
            if remove_os:
                data = np.fft.ifft(data)
                nx = data.shape[-1]
                data = np.delete(data, slice(nx//4, nx*3//4), -1)
                data = np.fft.fft(data)
            cal_data[cnt, :, :] = data

    return cal_data


def compress_twix(infile, outfile, remove_os=False, cc_mode=False, ncc=None, cc_tol=0.05, zfp=False, zfp_tol=1e-5, zfp_prec=None, rm_fidnav=False):
    
    twix = twixtools.read_twix(infile)

    filters = tables.Filters(complevel=5, complib='zlib') # lossless compression settings

    
    mtx = None
    noise_mtx = None
    noise_dmtx = None
    if cc_mode or zfp:
        # # calibrate noise decorrelation matrix for better compression
        # noise = list()
        # for mdb in twix[1]['mdb']:
        #     if mdb.is_flag_set('NOISEADJSCAN'):
        #         noise.append(mdb.data)
        # if len(noise)>0:
        #     noise_dmtx, noise_mtx = calculate_prewhitening(np.asarray(noise).swapaxes(0,1))
        # del(noise)
        pass

    if cc_mode:
        # calibrate coil compression based on last scan in list (image scan)
        # use the calibration coil weights for all data that fits
        cal_data = get_cal_data(twix[-1], cc_mode, remove_os)
        mtx, s = calibrate_mtx(cal_data, cc_mode)
        del(cal_data)
        if ncc is None:
            ncc = 1 + np.argwhere(np.cumsum(s)/s.sum() > (1-cc_tol))[0,0]
        mtx = mtx[:,:ncc, :]
        print('coil compression from %d channels to %d virtual channels'%(mtx.shape[-1], ncc))
    
    t_start = time.time()
    with tables.open_file(outfile, mode="w") as f:
        
        
        f.root._v_attrs.original_filename = os.path.basename(infile)
        f.root._v_attrs.cc_mode = cc_mode
        f.root._v_attrs.zfp = zfp
        if zfp_tol is None:
            f.root._v_attrs.zfp_tol = -1
        else:
            f.root._v_attrs.zfp_tol = zfp_tol
        if zfp_prec is None:
            f.root._v_attrs.zfp_prec = -1
        else:
            f.root._v_attrs.zfp_prec = zfp_prec
            
        f.create_carray(f.root,"multi_header", obj=np.frombuffer(twix[0].tobytes(), 'S1'), filters=filters)
        
        if mtx is not None:
            # save mtx for coil compression
            f.create_carray(f.root,"mtx", obj=mtx, filters=filters)
        if noise_dmtx is not None:
            f.create_carray(f.root,"noise_dmtx", obj=noise_dmtx, filters=filters)
            f.create_carray(f.root,"noise_mtx", obj=noise_mtx, filters=filters)
        
        scanlist = []
        for meas_key, meas in enumerate(twix[1:]):
            scanlist.append("scan%d"%(meas_key))
            grp = f.create_group("/","scan%d"%(meas_key))
            f.create_carray(grp,"hdr_str", obj=meas['hdr_str'], filters=filters)

            # remove fidnav scans if necessary
            if rm_fidnav:
                for mdb_key, mdb in enumerate(meas['mdb']):
                    if mdb.is_flag_set('noname60'):
                        del(meas['mdb'][mdb_key])

            mdh_count = len(meas['mdb'])
            mdh_size = mdh_def.scan_hdr_type.itemsize            
            ch_hdr_size = mdh_def.channel_hdr_type.itemsize
            
            f.create_carray(grp,"mdh_info", shape=[mdh_count*mdh_size], atom=tables.UInt8Atom(), filters=filters)
            # not all mdh's have coil_info's (namely ACQEND & SYNCDATA) but for simplicity allocate space anyway (only a few bytes overhead)
            f.create_carray(grp,"coil_info", shape=[mdh_count*ch_hdr_size], atom=tables.UInt8Atom(), filters=filters)
            # similarly, just allocate the maximum number of possible coils (64 - in special cases 128 - increase further!?)
            f.create_carray(grp,"coil_list", shape=[mdh_count, 64], atom=tables.UInt8Atom(), filters=filters)
            # create list to track for which mdbs os removal is active
            f.create_carray(grp,"rm_os_active", obj=np.zeros(mdh_count, dtype=bool), filters=filters)
            # create list to track which mdbs have been coil compressed
            f.create_carray(grp,"cc_active", obj=np.zeros(mdh_count, dtype=bool), filters=filters)

            dt = tables.UInt8Atom(shape=())
            if zfp:
                f.create_vlarray(grp,"DATA", atom=dt)
            else:
                f.create_vlarray(grp,"DATA", atom=dt, filters=filters)
            
            syncscans = 0
            for mdb_key, mdb in enumerate(meas['mdb']):
                is_syncscan = mdb.is_flag_set('SYNCDATA')
                if rm_fidnav: # we have to update the scan counters
                    if not is_syncscan:
                        mdb.mdh['ulScanCounter'] = mdb_key + 1 - syncscans # scanCounter starts at 1
                    else:
                        syncscans += 1

                # write mdh
                grp.mdh_info[mdb_key*mdh_size:(mdb_key+1)*mdh_size] = np.frombuffer(mdb.mdh, dtype = 'uint8')
                
                if is_syncscan or mdb.is_flag_set('ACQEND'):
                    data = np.ascontiguousarray(mdb.data).view('uint8')
                else:
                    restrictions = get_restrictions(mdb.get_flags())
                    if restrictions == 'NO_COILCOMP':
                        data, grp.rm_os_active[mdb_key],_ = reduce_data(mdb.data, mdb.mdh, remove_os, cc_mode=False)
                    else:
                        data, grp.rm_os_active[mdb_key], grp.cc_active[mdb_key] = reduce_data(mdb.data, mdb.mdh, remove_os, cc_mode=cc_mode, mtx=mtx)

                    data = data.flatten()
                    if zfp:
                        data = pyzfp.compress(data.view('float32'), tolerance=zfp_tol, precision=zfp_prec, parallel=True)
                        data = np.frombuffer(data, dtype = 'uint8')
                    else:
                        data = data.view('uint8')
                    if len(mdb.channel_hdr) > 0:
                        mdb.channel_hdr[0]['ulScanCounter'] = mdb.mdh['ulScanCounter']
                        grp.coil_info[mdb_key*ch_hdr_size:(mdb_key+1)*ch_hdr_size] = np.frombuffer(mdb.channel_hdr[0], dtype = 'uint8')
                        for coil_key, coil_item in enumerate(mdb.channel_hdr):
                            grp.coil_list[mdb_key, coil_key] = coil_item['ulChannelId']

                # write data
                grp.DATA.append(data)
        
        f.root._v_attrs.scanlist = scanlist    
            
            # from joblib import Parallel, delayed
            # Parallel(n_jobs=2)(delayed(task)(mdb_key, mdb, is_byte, count, grp, remove_os, zfp, zfp_tol, zfp_prec, mtx) for mdb_key, (mdb, is_byte, count) in enumerate(zip(meas['mdb'], is_bytearray, data_counter)))
            
            
                    

    elapsed_time = (time.time() - t_start)
    print("compression finished in %d:%02d:%02d h"%(elapsed_time//3600, (elapsed_time%3600)//60, elapsed_time%60))
    print("compression factor = %.2f"%(os.path.getsize(infile)/os.path.getsize(outfile)))


def reconstruct_twix(infile, outfile=None):
    #wip: function takes no parameters, all necessary information needs to be included in hdf file
    def write_sync_bytes(f):
        syncbytes = (512-f.tell()%512)%512
        # print('syncbytes', syncbytes)
        f.write(b'\x00' * syncbytes)

    if outfile is None:
        with tables.open_file(infile, mode="r") as f:
            outfile = f.root._v_attrs.original_filename
    
    t_start = time.time()

    with tables.open_file(infile, mode="r") as f, open(outfile, 'wb') as fout:
        
        cc_mode = f.root._v_attrs.cc_mode
        zfp = f.root._v_attrs.zfp
        zfp_tol = f.root._v_attrs.zfp_tol
        if zfp_tol<0:
            zfp_tol = None
        zfp_prec = f.root._v_attrs.zfp_prec
        if zfp_prec<0:
            zfp_prec = None
        
        inv_mtx = None
        if cc_mode is not False and hasattr(f.root,'mtx'):
            mtx = f.root.mtx[()]

            inv_mtx = np.zeros_like(mtx).swapaxes(1,-1)
            for x in range(mtx.shape[0]):
                inv_mtx[x,:,:] = np.linalg.pinv(mtx[x,:,:])
            del(mtx)
        
        # allocate space for multi-header
        fout.write(b'\x00' * 10240)
        
        scan_pos = list()
        scan_len = list()
        
        scanlist = f.root._v_attrs.scanlist
        for key, scan in enumerate(scanlist):
            # keep track of byte pos
            scan_pos.append(fout.tell())

            # get compression info
            rm_os_active = getattr(f.root,scan).rm_os_active[()]
            cc_active = getattr(f.root,scan).cc_active[()]

            # write header
            getattr(f.root,scan).hdr_str[()].tofile(fout)
                
            mdh_info = np.frombuffer(getattr(f.root,scan).mdh_info[()],dtype=mdh_def.scan_hdr_type)
            for mdh_key, raw_mdh in enumerate(mdh_info):

                # write mdh
                mdh = np.frombuffer(raw_mdh, mdh_def.scan_hdr_type)[0]
                mdh.tofile(fout)

                # write data
                is_bytearray = mdh_def.is_flag_set(mdh, 'ACQEND') or mdh_def.is_flag_set(mdh, 'SYNCDATA')
                if is_bytearray:
                    data = getattr(f.root,scan).DATA[mdh_key]
                    # data = f[scan]['BYTEARRAY'][data_count['BYTEARRAY']]
                    data.tofile(fout)
                else:
                    n_sampl = mdh['ushSamplesInScan']
                    n_coil = mdh['ushUsedChannels']
                    n_data_sampl = n_sampl
                    if rm_os_active[mdh_key]:
                        n_data_sampl //= 2
                    n_data_coils = n_coil
                    if cc_mode and cc_active[mdh_key]:
                        n_data_coils = inv_mtx.shape[-1]

                    data = getattr(f.root,scan).DATA[mdh_key]
                    if zfp:
                        data = memoryview(data)
                        data = pyzfp.decompress(data, [n_data_coils*2*n_data_sampl], np.dtype('float32'), tolerance=zfp_tol, precision=zfp_prec)
                    
                    data = np.ascontiguousarray(data).view('complex64')
                    data = data.reshape(n_data_coils, n_data_sampl)
                    
                    if cc_mode and cc_active[mdh_key]:
                        data = expand_data(data, mdh, rm_os_active[mdh_key], cc_mode=cc_mode, inv_mtx=inv_mtx)
                    else:
                        data = expand_data(data, mdh, rm_os_active[mdh_key], cc_mode=False)

                    data = data.reshape((n_coil, -1))
                    
                    coil_info = np.frombuffer(getattr(f.root,scan).coil_info[()],dtype=mdh_def.channel_hdr_type)
                    coil_hdr = np.frombuffer(coil_info[mdh_key], mdh_def.channel_hdr_type)[0].copy()
                    coil_idx = getattr(f.root,scan).coil_list[mdh_key]
                    for c in range(data.shape[0]):
                        #write channel header
                        coil_hdr['ulChannelId'] = coil_idx[c]
                        coil_hdr.tofile(fout)
                        # write data
                        data[c].tofile(fout)

            # update scan_len
            scan_len.append(fout.tell() - scan_pos[-1])

            # add sync bytes between scans
            write_sync_bytes(fout)

        # now write preallocated MultiRaidFileHeader
        n_scans = len(scan_pos)
        multi_header = np.frombuffer(f.root.multi_header[()], hdr_def.MultiRaidFileHeader)[0]
        # write NScans
        multi_header['hdr']['count_'] = n_scans
        # write scan_pos & scan_len for each scan
        for i, (pos_, len_) in enumerate(zip(scan_pos, scan_len)):
            # print('scan', i, ' len_ old: ', multi_header['entry'][i]['len_'], ' new:', len_)
            # print('scan', i, ' off_ old: ', multi_header['entry'][i]['off_'], ' new:', pos_)
            multi_header['entry'][i]['len_'] = len_
            multi_header['entry'][i]['off_'] = pos_

        # write MultiRaidFileHeader
        fout.seek(0)
        multi_header.tofile(fout)
    
    elapsed_time = (time.time() - t_start)
    print("decompression finished in %d:%02d:%02d h"%(elapsed_time//3600, (elapsed_time%3600)//60, elapsed_time%60))


class TwixFile(argparse.FileType):
    def __call__(self, string):
        base, ext = os.path.splitext(string)

        if ext == '':
            string = string + '.dat'  # .dat is default file extension
        else:
            if (str.lower(ext) != '.dat'):
                parser.error('twix file %s should have a .dat extension' % (string))

        returnFile = super(TwixFile, self).__call__(string)
        returnFile.close()
        returnFile = os.path.abspath(returnFile.name)
        return returnFile


class HDF5File(argparse.FileType):
    def __call__(self, string):
        base, ext = os.path.splitext(string)

        if ext == '':
            string = string + '.h5'  # .dat is default file extension
        else:
            if (str.lower(ext) != '.h5' and str.lower(ext) != '.hdf5'):
                parser.error('hdf5 file %s should have a .h5 extension' % (string))

        returnFile = super(HDF5File, self).__call__(string)
        returnFile.close()
        returnFile = os.path.abspath(returnFile.name)
        return returnFile

def md5(fname):
    import hashlib
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compresses and decompresses Siemens raw data files (.dat)")
    
    parser.add_argument("-d", "--decompress", action="store_true")

    args, rem_args = parser.parse_known_args()
    
    if args.decompress:
        parser.add_argument('-i', '--infile', '--in', type=HDF5File('r'),
                            help='Input HDF5 file', required=True)
        parser.add_argument('-o', '--outfile', '--out', type=TwixFile('x'),
                            help='Output twix .dat file', required=False)
    else:
        parser.add_argument('-i', '--infile', '--in', type=TwixFile('r'),
                            help='Input twix .dat file', required=True)
        parser.add_argument('-o', '--outfile', '--out', type=HDF5File('w'),
                            help='Output HDF5 file', required=True)
        parser.add_argument("--remove_os", action="store_true")
        parser.add_argument("--remove_fidnav", action="store_true")
        
        group_cc = parser.add_mutually_exclusive_group()
        group_cc.add_argument("--scc", action="store_true")
        group_cc.add_argument("--gcc", action="store_true")

        group_ncc = parser.add_mutually_exclusive_group()
        group_ncc.add_argument("-n", "--ncc", "-n_compressed_coils", type=int)
        group_ncc.add_argument("-t", "--cc_tol", default=0.05, type=float)

        parser.add_argument("--zfp", action="store_true")
        group_zfp = parser.add_mutually_exclusive_group()
        group_zfp.add_argument("--zfp_tol", default=1e-7, type=float)
        group_zfp.add_argument("--zfp_prec", default=None, type=float)

        parser.add_argument("--testmode", action="store_true")
        parser.add_argument("--profile", action="store_true")

    args = parser.parse_args()

    if args.decompress:
        if args.profile:
            cProfile.run('reconstruct_twix(args.infile, args.outfile)','stats_decompr')
            p = pstats.Stats('stats_decompr')
            p.strip_dirs().sort_stats('cumulative').print_stats(15)  
        else:
            reconstruct_twix(args.infile, args.outfile)
    else:
        if args.ncc is not None:
            args.cc_tol = None

        if args.zfp_prec is not None:
            args.zfp_tol = None

        cc_mode = False
        if args.scc:
            cc_mode = 'scc'
        elif args.gcc:
            cc_mode = 'gcc'
    
        if args.profile:
            cProfile.run('compress_twix(args.infile, args.outfile, remove_os=args.remove_os, cc_mode=cc_mode, ncc=args.ncc, cc_tol=args.cc_tol, zfp=args.zfp, zfp_tol=args.zfp_tol, zfp_prec=args.zfp_prec, rm_fidnav=args.remove_fidnav)','stats_compr')
            p = pstats.Stats('stats_compr')
            p.strip_dirs().sort_stats('cumulative').print_stats(15)    
        else:    
            compress_twix(args.infile, args.outfile, remove_os=args.remove_os, cc_mode=cc_mode, ncc=args.ncc, cc_tol=args.cc_tol, zfp=args.zfp, zfp_tol=args.zfp_tol, zfp_prec=args.zfp_prec, rm_fidnav=args.remove_fidnav)
        
        if args.testmode:
            with tables.open_file(args.outfile, mode="r") as f:
                for group in f.walk_groups():
                    print(group)
                for group in f.walk_groups("/scan0"):
                    for array in f.list_nodes(group, classname='Array'):
                        print(array)    

            # import re
            # out_name = re.search('meas_MID(\d+)_FID(\d+)', args.infile)
            
            # if out_name is None:
            #     import tempfile
            #     out_name = tempfile.mktemp(suffix='.dat')
            # else:
            #     out_name = out_name.group()
            out_name = os.path.splitext(args.infile)[0]
            if args.remove_fidnav:
                out_name += '_rmfidnav'
            if args.remove_os:
                out_name += '_rmOS'
            if args.zfp:
                if args.zfp_tol:
                    out_name += '_zfptol' + str(args.zfp_tol)
                else:
                    out_name += '_zfpprec' + str(args.zfp_prec)
            if args.scc:
                if args.ncc:
                    out_name += '_scc' + str(args.ncc)
                else:
                    out_name += '_scctol' + str(args.cc_tol)
            elif args.gcc:
                if args.ncc:
                    out_name += '_gcc' + str(args.ncc)
                else:
                    out_name += '_gcctol' + str(args.cc_tol)
            
            if out_name == os.path.splitext(args.infile)[0]:
                out_name += '_new'
            
            out_name += '.dat'
            if args.profile:
                cProfile.run('reconstruct_twix(args.outfile, out_name)','stats_decompr')
                p = pstats.Stats('stats_decompr')
                p.strip_dirs().sort_stats('cumulative').print_stats(15)  
            else:
                reconstruct_twix(args.outfile, out_name)
            inputsz = os.path.getsize(args.infile)
            comprsz = os.path.getsize(args.outfile)
            reconsz = os.path.getsize(out_name)
            print('original size = ', inputsz, ' compressed size =', comprsz, ' comp. factor =', inputsz/comprsz, ' reconstructed size =', reconsz)
            # if md5(args.infile) == md5(out_name):
            #     print('md5 hashes match, perfect reconstruction (lossless compression)')
            # # test read:
            # try:
            #     twix = twixtools.read_twix(out_name)
            #     print('\nreconstructed twix file successfully parsed')
            # except:
            #     print('\nerror parsing reconstructed twix file')
            # os.remove(out_name)
            print('reconstructed .dat file:', out_name)
