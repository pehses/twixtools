#!/usr/bin/env python

from __future__ import print_function   # for python 2.7 compatibility

import os
import time
import argparse
import numpy as np
import h5py
import pyzfp

import twixtools
import twixtools.mdh_def as mdh_def
import twixtools.hdr_def as hdr_def

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
    
    nc, nx = data.shape
    x_in_timedomain = True
    
    if remove_os or cc_mode=='gcc':
        reflect_data = bool(mdh['aulEvalInfoMask'][0] & (1 << 24))
        if reflect_data:
            data = data[:,::-1]
    else:
        reflect_data = False
    
    if cc_mode and inv_mtx is not None:
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


def create_restriction_categories(cc_mode=False):
    restriction_categories = {'NO_RESTRICTIONS': None}  # all data that does not fit into other categories
    restriction_categories['BYTEARRAY'] = {'ACQEND', 'SYNCDATA'}
    if cc_mode is not None and cc_mode is not False:
        restriction_categories['NO_COILCOMP'] = set()
        # restriction_categories['NO_COILCOMP'].add('RAWDATACORRECTION')
        if cc_mode == 'gcc':
            restriction_categories['NO_COILCOMP'].add('NOISEADJSCAN')
            restriction_categories['NO_COILCOMP'].add('noname60') # our own fidnav scan
    else:
        restriction_categories['NO_COILCOMP'] = ''
    print(restriction_categories)
    return restriction_categories


def get_restrictions(mdh_flags, restriction_categories):
    flags = {key for key, item in mdh_flags.items() if item}
    restrictions = 'NO_RESTRICTIONS'
    if flags.intersection(restriction_categories['BYTEARRAY']):
        restrictions = 'BYTEARRAY'
    elif flags.intersection(restriction_categories['NO_COILCOMP']):
        restrictions = 'NO_COILCOMP'
    if restrictions=='BYTEARRAY':
        data_category = 'BYTEARRAY'
    else:
        data_category = 'COMPLEXDATA'

    return restrictions, data_category


def scc_calibrate_mtx(data):
    [ny, nc, nx] = np.shape(data)
    data = np.moveaxis(data,1,-1)
    data = data.flatten().reshape((-1, nc))
    U, s, V = np.linalg.svd(data, full_matrices=False)
    return V[np.newaxis, :, :], s


def gcc_calibrate_mtx(data):
    [ny, nc, nx] = np.shape(data)
    data = np.moveaxis(data,1,-1)
    im = np.fft.ifft(data, axis=1)
    mtx = np.zeros((nx, nc, nc),dtype = 'complex64')
    s = np.zeros(nc, dtype='float32')
    for x in range(nx):
        U, s_, V = np.linalg.svd(im[:,x,:], full_matrices=False)
        mtx[x,:,:] = V
        s += s_
    return mtx, s


def calibrate_mtx(data, cc_mode):
    if cc_mode == 'scc':
        mtx, s = scc_calibrate_mtx(data)
    elif cc_mode == 'gcc':
        mtx, s = gcc_calibrate_mtx(data)
    else:
        print('unknown cc_mode "%s"'%(cc_mode))
        raise ValueError
    return mtx, s


def get_cal_data(meas, cc_mode, remove_os, restriction_categories):
    cal_list, cal_isima, cal_nx, cal_nc = list(), list(), list(), list()
    for mdb_key, mdb in enumerate(meas['mdb']):
        restrictions, data_category =  get_restrictions(mdb.get_flags(), restriction_categories)
        if restrictions=='NO_RESTRICTIONS': # and not mdb.is_flag_set('RAWDATACORRECTION'): # wip!?
            cal_list.append(int(mdb_key))
            cal_isima.append(mdh_def.is_image_scan(mdb.mdh))
            cal_nx.append(int(mdb.data.shape[-1]))
            cal_nc.append(int(mdb.data.shape[0]))
    
    mtx = None
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
        max_calib_samples = 1024
        # first pick image scans
        if cal_list.size<=max_calib_samples:
            # pick all
            pick_mask = np.ones(cal_list.shape, dtype=bool)
        else:
            # pick all image scans
            pick_mask = cal_isima.copy()
            missing_picks = max_calib_samples - pick_mask.sum()
            if missing_picks>0:
                # suplement with non-image scans
                tmp = np.full(pick_mask.size()-pick_mask.sum(), False)
                tmp[:missing_picks] = True
                np.random.shuffle(tmp)
                pick_mask[~pick_mask] = tmp
            
        cal_list = cal_list[pick_mask].tolist()
        if remove_os:
            cal_data = np.zeros((len(cal_list), mode_c, mode_x//2), dtype=np.complex64)
        else:
            cal_data = np.zeros((len(cal_list), mode_c, mode_x), dtype=np.complex64)
        for cnt, mdb_idx in enumerate(cal_list):
            data = meas['mdb'][mdb_idx].data
            if cc_mode=='gcc' and meas['mdb'][mdb_idx].is_flag_set('REFLECT'):
                data = data[:,::-1]
            if remove_os:
                data = np.fft.ifft(data)
                nx = data.shape[-1]
                data = np.delete(data, slice(nx//4, nx*3//4), -1)
                data = np.fft.fft(data)
                data = np.complex64(data)
            cal_data[cnt, :, :] = data

    return cal_data


def compress_twix(infile, outfile, remove_os=False, cc_mode=False, ncc=None, cc_tol=0.05, zfp=False, zfp_tol=1e-5, zfp_prec=None, rm_fidnav=False):
    
    twix = twixtools.read_twix(infile)
    restriction_categories = create_restriction_categories(cc_mode)
    if cc_mode:
        mtx = None
        # calibrate coil compression based on last scan in list (image scan)
        # use the calibration coil weights for all data that fits
        cal_data = get_cal_data(twix[-1], cc_mode, remove_os, restriction_categories)
        mtx, s = calibrate_mtx(cal_data, cc_mode)
        del(cal_data)
        if ncc is None:
            ncc = 1 + np.argwhere(np.cumsum(s)/s.sum() > (1-cc_tol))[0,0]
        mtx = mtx[:,:ncc, :]
        print('coil compression from %d channels to %d virtual channels'%(mtx.shape[-1], ncc))
    else:
        mtx = None
    
    t_start = time.time()
    with h5py.File(outfile, "w") as f:
        
        f.attrs["original_filename"] = os.path.basename(infile)
        f.attrs['remove_os'] = remove_os
        f.attrs['cc_mode'] = cc_mode
        f.attrs['zfp'] = zfp
        if zfp_tol is None:
            f.attrs['zfp_tol'] = -1
        else:
            f.attrs['zfp_tol'] = zfp_tol
        if zfp_prec is None:
            f.attrs['zfp_prec'] = -1
        else:
            f.attrs['zfp_prec'] = zfp_prec
        f.create_dataset("multi_header", data=np.frombuffer(twix[0].tobytes(), 'S1'), compression="gzip", compression_opts=9)
        
        if mtx is not None:
            # save mtx for coil compression
            if zfp:
                f.attrs['mtx_shape'] = mtx.shape
                f.create_dataset("mtx", data=pyzfp.compress(mtx.flatten().view('float32'), tolerance=zfp_tol, precision=zfp_prec, parallel=True))
            else:
                f.create_dataset("mtx", data=mtx, compression="gzip", compression_opts=9)

        for meas_key, meas in enumerate(twix[1:]):
            grp = f.create_group("scan%d"%(meas_key))
            grp.create_dataset("hdr_str", data=meas['hdr_str'], compression="gzip", compression_opts=9)

            # first gather some information about the measurement
            data_len = {'BYTEARRAY': 0, 'COMPLEXDATA': 0}
            mdh_count = len(meas['mdb'])
            for mdb in meas['mdb']:
                if mdb.is_flag_set('ACQEND') or mdb.is_flag_set('SYNCDATA'):
                    data_len['BYTEARRAY']+=1
            data_len['COMPLEXDATA'] = mdh_count - data_len['BYTEARRAY']

            grp.create_dataset('mdh_info', shape=[mdh_count], dtype=mdh_def.scan_hdr_type, compression="gzip", compression_opts=9)
            # not all mdh's have coil_info's (namely ACQEND & SYNCDATA) but for simplicity allocate space anyway (only a few bytes overhead)
            grp.create_dataset('coil_info', shape=[data_len['COMPLEXDATA']], dtype=mdh_def.channel_hdr_type, compression="gzip", compression_opts=9)
            # similarly, just allocate the maximum number of possible coils (64 - in special cases 128 - increase further!?)
            grp.create_dataset('coil_list', shape=[data_len['COMPLEXDATA'], 64], dtype=np.uint8, compression="gzip", compression_opts=9)
            # create list to track for which mdbs os removal is active
            grp.create_dataset('rm_os_active', data=np.zeros(data_len['COMPLEXDATA'], dtype=bool), compression="gzip", compression_opts=9)
            # create list to track which mdbs have been coil compressed
            grp.create_dataset('cc_active', data=np.zeros(data_len['COMPLEXDATA'], dtype=bool), compression="gzip", compression_opts=9)

            dt = h5py.vlen_dtype(np.dtype('S1'))
            # dt = h5py.vlen_dtype(np.dtype('uint8'))
            grp.create_dataset('BYTEARRAY', shape=[data_len['BYTEARRAY']], dtype=dt, compression="gzip", compression_opts=9)

            if zfp:
                dt = h5py.vlen_dtype(np.dtype('uint8'))
                grp.create_dataset('COMPLEXDATA', shape=[data_len['COMPLEXDATA']], dtype=dt)
            else:
                dt = h5py.vlen_dtype(np.dtype('complex64'))
                grp.create_dataset('COMPLEXDATA', shape=[data_len['COMPLEXDATA']], dtype=dt, compression="gzip", compression_opts=9)

            data_count = {'BYTEARRAY': 0, 'COMPLEXDATA': 0}
            for mdb_key, mdb in enumerate(meas['mdb']):
                if rm_fidnav and mdb.is_flag_set('noname60'):
                    # ignore fidnav scans
                    continue

                grp['mdh_info'][mdb_key] = mdb.mdh

                if mdb.is_flag_set('ACQEND') or mdb.is_flag_set('SYNCDATA'): # is_bytearray
                    grp['BYTEARRAY'][data_count['BYTEARRAY']] = mdb.data
                    data_count['BYTEARRAY'] += 1
                else:
                    cur_count = data_count['COMPLEXDATA']
                    restrictions, _ = get_restrictions(mdb.get_flags(), restriction_categories)
                    if restrictions=='NO_COILCOMP':
                        data, grp['rm_os_active'][cur_count], grp['cc_active'][cur_count] = reduce_data(mdb.data, mdb.mdh, remove_os, cc_mode=False)
                    else:
                        data, grp['rm_os_active'][cur_count], grp['cc_active'][cur_count] = reduce_data(mdb.data, mdb.mdh, remove_os, cc_mode=cc_mode, mtx=mtx)

                    data = data.flatten()
                    
                    if zfp:
                        data = pyzfp.compress(data.view('float32'), tolerance=zfp_tol, precision=zfp_prec, parallel=True)
                    
                    grp['COMPLEXDATA'][cur_count] = data
                    
                        # grp['COMPLEXDATA'][pos['COMPLEXDATA']:pos['COMPLEXDATA']+data.size] = data.flatten()
                        # pos['COMPLEXDATA'] += data.size
                    if len(mdb.channel_hdr) > 0:
                        grp['coil_info'][cur_count] = mdb.channel_hdr[0]
                        for coil_key, coil_item in enumerate(mdb.channel_hdr):
                            grp['coil_list'][cur_count, coil_key] = coil_item['ulChannelId']
                    
                    data_count['COMPLEXDATA'] += 1

    elapsed_time = (time.time() - t_start)
    print("compression finished in %d:%02d:%02d h"%(elapsed_time//3600, (elapsed_time%3600)//60, elapsed_time%60))
    print("compression factor = %.2f"%(os.path.getsize(args.infile)/os.path.getsize(args.outfile)))


def reconstruct_twix(infile, outfile=None):
    #wip: function takes no parameters, all necessary information needs to be included in hdf file
    def write_sync_bytes(f):
        syncbytes = (512-(f.tell())%512)%512
        f.write(b'\x00' * syncbytes)

    if outfile is None:
        with h5py.File(infile, "r") as f:
            outfile = f.attrs["original_filename"]
    
    t_start = time.time()

    with h5py.File(infile, "r") as f, open(outfile, 'xb') as fout:
        
        remove_os = f.attrs['remove_os']
        cc_mode = f.attrs['cc_mode']
        zfp = f.attrs['zfp']
        zfp_tol = f.attrs['zfp_tol']
        if zfp_tol<0:
            zfp_tol = None
        zfp_prec = f.attrs['zfp_prec']
        if zfp_prec<0:
            zfp_prec = None
        
        inv_mtx = None
        if cc_mode is not False and 'mtx' in f.keys():
            mtx = f['mtx'][()]
            if zfp:
                mtx_shape = f.attrs['mtx_shape']
                mtx = pyzfp.decompress(mtx, [2*mtx_shape.prod()], np.dtype('float32'), tolerance=zfp_tol, precision=zfp_prec).view('complex64').reshape(mtx_shape)

            inv_mtx = np.zeros_like(mtx).swapaxes(1,-1)
            for x in range(mtx.shape[0]):
                inv_mtx[x,:,:] = np.linalg.pinv(mtx[x,:,:])
            del(mtx)

        restriction_categories = create_restriction_categories(cc_mode)
        
        # allocate space for multi-header
        fout.write(b'\x00' * 10240)
        
        scan_pos = list()
        scan_len = list()
        
        scanlist = [key for key in f.keys() if "scan" in key]
        for key, scan in enumerate(scanlist):
            # keep track of byte pos
            scan_pos.append(fout.tell())

            # get compression info
            rm_os_active = f[scan]['rm_os_active'][()]
            cc_active = f[scan]['cc_active'][()]

            # write header
            f[scan]['hdr_str'][()].tofile(fout)

            data_count = {'BYTEARRAY': 0, 'COMPLEXDATA': 0}
                
            for raw_mdh in f[scan]['mdh_info']:
                mdh = np.frombuffer(raw_mdh, mdh_def.scan_hdr_type)[0]
                n_sampl = mdh['ushSamplesInScan']
                n_coil = mdh['ushUsedChannels']
                mdh_flags = mdh_def.get_flags(mdh)
                is_bytearray = mdh_def.is_flag_set(mdh, 'ACQEND') or mdh_def.is_flag_set(mdh, 'SYNCDATA')
                
                # write mdh
                mdh.tofile(fout)
                
                # write data
                if is_bytearray:
                    data = f[scan]['BYTEARRAY'][data_count['BYTEARRAY']]
                    data.tofile(fout)
                    data_count['BYTEARRAY']+=1
                else:
                    cur_count = data_count['COMPLEXDATA']
                    n_data_sampl = n_sampl
                    if rm_os_active[cur_count]:
                        n_data_sampl //= 2
                    n_data_coils = n_coil
                    if cc_mode and cc_active[cur_count]:
                        n_data_coils = inv_mtx.shape[-1]

                    data = f[scan]['COMPLEXDATA'][cur_count]    
                    if zfp:
                        data = pyzfp.decompress(data, [n_data_coils*2*n_data_sampl], np.dtype('float32'), tolerance=zfp_tol, precision=zfp_prec).view('complex64')
                    data = data.reshape(n_data_coils, n_data_sampl)
                    
                    if cc_mode and cc_active[cur_count]:
                        data = expand_data(data, mdh, rm_os_active[cur_count], cc_mode=cc_mode, inv_mtx=inv_mtx)
                    else:
                        data = expand_data(data, mdh, rm_os_active[cur_count], cc_mode=False)

                    data = data.reshape((n_coil, -1))

                    coil_hdr = np.frombuffer(f[scan]['coil_info'][cur_count], mdh_def.channel_hdr_type)[0].copy()
                    coil_idx = f[scan]['coil_list'][cur_count]
                    for c in range(data.shape[0]):
                        #write channel header
                        coil_hdr['ulChannelId'] = coil_idx[c]
                        coil_hdr.tofile(fout)
                        # write data
                        data[c].tofile(fout)
                    
                    data_count['COMPLEXDATA'] += 1

            # update scan_len
            scan_len.append(fout.tell() - scan_pos[-1])
            if mdh_def.is_flag_set(mdh, 'ACQEND'):
                # scan_len fix for acqend
                scan_len[-1] -= 192

            # add sync bytes between scans
            write_sync_bytes(fout)

        # now write preallocated MultiRaidFileHeader
        n_scans = len(scan_pos)
        multi_header = np.frombuffer(f["multi_header"][()], hdr_def.MultiRaidFileHeader)[0]
        # write NScans
        multi_header['hdr']['count_'] = n_scans
        # write scan_pos & scan_len for each scan
        for i, (pos_, len_) in enumerate(zip(scan_pos, scan_len)):
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
        parser.add_argument('-o', '--outfile', '--out', type=TwixFile('w'),
                            help='Output twix .dat file', required=False)
    else:
        parser.add_argument('-i', '--infile', '--in', type=TwixFile('r'),
                            help='Input twix .dat file', required=True)
        parser.add_argument('-o', '--outfile', '--out', type=HDF5File('w'),
                            help='Output HDF5 file', required=True)
        parser.add_argument("--remove_os", action="store_true")
        parser.add_argument("--remove_fidnav", action="store_true")
        
        parser.add_argument('--cc_mode', choices=[False, 'scc', 'gcc'], default=False)

        group_cc = parser.add_mutually_exclusive_group()
        group_cc.add_argument("-n", "--ncc", "-n_compressed_coils", type=int)
        group_cc.add_argument("-t", "--cc_tol", default=0.05, type=float)

        parser.add_argument("--zfp", action="store_true")
        group_zfp = parser.add_mutually_exclusive_group()
        group_zfp.add_argument("--zfp_tol", default=1e-7, type=float)
        group_zfp.add_argument("--zfp_prec", default=None, type=float)


        parser.add_argument("--testmode", action="store_true")

    args = parser.parse_args()
    
    if args.ncc is not None:
        args.cc_tol = None

    if args.zfp_prec is not None:
        args.zfp_tol = None

    if args.decompress:
        reconstruct_twix(args.infile, args.outfile)
    else:
        compress_twix(args.infile, args.outfile, remove_os=args.remove_os, cc_mode=args.cc_mode, ncc=args.ncc, cc_tol=args.cc_tol, zfp=args.zfp, zfp_tol=args.zfp_tol, zfp_prec=args.zfp_prec, rm_fidnav=args.remove_fidnav)
        if args.testmode:
            with h5py.File(args.outfile, "r") as f:
                print(f.keys())
                print(f['scan0'].keys())

            import tempfile
            tmp_name = tempfile.mktemp(suffix='.dat')
            reconstruct_twix(args.outfile, tmp_name)
            inputsz = os.path.getsize(args.infile)
            comprsz = os.path.getsize(args.outfile)
            reconsz = os.path.getsize(tmp_name)
            print('original size = ', inputsz, ' compressed size =', comprsz, ' comp. factor =', inputsz/comprsz, ' reconstructed size =', reconsz)
            if md5(args.infile) == md5(tmp_name):
                print('md5 hashes match, perfect reconstruction (lossless compression)')
            # test read:
            try:
                twix = twixtools.read_twix(tmp_name)
                print('\nreconstructed twix file successfully parsed')
            except:
                print('\nerror parsing reconstructed twix file')
            # os.remove(tmp_name)
            print('reconstructed .dat file:', tmp_name)
