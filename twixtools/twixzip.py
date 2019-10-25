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

def reduce_data(data, mdh, remove_os=False, cc_mode=False, mtx=None):

    if data.dtype == np.dtype("S1"):
        # nothing to do in case of bytearray
        return data, False, False 

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
        data = np.fft.ifft(data)
        nx = data.shape[-1]
        data = np.delete(data, slice(nx//4, nx*3//4), -1)
        data = np.fft.fft(data)
        data = np.complex64(data)
    
    if cc_active:
        if cc_mode == 'scc':
            data = mtx @ data
        elif cc_mode == 'gcc':
            nx = data.shape[-1]
            if nx!=mtx.shape[0]:
                # nx mismatch; deactivate cc mode
                cc_active = False
            else:
                ncc = mtx.shape[1]
                data = np.fft.ifft(data, axis=-1)
                for x in range(nx):
                    data[:ncc,x] = mtx[x,:,:] @ data[:,x]
                data = np.fft.fft(data[:ncc,:], axis=-1)

    if reflect_data:
        return data[:,::-1], rm_os_active, cc_active
    else:
        return data, rm_os_active, cc_active


def expand_data(data, mdh, remove_os=False, cc_mode=False, inv_mtx=None):
      
    if data.dtype == np.dtype("S1"):
        return data # nothing to do in case of bytearray
    
    if remove_os or cc_mode=='gcc':
        reflect_data = bool(mdh['aulEvalInfoMask'][0] & (1 << 24))
        if reflect_data:
            data = data[:,::-1]
    else:
        reflect_data = False
    
    if cc_mode and inv_mtx is not None:
        if cc_mode=='scc':
            try:
                data = inv_mtx @ data
            except:
                print('error during inv_mtx @ data')
                print('mdh flags:', mdh_def.get_active_flags(mdh))
                print('data shape: ', data.shape)
                print('inv_mtx shape: ', inv_mtx.shape)
        else: # 'gcc'
            nx = data.shape[-1]
            if nx!=inv_mtx.shape[0]:
                # nx mismatch; deactivate cc mode
                cc_active = False
            else:
                nc, ncc = inv_mtx.shape[1:]
                data = np.fft.ifft(data, axis=-1)
                # pad missing channels in data with zeros
                data = np.pad(data, [(0, nc-ncc), (0, 0)])
                for x in range(nx):
                    data[:,x] = inv_mtx[x,:,:] @ data[:ncc,x]
                data = np.fft.fft(data, axis=-1)

    if remove_os:
        data = np.fft.ifft(data)
        nx = data.shape[-1]
        data = np.insert(data, nx//2, np.zeros((nx, 1), dtype=data.dtype), -1)
        data = np.fft.fft(data)
        data = np.complex64(data)
        
    if reflect_data:
        data = data[:,::-1]
    else: 
        return data


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
    return V, s


def gcc_calibrate_mtx(data):
    [ny, nc, nx] = np.shape(data)
    data = np.moveaxis(data,1,-1)
    im = np.fft.ifft(data, axis=1)
    mtx = np.zeros((nx, nc, nc),dtype = 'complex64')
    s = np.zeros((nx, nc),dtype = 'float32')
    for x in range(nx):
        U, s_, V = np.linalg.svd(im[:,x,:], full_matrices=False)
        mtx[x,:,:] = V
        s[x, :] = s_
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


def compress_twix(infile, outfile, remove_os=False, cc_mode=False, ncc=None, cc_tol=0.05, zfp=False, zfp_tol=1e-5):
    
    twix = twixtools.read_twix(infile)
    restriction_categories = create_restriction_categories(cc_mode)
    
    
    t_start = time.time()
    with h5py.File(outfile, "w") as f:
        
        f.attrs["original_filename"] = os.path.basename(infile)
        f.attrs['remove_os'] = remove_os
        f.attrs['cc_mode'] = cc_mode
        f.attrs['zfp'] = zfp
        f.attrs['zfp_tol'] = zfp_tol
        
        f.create_dataset("multi_header", data=np.frombuffer(twix[0].tobytes(), 'S1'), compression="gzip", compression_opts=9)
        
        for meas_key, meas in enumerate(twix[1:]):
            grp = f.create_group("scan%d"%(meas_key))
            grp.create_dataset("hdr_str", data=meas['hdr_str'], compression="gzip", compression_opts=9)

            # first gather some information about the measurement
            data_size = {'BYTEARRAY': 0, 'COMPLEXDATA': 0}
            data_count = {'BYTEARRAY': 0, 'COMPLEXDATA': 0}
            mdh_count = 0            
            cal_list, cal_isima, cal_nx, cal_nc = list(), list(), list(), list()
            for mdb_key, mdb in enumerate(meas['mdb']):
                restrictions, data_category =  get_restrictions(mdb.get_flags(), restriction_categories)
                mdh_count += 1
                data_count[data_category] += 1
                if remove_os and not data_category=='BYTEARRAY' and not mdb.is_flag_set('NOISEADJSCAN'):  # we don't want to remove os in noise data
                    data_size[data_category] += mdb.data.size//2
                else:
                    data_size[data_category] += mdb.data.size
                if cc_mode:
                    if restrictions=='NO_RESTRICTIONS': # and not mdb.is_flag_set('RAWDATACORRECTION'): # wip!?
                        cal_list.append(int(mdb_key))
                        cal_isima.append(mdh_def.is_image_scan(mdb.mdh))
                        cal_nx.append(int(mdb.data.shape[-1]))
                        cal_nc.append(int(mdb.data.shape[0]))
            
            mtx = None
            if cc_mode and len(cal_list)>0:
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
                
                mtx, s = calibrate_mtx(cal_data, cc_mode)
                if cc_mode == 'scc':
                    if ncc is None:
                        ncc = 1 + np.argwhere(np.cumsum(s)/s.sum() > (1-cc_tol))[0,0]
                    # reduce mtx
                    mtx = mtx[:ncc, :]
                else:
                    if ncc is None:
                        s = s.sum(axis=0)
                        ncc = 1 + np.argwhere(np.cumsum(s)/s.sum() > (1-cc_tol))[0,0]
                    # reduce mtx
                    mtx = mtx[:,:ncc, :]

                grp.create_dataset('mtx', data=mtx, compression="gzip", compression_opts=9)

            grp.create_dataset('mdh_info', shape=[mdh_count], dtype=mdh_def.scan_hdr_type, compression="gzip", compression_opts=9)
            # not all mdh's have coil_info's (namely ACQEND & SYNCDATA) but for simplicity allocate space anyway (only a few bytes overhead)
            grp.create_dataset('coil_info', shape=[mdh_count], dtype=mdh_def.channel_hdr_type, compression="gzip", compression_opts=9)
            # similarly, just allocate the maximum number of possible coils (64 - in special cases 128 - increase further!?)
            grp.create_dataset('coil_list', shape=[mdh_count, 64], dtype=np.uint8, compression="gzip", compression_opts=9)
            # create list to track for which mdbs os removal is active
            grp.create_dataset('rm_os_active', data=np.zeros(mdh_count, dtype=bool), compression="gzip", compression_opts=9)
            # create list to track which mdbs have been coil compressed
            grp.create_dataset('cc_active', data=np.zeros(mdh_count, dtype=bool), compression="gzip", compression_opts=9)


            pos = {'BYTEARRAY': 0, 'COMPLEXDATA': 0}
            grp.create_dataset('BYTEARRAY', shape=[data_size['BYTEARRAY']], dtype='S1', compression="gzip", compression_opts=9)

            if zfp:
                dt = h5py.vlen_dtype(np.dtype('uint8'))
                grp.create_dataset('COMPLEXDATA', shape=[mdh_count], dtype=dt)
            else:
                grp.create_dataset('COMPLEXDATA', shape=[data_size['COMPLEXDATA']], dtype=np.complex64, compression="gzip", compression_opts=9)


            for mdb_key, mdb in enumerate(meas['mdb']):
                restrictions, _ = get_restrictions(mdb.get_flags(), restriction_categories)

                if restrictions=='NO_COILCOMP':
                    data, grp['rm_os_active'][mdb_key], grp['cc_active'][mdb_key] = reduce_data(mdb.data, mdb.mdh, remove_os, cc_mode=False)
                else:
                    data, grp['rm_os_active'][mdb_key], grp['cc_active'][mdb_key] = reduce_data(mdb.data, mdb.mdh, remove_os, cc_mode=cc_mode, mtx=mtx)
                
                grp['mdh_info'][mdb_key] = mdb.mdh

                if data.dtype == np.dtype("S1"):
                    grp['BYTEARRAY'][pos['BYTEARRAY']:pos['BYTEARRAY']+data.size] = data.flatten()
                    pos['BYTEARRAY'] += data.size
                else:
                    if zfp:
                        grp['COMPLEXDATA'][mdb_key] = pyzfp.compress(data.flatten().view('float32'), tolerance=zfp_tol, parallel=True)
                    else:
                        grp['COMPLEXDATA'][pos['COMPLEXDATA']:pos['COMPLEXDATA']+data.size] = data.flatten()
                        pos['COMPLEXDATA'] += data.size
                    if len(mdb.channel_hdr) > 0:
                        grp['coil_info'][mdb_key] = mdb.channel_hdr[0]
                        for coil_key, coil_item in enumerate(mdb.channel_hdr):
                            grp['coil_list'][mdb_key, coil_key] = coil_item['ulChannelId']

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

            pos = {'BYTEARRAY': 0, 'COMPLEXDATA': 0}

            inv_mtx = None
            if cc_mode is not False and 'mtx' in f[scan].keys():
                if cc_mode=='gcc':
                    mtx = f[scan]['mtx'][()]
                    inv_mtx = np.zeros_like(mtx).swapaxes(1,-1)
                    for x in range(mtx.shape[0]):
                        inv_mtx[x,:,:] = np.linalg.pinv(mtx[x,:,:])
                    del(mtx)
                else:
                    inv_mtx = np.linalg.pinv(f[scan]['mtx'][()])
                
            for acq_no in range(f[scan]['mdh_info'].size):
                mdh = np.frombuffer(f[scan]['mdh_info'][acq_no], mdh_def.scan_hdr_type)[0]
                coil_hdr = np.frombuffer(f[scan]['coil_info'][acq_no], mdh_def.channel_hdr_type)[0].copy()
                coil_idx = f[scan]['coil_list'][acq_no]
                n_sampl = mdh['ushSamplesInScan']
                n_coil = mdh['ushUsedChannels']
                mdh_flags = mdh_def.get_flags(mdh)
                _, data_category = get_restrictions(mdh_flags, restriction_categories)
                
                # write mdh
                mdh.tofile(fout)
                
                # write data
                if data_category == 'BYTEARRAY':
                    dma_len = np.uint32(mdh['ulFlagsAndDMALength'] % (2**25)) - 192
                    data = f[scan]['BYTEARRAY'][pos['BYTEARRAY']:pos['BYTEARRAY']+dma_len]
                    pos['BYTEARRAY']+=dma_len
                    data.tofile(fout)
                else:
                    n_data_sampl = n_sampl
                    if rm_os_active[acq_no]:
                        n_data_sampl //= 2
                    n_data_coils = n_coil
                    if cc_mode and cc_active[acq_no]:
                        n_data_coils = inv_mtx.shape[-1]
                        
                    if zfp:
                        data = pyzfp.decompress(f[scan]['COMPLEXDATA'][acq_no], [n_data_coils*2*n_data_sampl], np.dtype('float32'), tolerance=zfp_tol).view('complex64').reshape(n_data_coils, n_data_sampl)
                    else:
                        len_ = int(n_data_coils * n_data_sampl)
                        data = f[scan]['COMPLEXDATA'][pos['COMPLEXDATA']:pos['COMPLEXDATA']+len_].reshape(n_data_coils, n_data_sampl)
                        pos['COMPLEXDATA']+=len_

                    if cc_mode and cc_active[acq_no]:
                        data = expand_data(data, mdh, rm_os_active[acq_no], cc_mode=cc_mode, inv_mtx=inv_mtx)
                    else:
                        data = expand_data(data, mdh, rm_os_active[acq_no], cc_mode=False)

                    data = data.reshape((n_coil, -1))
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
        
        parser.add_argument('--cc_mode', choices=[False, 'scc', 'gcc'], default=False)

        group_cc = parser.add_mutually_exclusive_group()
        group_cc.add_argument("-n", "--ncc", "-n_compressed_coils", type=int)
        group_cc.add_argument("-t", "--cc_tol", default=0.05, type=float)

        parser.add_argument("--zfp", action="store_true")
        parser.add_argument("--zfp_tol", default=1e-7, type=float)


        parser.add_argument("--testmode", action="store_true")

    args = parser.parse_args()
    
    if args.ncc is not None:
        args.cc_tol = None

    print(args)

    if args.decompress:
        reconstruct_twix(args.infile, args.outfile)
    else:
        compress_twix(args.infile, args.outfile, remove_os=args.remove_os, cc_mode=args.cc_mode, ncc=args.ncc, cc_tol=args.cc_tol, zfp=args.zfp, zfp_tol=args.zfp_tol)
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
