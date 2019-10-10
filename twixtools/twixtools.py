from __future__ import print_function   # for python 2.7 compatibility
import os
import re
import numpy as np

from twixtools.twixprot import twixprot, parse_twix_hdr
from twixtools.helpers import idea_version_check, update_progress
from twixtools.mdb import Mdb, Mdb_base
from twixtools.hdr_def import MultiRaidFileHeader, SingleMeasInit


def read_twix(infile, read_prot=True, keep_syncdata_and_acqend=True):
    """Function for reading siemens twix raw data files."""
    if isinstance(infile, str):
        # assume that complete path is given
        if infile[-4:].lower() != '.dat':
            infile += '.dat'   # adds filetype ending to file
    else:
        # filename not a string, so assume that it is the MeasID
        measID = infile
        infile = [f for f in os.listdir('.') if re.search(
            r'^meas_MID0*' + str(measID) + '.*\.dat$', f)]
        if len(infile) == 0:
            print('error: .dat file with measID', measID, 'not found')
            raise ValueError
        elif len(infile) > 1:
            print('multiple files with measID', measID,
                  'found, choosing first occurence')
        infile = infile[0]

    infile = os.path.realpath(infile)

    fid = open(infile, 'rb')
    fid.seek(0, os.SEEK_END)
    fileSize = np.uint64(fid.tell())
    version_is_ve, NScans = idea_version_check(fid)

    out = list()
    # lazy software version check (VB or VD?)
    if version_is_ve:
        print('Software version: VD/VE (!?)')
        fid.seek(0, os.SEEK_SET)  # move pos to 9th byte in file
        raidfile_hdr = np.fromfile(fid, dtype=MultiRaidFileHeader, count=1)[0]
        out.append(raidfile_hdr)
        NScans = raidfile_hdr["hdr"]["count_"]
        measOffset = list()
        measLength = list()
        for k in range(NScans):
            measOffset.append(raidfile_hdr['entry'][k]['off_'])
            measLength.append(raidfile_hdr['entry'][k]['len_'])
    else:
        # in VB versions, the first 4 bytes indicate the beginning of the
        # raw data part of the file
        print('Software  : VB (!?)')
        # VB does not support multiple scans in one file:
        measOffset = [np.uint64(0)]
        measLength = [fileSize]

    for s in range(NScans):
        scanStart = measOffset[s]
        scanEnd = scanStart + measLength[s]
        pos = measOffset[s]
        fid.seek(pos, os.SEEK_SET)
        meas_init = np.fromfile(fid, dtype=SingleMeasInit, count=1)[0]
        hdr_len = meas_init["hdr_len"]
        out.append(dict())
        if read_prot:
            fid.seek(pos, os.SEEK_SET)
            hdr = parse_twix_hdr(fid)
            out[-1]['init'] = meas_init
            out[-1]['hdr'] = hdr
            fid.seek(pos, os.SEEK_SET)
            out[-1]['hdr_str'] = np.fromfile(fid, dtype="<S1", count=hdr_len)
            # prot = twixprot(fid, hdr_len, version_is_ve)
            # out[-1]['prot'] = prot
        out[-1]['mdb'] = list()

        pos = measOffset[s] + np.uint64(hdr_len)
        scanStart = pos
        print('\nscan ', s)
        update_progress(pos - scanStart, scanEnd - scanStart, True)
        while pos + 128 < scanEnd:  # fail-safe not to miss ACQEND
            update_progress(pos - scanStart, scanEnd - scanStart, False)
            fid.seek(pos, os.SEEK_SET)
            mdb = Mdb(fid, version_is_ve)

            # jump to mdh of next scan
            pos += mdb.dma_len

            if not keep_syncdata_and_acqend:
                if mdb.is_flag_set('SYNCDATA'):
                    continue
                elif mdb.is_flag_set('ACQEND'):
                    break

            out[-1]['mdb'].append(mdb)

            if mdb.is_flag_set('ACQEND'):
                break
            

    # fid.close()

    return out


do_not_zfp_compress = ['SYNCDATA', 'ACQEND']
do_not_remove_os = ['SYNCDATA', 'ACQEND', '__fidnav___']
do_not_scc_compress = ['SYNCDATA', 'ACQEND', 'RTFEEDBACK']
do_not_gcc_compress = ['SYNCDATA', 'ACQEND', 'RTFEEDBACK', 'NOISEADJSCAN', '__fidnav___']
#'REFLECT'

def write_twix(scanlist, outfile, version_is_ve=True):
    
    def write_sync_bytes(fid):
        syncbytes = (512-(fid.tell())%512)%512
        fid.write(b'\x00' * syncbytes)

    if isinstance(scanlist, dict):
        scanlist = [scanlist]

    with open(outfile, 'xb') as fid:
        if version_is_ve:
            # allocate space for multi-header
            fid.write(b'\x00' * 10240)
        
        scan_pos = list()
        scan_len = list()
        for key, scan in enumerate(scanlist):
            
            if not isinstance(scan, dict):
                continue

            # keep track of byte pos
            scan_pos.append(fid.tell())

            # write header
            scan['hdr_str'].tofile(fid)

            acq_end_len = 0

            # write mdbs
            for mdb in scan['mdb']:
                # write mdh
                mdb.mdh.tofile(fid)
                data = np.atleast_2d(mdb.data)
                if version_is_ve:
                    if mdb.is_flag_set('SYNCDATA'):
                        data.tofile(fid)
                    elif mdb.is_flag_set('ACQEND'):
                        data.tofile(fid)
                        acq_end_len = 192
                    else:
                        for c in range(data.shape[0]):
                            #write channel header
                            mdb.channel_hdr[c].tofile(fid)
                            # write data
                            data[c].tofile(fid)
                else: # WIP: VB version
                    mdb.mdh.tofile(fid)
                    # write data
                    data[c].tofile(fid)

            # update scan_len
            scan_len.append(fid.tell() - scan_pos[-1] - acq_end_len)

            # add sync bytes between scans
            write_sync_bytes(fid)

        # now write preallocated MultiRaidFileHeader
        if version_is_ve:
            n_scans = len(scan_pos)
            if isinstance(scanlist[0], np.void):
                # we have a template
                multi_header = scanlist[0].copy()
            else:
                # start from scratch
                multi_header = np.zeros(1, dtype=MultiRaidFileHeader)[0]
                for k in range(n_scans):
                    multi_header['entry']['patName_'] = b'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
                    multi_header['entry']['protName_'] = b'noname'

            # write NScans
            multi_header['hdr']['count_'] = n_scans
            # write scan_pos & scan_len for each scan
            for i, (pos_, len_) in enumerate(zip(scan_pos, scan_len)):
                multi_header['entry'][i]['len_'] = len_
                multi_header['entry'][i]['off_'] = pos_
                
            # write MultiRaidFileHeader
            fid.seek(0)
            multi_header.tofile(fid)
