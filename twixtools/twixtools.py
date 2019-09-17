from __future__ import print_function   # for python 2.7 compatibility
import os
import re
import numpy as np

from twixtools.twixprot import twixprot, parse_twix_hdr
from twixtools.helpers import idea_version_check, update_progress
from twixtools.mdb import Mdb, Mdb_base


def read_twix(infile, read_prot=True, keep_syncdata_and_acqend=False):
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

    # lazy software version check (VB or VD?)
    if version_is_ve:
        print('Software version: VD/VE (!?)')
        fid.seek(8, os.SEEK_SET)  # move pos to 9th byte in file
        measID, fileID = np.fromfile(fid, dtype=np.uint32, count=2)
        measOffset = list()
        measLength = list()
        for _ in range(NScans):
            offset, length = np.fromfile(fid, dtype=np.uint64, count=2)
            measOffset.append(offset)
            measLength.append(length)
            fid.seek(152 - 16, os.SEEK_CUR)
    else:
        # in VB versions, the first 4 bytes indicate the beginning of the
        # raw data part of the file
        print('Software  : VB (!?)')
        # VB does not support multiple scans in one file:
        measOffset = [np.uint64(0)]
        measLength = [fileSize]

    out = list()
    for s in range(NScans):
        scanStart = measOffset[s]
        scanEnd = scanStart + measLength[s]
        pos = measOffset[s]
        fid.seek(pos, os.SEEK_SET)
        hdr_len = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        out.append(dict())
        out[-1]['mdb'] = list()
        if read_prot:
            prot = twixprot(fid, hdr_len, version_is_ve)
            fid.seek(pos, os.SEEK_SET)
            hdr = parse_twix_hdr(fid)
            out[-1]['prot'] = prot
            out[-1]['hdr'] = hdr

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


def write_header(hdr, fid, mdb_bytesize):
    pass


def write_twix(scanlist, outfile):
    if type(scanlist) == dict:
        scanlist = [scanlist]
    with open(outfile, 'xb') as fid:
        for scan in scanlist:
            mdb_bytesize = 0
            for mdb in scan['mdb']:
                if issubclass(type(mdb), Mdb_base):
                    bytesize += 0 #wip
            write_header(scan['hdr'], fid, mdb_bytesize)
            for mdb in scan['mdb']:
                if issubclass(type(mdb), Mdb_base):
                    mdb.write_to_file(fid)
