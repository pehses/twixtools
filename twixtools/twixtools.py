"""
twixtools: provides reading and limited writing capability of Siemens MRI raw
           data files (.dat).

@author: Philipp Ehses (philipp.ehses@dzne.de)
"""

import os
import re
import numpy as np
from tqdm import tqdm

import twixtools.twixprot as twixprot
import twixtools.helpers as helpers
import twixtools.mdb
import twixtools.hdr_def as hdr_def
import twixtools.geometry


def read_twix(infile, include_scans=None, parse_prot=True, parse_data=True,
              parse_geometry=True, verbose=True, keep_syncdata_and_acqend=False):
    """Function for reading siemens twix raw data files.

    Parameters
    ----------
    infile : filename or measurement id of .dat file
    include_scans: list of scan numbers or None, optional
        By default, all scans in a multi-raid file are parsed.
    parse_prot : bool, optional
        By default, the protocol information is parsed.
        (this is also highly recommended)
    parse_data: bool, optional
        Set to False to parse only protocol information.
    parse_geometry: bool, optional
        Set to False to skip creation of transformation matrix from data
        coordinates to physical coordinates.
    verbose: bool, optional
        Switch progress and other messages on or off.
    keep_syncdata_and_acqend : bool, optional
        By default, syncdata and acqend blocks are not included in the mdb list.
        These blocks are helpful for twix writing, but unnecessary otherwise.

    Returns
    -------
    out: list of twix scan(s)
        The twix scans themselves consist of a dict with these elements:
            - hdr: dict of parsed ascconv and XProtocol header information
            - hdr_str: header bytearray (used by write_twix)
            - mdb: list of measurement data blocks -- here is the MRI data
              use `help(twixtools.mdb.Mdb)` for more information
            - geometry: dict containing geometry information about the scan
              use `help(twixtools.geometry)` for more information
    """
    if isinstance(infile, str):
        # assume that complete path is given
        # check if filepath contains an extension
        if len(os.path.splitext(infile)[1]) == 0:
            infile += '.dat'   # adds filetype ending to file
    else:
        # filename not a string, so assume that it is the MeasID
        measID = infile
        infile = [f for f in os.listdir('.') if re.search(
            r'^meas_MID0*' + str(measID) + r'.*\.dat$', f)]
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
    version_is_ve, NScans = helpers.idea_version_check(fid)

    out = list()
    # lazy software version check (VB or VD?)
    if version_is_ve:
        if verbose:
            print('Software version: VD/VE (!?)')
        fid.seek(0, os.SEEK_SET)  # move pos to 9th byte in file
        raidfile_hdr = np.fromfile(fid, dtype=hdr_def.MultiRaidFileHeader,
                                   count=1)[0]
        # WARNING:
        # it is probably no longer necessary to append the raidfile_hdr for
        # lossless twix writing (as long as there are no changes to the twix
        # format!), so we don't need the following line
        # out.append(raidfile_hdr)
        NScans = raidfile_hdr["hdr"]["count_"]
        measOffset = list()
        measLength = list()
        for k in range(NScans):
            if raidfile_hdr['entry'][k]['len_'] == 0:
                if NScans == 1:
                    raidfile_hdr['entry'][0]['len_'] = fileSize - raidfile_hdr['entry'][0]['off_']
                    print('WARNING: raidfile_hdr of single-raidfile has length 0; determing length from file size')
                else:
                    print(f'WARNING: raidfile_hdr entry {k} has length 0')
            measOffset.append(raidfile_hdr['entry'][k]['off_'])
            measLength.append(raidfile_hdr['entry'][k]['len_'])
    else:
        # in VB versions, the first 4 bytes indicate the beginning of the
        # raw data part of the file
        if verbose:
            print('Software  : VB (!?)')
        # VB does not support multiple scans in one file:
        measOffset = [np.uint64(0)]
        measLength = [fileSize]

    if include_scans is not None:
        # force include_scans to be iterable:
        if not hasattr(include_scans, '__iter__'):
            include_scans = [include_scans]
        # handle negative indexing:
        include_scans = [range(NScans)[k] for k in include_scans]

    if verbose:
        print('')

    for s in range(NScans):
        if include_scans is not None and s not in include_scans:
            # skip scan if it is not requested
            continue

        scanStart = measOffset[s]
        scanEnd = scanStart + measLength[s]
        pos = measOffset[s]
        fid.seek(pos, os.SEEK_SET)
        meas_init = np.fromfile(fid, dtype=hdr_def.SingleMeasInit, count=1)[0]
        hdr_len = meas_init["hdr_len"]
        out.append(dict())
        out[-1]['mdb'] = list()
        if parse_prot:
            fid.seek(pos, os.SEEK_SET)
            hdr = twixprot.parse_twix_hdr(fid)
            out[-1]['hdr'] = hdr
        fid.seek(pos, os.SEEK_SET)
        out[-1]['hdr_str'] = np.fromfile(fid, dtype="<S1", count=hdr_len)

        if version_is_ve:
            out[-1]['raidfile_hdr'] = raidfile_hdr['entry'][s]

        if parse_geometry:
            if not parse_prot:
                print('WARNING: geometry parsing requires protocol parsing, skipping geometry parsing')
                parse_geometry = False
            else:
                out[-1]['geometry'] = twixtools.geometry.Geometry(out[-1])

        # if data is not requested (headers only)
        if not parse_data:
            continue

        pos = measOffset[s] + np.uint64(hdr_len)

        if verbose:
            print('Scan ', s)
            progress_bar = tqdm(total=scanEnd - pos, unit='B', unit_scale=True, unit_divisor=1024)
        while pos + 128 < scanEnd:  # fail-safe not to miss ACQEND
            fid.seek(pos, os.SEEK_SET)
            try:
                mdb = twixtools.mdb.Mdb(fid, version_is_ve)
            except ValueError:
                print(f"WARNING: Mdb parsing encountered an error at file position {pos}/{scanEnd}, stopping here.")

            # jump to mdh of next scan
            pos += mdb.dma_len
            if verbose:
                progress_bar.update(mdb.dma_len)

            if not keep_syncdata_and_acqend:
                if mdb.is_flag_set('SYNCDATA'):
                    continue
                elif mdb.is_flag_set('ACQEND'):
                    break

            out[-1]['mdb'].append(mdb)

            if mdb.is_flag_set('ACQEND'):
                break

        if verbose:
            progress_bar.close()

    fid.close()

    return out


def write_twix(scanlist, outfile, version_is_ve=True):
    """Function for writing siemens twix raw data files.

    Parameters
    ----------
    scanlist: list of twix scan(s)
    outfile: output filename for twix file (.dat)
    version_is_ve: bool that determines what whether to write a VA/VB
        or VD/VE compatible twix file.
        IMPORTANT: This tool does not allow for conversion between versions.
        This bool should be set to the original twix file version!
        IMPORTANT: `write_twix` currently only supports VE twix files!
    """

    def write_sync_bytes(fid):
        syncbytes = (512-(fid.tell()) % 512) % 512
        fid.write(b'\x00' * syncbytes)

    if isinstance(scanlist, dict):
        scanlist = [scanlist]

    with open(outfile, 'xb') as fid:
        if version_is_ve:
            # allocate space for multi-header
            fid.write(b'\x00' * 10240)

        scan_pos = list()
        scan_len = list()
        for scan in scanlist:

            if not isinstance(scan, dict):
                continue

            # keep track of byte pos
            scan_pos.append(fid.tell())

            # write header
            scan['hdr_str'].tofile(fid)

            # make sure that scan counters are consecutive integers
            fix_scancounters(scan['mdb'])

            # write mdbs
            for mdb in scan['mdb']:
                # write mdh
                fid.write(bytearray(mdb.mdh))
                if version_is_ve:
                    if mdb.is_flag_set('SYNCDATA')\
                            or mdb.is_flag_set('ACQEND'):
                        fid.write(mdb.data)
                        if mdb.is_flag_set('ACQEND')\
                                and mdb is not scan['mdb'][-1]:
                            print("WARNING: Early ACQEND detected, skipping some data during write.")
                            break
                    else:
                        data = np.atleast_2d(mdb.data)
                        for c in range(data.shape[0]):
                            # write channel header
                            fid.write(bytearray(mdb.channel_hdr[c]))
                            # write data
                            data[c].tofile(fid)
                else:  # WIP: VB version
                    data = np.atleast_2d(mdb.data)
                    for c in range(data.shape[0]):
                        fid.write(bytearray(mdb.mdh))
                        data[c].tofile(fid)

            if not mdb.is_flag_set('ACQEND'):
                print("ACQEND missing at the end of the mdb list. Generating new one.")
                acqend = twixtools.mdb.Mdb_local()
                acqend.add_flag('ACQEND')
                acqend.mdh.ScanCounter = mdb.mdh.ScanCounter + 1
                acqend.mdh.TimeMeasUID = mdb.mdh.MeasUID
                acqend.mdh.TimeStamp = mdb.mdh.TimeStamp
                acqend.mdh.PMUTimeStamp = mdb.mdh.PMUTimeStamp
                fid.write(bytearray(acqend.mdh))
                fid.write(acqend.data)

            # update scan_len
            scan_len.append(fid.tell() - scan_pos[-1])

            # add sync bytes between scans
            write_sync_bytes(fid)

        # now write preallocated MultiRaidFileHeader
        if version_is_ve:
            multi_header = construct_multiheader(scanlist)
            # write scan_pos & scan_len for each scan (in case they changed)
            for key, (pos_, len_) in enumerate(zip(scan_pos, scan_len)):
                multi_header['entry'][key]['len_'] = len_
                multi_header['entry'][key]['off_'] = pos_
            # write MultiRaidFileHeader
            fid.seek(0)
            multi_header.tofile(fid)


def construct_multiheader(scanlist):
    multi_header = np.zeros(1, dtype=hdr_def.MultiRaidFileHeader)[0]
    n_scans = len(scanlist)
    for key, scan in enumerate(scanlist):
        if 'raidfile_hdr' in scan:
            # we have a template
            multi_header['entry'][key] = scan['raidfile_hdr'].copy()
        else:  # not really necessary anymore, but may be helpful in future
            # start from scratch
            pat = b'x' * 20
            prot = b'noname'
            meas_id = np.uint32(0)
            file_id = np.uint32(0)
            if 'file_id' in scan:
                file_id = scan['file_id']
            if 'hdr' in scan and 'Config' in scan['hdr']:
                config = scan['hdr']['Config']
                if 'tPatientName' in config:
                    pat = config['tPatientName'].encode()
                    if all([item == 'x' for item in pat.decode()]):
                        # fix number of 'x' for anonymized names (lucky? guess)
                        pat = b'x' * min(64, len(pat)+9)
                if 'SequenceDescription' in config:
                    prot = config['SequenceDescription'].encode()
                if 'SequenceDescription' in config:
                    meas_id = np.uint32(config['MeasUID'])

            multi_header['entry']['measId_'][key] = meas_id
            multi_header['entry']['fileId_'][key] = file_id
            multi_header['entry']['patName_'][key] = pat
            multi_header['entry']['protName_'][key] = prot

    # write NScans
    multi_header['hdr']['count_'] = n_scans

    return multi_header


def fix_scancounters(mdb_list, start_cnt=1):
    ''' makes sure that all ScanCounter in mdb_list are consecutive integers
    This is necessary if mdbs are added/removed to/from a mdb_list.
    '''
    cnt = start_cnt
    for mdb in mdb_list:
        if mdb.is_flag_set('SYNCDATA'):  # ignore SYNCDATA
            continue
        mdb.mdh.ScanCounter = cnt
        for cha in mdb.channel_hdr:
            cha.ScanCounter = cnt
        cnt += 1


def del_from_mdb_list(mdb_list, function):
    ''' helper function to safely remove multiple items from mdb_list at once
    Parameters
    ----------
    mdb_list: input list of mdbs
    function: function used to filter mdbs

    Example
    --------
    Remove all mdbs from mdb_list which have the flag 'noname60' set to True.
    >>> del_from_mdb_list(mdb_list, lambda mdb: mdb.is_flag_set('noname60'))
    '''

    ind2remove = [key for key, mdb in enumerate(mdb_list) if function(mdb)]

    for key in sorted(ind2remove, reverse=True):
        del mdb_list[key]

    return
