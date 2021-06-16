import numpy as np


mdhLC = [("ushLine", "<u2"),
         ("ushAcquisition", "<u2"),
         ("ushSlice", "<u2"),
         ("ushPartition", "<u2"),
         ("ushEcho", "<u2"),
         ("ushPhase", "<u2"),
         ("ushRepetition", "<u2"),
         ("ushSet", "<u2"),
         ("ushSeg", "<u2"),
         ("ushIda", "<u2"),
         ("ushIdb", "<u2"),
         ("ushIdc", "<u2"),
         ("ushIdd", "<u2"),
         ("ushIde", "<u2")]

mdhCutOff = [("ushPre", "<u2"),
             ("ushPost", "<u2")]

mdhSlicePosVec = [("flSag", "<f4"),
                  ("flCor", "<f4"),
                  ("flTra", "<f4")]

mdhSliceData = [("sSlicePosVec", mdhSlicePosVec),
                ("aflQuaternion", "<f4", 4)]


# This is the VB line header
vb17_header = [("ulFlagsAndDMALength", "<u4"),
               ("lMeasUID", "<i4"),
               ("ulScanCounter", "<u4"),
               ("ulTimeStamp", "<u4"),
               ("ulPMUTimeStamp", "<u4"),
               ("aulEvalInfoMask", "<u4", 2),
               ("ushSamplesInScan", "<u2"),
               ("ushUsedChannels", "<u2"),
               ("sLC", mdhLC),
               ("sCutOff", mdhCutOff),
               ("ushKSpaceCentreColumn", "<u2"),
               ("ushCoilSelect", "<u2"),
               ("fReadOutOffcentre", "<f4"),
               ("ulTimeSinceLastRF", "<u4"),
               ("ushKSpaceCentreLineNo", "<u2"),
               ("ushKSpaceCentrePartitionNo", "<u2"),
               ("aushIceProgramPara", "<u2", 4),
               ("aushFreePara", "<u2", 4),
               ("sSliceData", mdhSliceData),
               ("ushChannelId", "<u2"),
               ("ushPTABPosNeg", "<u2")]


# VD/VE: One scan header for all channels
scan_header = [("ulFlagsAndDMALength", "<u4"),
               ("lMeasUID", "<i4"),
               ("ulScanCounter", "<u4"),
               ("ulTimeStamp", "<u4"),
               ("ulPMUTimeStamp", "<u4"),
               ("ushSystemType", "<u2"),
               ("ulPTABPosDelay", "<u2"),
               ("lPTABPosX", "<i4"),
               ("lPTABPosY", "<i4"),
               ("lPTABPosZ", "<i4"),
               ("ulReserved1", "<i4"),
               ("aulEvalInfoMask", "<u4", 2),
               ("ushSamplesInScan", "<u2"),
               ("ushUsedChannels", "<u2"),
               ("sLC", mdhLC),
               ("sCutOff", mdhCutOff),
               ("ushKSpaceCentreColumn", "<u2"),
               ("ushCoilSelect", "<u2"),
               ("fReadOutOffcentre", "<f4"),
               ("ulTimeSinceLastRF", "<u4"),
               ("ushKSpaceCentreLineNo", "<u2"),
               ("ushKSpaceCentrePartitionNo", "<u2"),
               ("sSliceData", mdhSliceData),
               ("aushIceProgramPara", "<u2", 24),
               ("aushReservedPara", "<u2", 4),
               ("ushApplicationCounter", "<u2"),
               ("ushApplicationMask", "<u2"),
               ("ulCRC", "<u4")]


# VD/VE: One channel header per channel
channel_header = [("ulTypeAndChannelLength", "<u4"),
                  ("lMeasUID", "<i4"),
                  ("ulScanCounter", "<u4"),
                  ("ulReserved1", "<i4"),
                  ("ulSequenceTime", "<u4"),
                  ("ulUnused2", "<u4"),
                  ("ulChannelId", "<u2"),
                  ("ulUnused3", "<u2"),
                  ("ulCRC", "<u4")]

vb17_hdr_type = np.dtype(vb17_header)
scan_hdr_type = np.dtype(scan_header)
channel_hdr_type = np.dtype(channel_header)


mask_id = ['noname%d' % (k) for k in range(64)]
mask_id[0] = 'ACQEND'  # last scan
mask_id[1] = 'RTFEEDBACK'  # Realtime feedback scan
mask_id[2] = 'HPFEEDBACK'  # High perfomance feedback scan
mask_id[3] = 'ONLINE'  # processing should be done online
mask_id[4] = 'OFFLINE'  # processing should be done offline
mask_id[5] = 'SYNCDATA'  # readout contains synchroneous data
mask_id[8] = 'LASTSCANINCONCAT'  # Flag for last scan in concatination
mask_id[10] = 'RAWDATACORRECTION'  # Correct with the rawdata corr. factor
mask_id[11] = 'LASTSCANINMEAS'  # Flag for last scan in measurement
mask_id[12] = 'SCANSCALEFACTOR'  # Flag for scan specific additional scale
mask_id[13] = '2NDHADAMARPULSE'  # 2nd RF exitation of HADAMAR
mask_id[14] = 'REFPHASESTABSCAN'  # reference phase stabilization scan
mask_id[15] = 'PHASESTABSCAN'  # phase stabilization scan
mask_id[16] = 'D3FFT'  # execute 3D FFT
mask_id[17] = 'SIGNREV'  # sign reversal
mask_id[18] = 'PHASEFFT'  # execute phase fft
mask_id[19] = 'SWAPPED'  # swapped phase/readout direction
mask_id[20] = 'POSTSHAREDLINE'  # shared line
mask_id[21] = 'PHASCOR'  # phase correction data
mask_id[22] = 'PATREFSCAN'  # additional scan for PAT ref line/partition
mask_id[23] = 'PATREFANDIMASCAN'  # PAT ref that is also used as image scan
mask_id[24] = 'REFLECT'  # reflect line
mask_id[25] = 'NOISEADJSCAN'  # noise adjust scan
mask_id[26] = 'SHARENOW'  # lines may be shared between e.g. phases
mask_id[27] = 'LASTMEASUREDLINE'  # indicates last meas line of e.g. phases
mask_id[28] = 'FIRSTSCANINSLICE'  # first scan in slice; req for timestamps
mask_id[29] = 'LASTSCANINSLICE'  # last scan in slice; req for timestamps
mask_id[30] = 'TREFFECTIVEBEGIN'  # indicates the TReff begin (triggered)
mask_id[31] = 'TREFFECTIVEEND'  # indicates the TReff end (triggered)
mask_id[32] = 'REF_POSITION'  # indicates ref pos for move during scan acq.
mask_id[33] = 'SLC_AVERAGED'  # indicates averaged slice for sl. partial av
mask_id[34] = 'TAGFLAG1'  # adjust scans
mask_id[35] = 'CT_NORMALIZE'  # Marks corr maps scan for TimCT-Prescan norm
mask_id[36] = 'SCAN_FIRST'  # Marks the first scan of a particular map
mask_id[37] = 'SCAN_LAST'  # Marks the last scan of a particular map
mask_id[38] = 'SLICE_ACCEL_REFSCAN'  # single-band ref. scan for multi-band
mask_id[39] = 'SLICE_ACCEL_PHASCOR'  # additional phase corr. in multi-band
mask_id[40] = 'FIRST_SCAN_IN_BLADE'  # Marks the first line of a blade
mask_id[41] = 'LAST_SCAN_IN_BLADE'  # Marks the last line of a blade
mask_id[42] = 'LAST_BLADE_IN_TR'  # Marks all lin. of last BLADE in each TR
mask_id[44] = 'PACE'  # Distinguishes PACE scans from non PACE scans.
mask_id[45] = 'RETRO_LASTPHASE'  # Marks the last phase in a heartbeat
mask_id[46] = 'RETRO_ENDOFMEAS'  # Marks an ADC at end of measurement
mask_id[47] = 'RETRO_REPEATTHISHEARTBEAT'  # Repeat the current heartbeat
mask_id[48] = 'RETRO_REPEATPREVHEARTBEAT'  # Repeat the previous heartbeat
mask_id[49] = 'RETRO_ABORTSCANNOW'  # Just abort everything
mask_id[50] = 'RETRO_LASTHEARTBEAT'  # adc is from last heartbeat (a dummy)
mask_id[51] = 'RETRO_DUMMYSCAN'  # adc is just a dummy scan, throw it away
mask_id[52] = 'RETRO_ARRDETDISABLED'  # Disable all arrhythmia detection
mask_id[53] = 'B1_CONTROLLOOP'  # readout to be used for B1 Control Loop
mask_id[54] = 'SKIP_ONLINE_PHASCOR'  # scans not to be online phase corr.
mask_id[55] = 'SKIP_REGRIDDING'  # Marks scans not to be regridded
mask_id[56] = 'MDH_VOP'  # Marks scans to be used for VOP based RF monitoring
mask_id[61] = 'WIP_1'  # Mark scans for WIP application "type 1"
mask_id[62] = 'WIP_2'  # Mark scans for WIP application "type 1"
mask_id[63] = 'WIP_3'  # Mark scans for WIP application "type 1"


# helper function (copied from mdb class)

def unpack_bits(infomask):
    # numpy's unpackbits does not work correctly for some reason
    infomask = infomask.view(np.uint64)[0]
    return np.bitwise_and(
        infomask, 2**np.arange(8*infomask.nbytes)).astype(bool)


def is_flag_set(mdh, flag):
    mask = mdh['aulEvalInfoMask']
    bit = mask_id.index(flag)
    if bit < 32:
        return bool(mask[0] & (1 << bit))
    else:
        return bool(mask[1] & (1 << bit-32))


def get_flag(mdh, flag):
    return is_flag_set(mdh, flag)


def set_flag(mdh, flag, val):
    if val:
        add_flag(mdh, flag)
    else:
        remove_flag(mdh, flag)


def add_flag(mdh, flag):
    bit = mask_id.index(flag)
    if bit < 32:
        mdh['aulEvalInfoMask'][0] |= (1 << bit)
    else:
        mdh['aulEvalInfoMask'][1] |= (1 << bit-32)


def remove_flag(mdh, flag):
    bit = mask_id.index(flag)
    if bit < 32:
        mdh['aulEvalInfoMask'][0] &= ~(1 << bit)
    else:
        mdh['aulEvalInfoMask'][1] &= ~(1 << bit-32)


def get_flags(mdh):
    mask = unpack_bits(mdh['aulEvalInfoMask'])
    return dict(zip(mask_id, mask))


def get_active_flags(mdh):
    return [key for key, item in get_flags(mdh).items() if item]


def set_flags(mdh, flags):
    if isinstance(flags, list):
        for key in flags:
            set_flag(mdh, key, True)
    elif isinstance(flags, dict):
        for key, item in flags:
            set_flag(mdh, key, item)
    else:
        raise ValueError


def clear_all_flags(mdh):
    mdh['aulEvalInfoMask'][0] = 0
    mdh['aulEvalInfoMask'][1] = 0


def is_image_scan(mdh):
    disqualifier = [
        'ACQEND', 'RTFEEDBACK', 'HPFEEDBACK', 'SYNCDATA', 'REFPHASESTABSCAN',
        'PHASESTABSCAN', 'PHASCOR', 'NOISEADJSCAN', 'noname60']
    for name in disqualifier:
        if is_flag_set(mdh, name):
            return False
    # check for patref scan
    if is_flag_set(mdh, 'PATREFSCAN')\
            and not is_flag_set(mdh, 'PATREFANDIMASCAN'):
        return False
    return True
