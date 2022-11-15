import numpy as np


mdhLC = [("Lin", "<u2"),
         ("Ave", "<u2"),
         ("Sli", "<u2"),
         ("Par", "<u2"),
         ("Eco", "<u2"),
         ("Phs", "<u2"),
         ("Rep", "<u2"),
         ("Set", "<u2"),
         ("Seg", "<u2"),
         ("Ida", "<u2"),
         ("Idb", "<u2"),
         ("Idc", "<u2"),
         ("Idd", "<u2"),
         ("Ide", "<u2")]

mdhCutOff = [("Pre", "<u2"),
             ("Post", "<u2")]

mdhSlicePosVec = [("Sag", "<f4"),
                  ("Cor", "<f4"),
                  ("Tra", "<f4")]

mdhSliceData = [("SlicePos", mdhSlicePosVec),
                ("Quaternion", "<f4", 4)]


# This is the VB line header
vb17_header = [("FlagsAndDMALength", "<u4"),
               ("MeasUID", "<i4"),
               ("ScanCounter", "<u4"),
               ("TimeStamp", "<u4"),
               ("PMUTimeStamp", "<u4"),
               ("EvalInfoMask", "<u8"),
               ("SamplesInScan", "<u2"),
               ("UsedChannels", "<u2"),
               ("Counter", mdhLC),
               ("CutOff", mdhCutOff),
               ("CenterCol", "<u2"),
               ("CoilSelect", "<u2"),
               ("ReadOutOffcentre", "<f4"),
               ("TimeSinceLastRF", "<u4"),
               ("CenterLin", "<u2"),
               ("CenterPar", "<u2"),
               ("IceProgramPara", "<u2", 4),
               ("FreePara", "<u2", 4),
               ("SliceData", mdhSliceData),
               ("ChannelId", "<u2"),
               ("PTABPosNeg", "<u2")]


# VD/VE: One scan header for all channels
scan_header = [("FlagsAndDMALength", "<u4"),
               ("MeasUID", "<i4"),
               ("ScanCounter", "<u4"),
               ("TimeStamp", "<u4"),
               ("PMUTimeStamp", "<u4"),
               ("SystemType", "<u2"),
               ("PTABPosDelay", "<u2"),
               ("PTABPosX", "<i4"),
               ("PTABPosY", "<i4"),
               ("PTABPosZ", "<i4"),
               ("Reserved1", "<i4"),
               ("EvalInfoMask", "<u8"),
               ("SamplesInScan", "<u2"),
               ("UsedChannels", "<u2"),
               ("Counter", mdhLC),
               ("CutOff", mdhCutOff),
               ("CenterCol", "<u2"),
               ("CoilSelect", "<u2"),
               ("ReadOutOffcentre", "<f4"),
               ("TimeSinceLastRF", "<u4"),
               ("CenterLin", "<u2"),
               ("CenterPar", "<u2"),
               ("SliceData", mdhSliceData),
               ("IceProgramPara", "<u2", 24),
               ("ReservedPara", "<u2", 4),
               ("ApplicationCounter", "<u2"),
               ("ApplicationMask", "<u2"),
               ("CRC", "<u4")]


# VD/VE: One channel header per channel
channel_header = [("TypeAndChannelLength", "<u4"),
                  ("MeasUID", "<i4"),
                  ("ScanCounter", "<u4"),
                  ("Reserved1", "<i4"),
                  ("SequenceTime", "<u4"),
                  ("Unused2", "<u4"),
                  ("ChannelId", "<u2"),
                  ("Unused3", "<u2"),
                  ("CRC", "<u4")]

vb17_hdr_type = np.dtype(vb17_header)
scan_hdr_type = np.dtype(scan_header)
channel_hdr_type = np.dtype(channel_header)


mask_id = (
    'ACQEND',  # last scan
    'RTFEEDBACK',  # Realtime feedback scan
    'HPFEEDBACK',  # High perfomance feedback scan
    'ONLINE',  # processing should be done online
    'OFFLINE',  # processing should be done offline
    'SYNCDATA',  # readout contains synchroneous data
    'noname6',
    'noname7',
    'LASTSCANINCONCAT',  # Flag for last scan in concatination
    'noname9',
    'RAWDATACORRECTION',  # Correct with the rawdata corr. factor
    'LASTSCANINMEAS',  # Flag for last scan in measurement
    'SCANSCALEFACTOR',  # Flag for scan specific additional scale
    '2NDHADAMARPULSE',  # 2nd RF exitation of HADAMAR
    'REFPHASESTABSCAN',  # reference phase stabilization scan
    'PHASESTABSCAN',  # phase stabilization scan
    'D3FFT',  # execute 3D FFT
    'SIGNREV',  # sign reversal
    'PHASEFFT',  # execute phase fft
    'SWAPPED',  # swapped phase/readout direction
    'POSTSHAREDLINE',  # shared line
    'PHASCOR',  # phase correction data
    'PATREFSCAN',  # additional scan for PAT ref line/partition
    'PATREFANDIMASCAN',  # PAT ref that is also used as image scan
    'REFLECT',  # reflect line
    'NOISEADJSCAN',  # noise adjust scan
    'SHARENOW',  # lines may be shared between e.g. phases
    'LASTMEASUREDLINE',  # indicates last meas line of e.g. phases
    'FIRSTSCANINSLICE',  # first scan in slice; req for timestamps
    'LASTSCANINSLICE',  # last scan in slice; req for timestamps
    'TREFFECTIVEBEGIN',  # indicates the TReff begin (triggered)
    'TREFFECTIVEEND',  # indicates the TReff end (triggered)
    'REF_POSITION',  # indicates ref pos for move during scan acq.
    'SLC_AVERAGED',  # indicates averaged slice for sl. partial av
    'TAGFLAG1',  # adjust scans
    'CT_NORMALIZE',  # Marks corr maps scan for TimCT-Prescan norm
    'SCAN_FIRST',  # Marks the first scan of a particular map
    'SCAN_LAST',  # Marks the last scan of a particular map
    'SLICE_ACCEL_REFSCAN',  # single-band ref. scan for multi-band
    'SLICE_ACCEL_PHASCOR',  # additional phase corr. in multi-band
    'FIRST_SCAN_IN_BLADE',  # Marks the first line of a blade
    'LAST_SCAN_IN_BLADE',  # Marks the last line of a blade
    'LAST_BLADE_IN_TR',  # Marks all lin. of last BLADE in each TR
    'PACE',  # Distinguishes PACE scans from non PACE scans.
    'RETRO_LASTPHASE',  # Marks the last phase in a heartbeat
    'RETRO_ENDOFMEAS',  # Marks an ADC at end of measurement
    'RETRO_REPEATTHISHEARTBEAT',  # Repeat the current heartbeat
    'RETRO_REPEATPREVHEARTBEAT',  # Repeat the previous heartbeat
    'RETRO_ABORTSCANNOW',  # Just abort everything
    'RETRO_LASTHEARTBEAT',  # adc is from last heartbeat (a dummy)
    'RETRO_DUMMYSCAN',  # adc is just a dummy scan, throw it away
    'RETRO_ARRDETDISABLED',  # Disable all arrhythmia detection
    'B1_CONTROLLOOP',  # readout to be used for B1 Control Loop
    'SKIP_ONLINE_PHASCOR',  # scans not to be online phase corr.
    'SKIP_REGRIDDING',  # Marks scans not to be regridded
    'MDH_VOP',  # Marks scans to be used for VOP based RF monitoring
    'noname57',
    'noname58',
    'noname59',
    'noname60',
    'WIP_1',  # Mark scans for WIP application "type 1"
    'WIP_2',  # Mark scans for WIP application "type 1"
    'WIP_3'   # Mark scans for WIP application "type 1"
)

# create dict for faster access by name
mask_dict = {item: key for key, item in enumerate(mask_id)}


# helper function (copied from mdb class)

def unpack_bits(infomask):
    # numpy's unpackbits does not work correctly for some reason
    return np.bitwise_and(
        infomask, np.left_shift(2, np.arange(-1, 8*infomask.nbytes-1, dtype=np.uint64))).astype(bool)


def is_flag_set(mdh, flag):
    return bool(int(mdh['EvalInfoMask']) & 1 << mask_dict[flag])


def get_flag(mdh, flag):
    return is_flag_set(mdh, flag)


def set_flag(mdh, flag, val):
    if val:
        add_flag(mdh, flag)
    else:
        remove_flag(mdh, flag)


def add_flag(mdh, flag):
    mdh['EvalInfoMask'] |= np.uint64(1 << mask_dict[flag])


def remove_flag(mdh, flag):
    mdh['EvalInfoMask'] &= ~np.uint64(1 << mask_dict[flag])


def get_flags(mdh):
    mask = unpack_bits(mdh['EvalInfoMask'])
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
    mdh['EvalInfoMask'] = 0


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
