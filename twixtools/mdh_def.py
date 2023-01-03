import numpy as np
import ctypes


class MyStruct(ctypes.LittleEndianStructure):

    def __eq__(self, other):
        for field in self._fields_:
            attr_name = field[0]
            a, b = object.__getattribute__(
                self, attr_name), object.__getattribute__(other, attr_name)
            is_array = isinstance(a, ctypes.Array)
            if is_array and a[:] != b[:] or not is_array and a != b:
                return False
        return True

    def __ne__(self, other):
        for field in self._fields_:
            attr_name = field[0]
            a, b = object.__getattribute__(
                self, attr_name), object.__getattribute__(other, attr_name)
            is_array = isinstance(a, ctypes.Array)
            if is_array and a[:] != b[:] or not is_array and a != b:
                return True
        return False

    def __getattribute__(self, name):
        val = object.__getattribute__(self, name)
        if isinstance(val, ctypes.Array):
            return np.ctypeslib.as_array(val)
        else:
            return val

    def __str__(self):
        out = list()
        for field in self._fields_:
            name = field[0]
            val = getattr(self, name)
            if isinstance(val, ctypes.Structure):
                # add one more recursion depth:
                lst = []
                for item in val._fields_:
                    val2 = getattr(val, item[0])
                    if isinstance(val2, ctypes.Structure):
                        lst.append(str(val2))
                    else:
                        lst.append(val2)
                out.append(lst)
            else:
                out.append(val)
        return str(out)


class LineCounter(MyStruct):
    _pack_ = 1
    _fields_ = [
        ("Lin", ctypes.c_uint16),
        ("Ave", ctypes.c_uint16),
        ("Sli", ctypes.c_uint16),
        ("Par", ctypes.c_uint16),
        ("Eco", ctypes.c_uint16),
        ("Phs", ctypes.c_uint16),
        ("Rep", ctypes.c_uint16),
        ("Set", ctypes.c_uint16),
        ("Seg", ctypes.c_uint16),
        ("Ida", ctypes.c_uint16),
        ("Idb", ctypes.c_uint16),
        ("Idc", ctypes.c_uint16),
        ("Idd", ctypes.c_uint16),
        ("Ide", ctypes.c_uint16)]


class CutOff(MyStruct):
    _pack_ = 1
    _fields_ = [
        ("Pre", ctypes.c_uint16),
        ("Post", ctypes.c_uint16)]


class SlicePos(MyStruct):
    _pack_ = 1
    _fields_ = [
        ("Sag", ctypes.c_float),
        ("Cor", ctypes.c_float),
        ("Tra", ctypes.c_float)]


class SliceData(MyStruct):
    _pack_ = 1
    _fields_ = [
        ("SlicePos", SlicePos),
        ("Quaternion", ctypes.c_float * 4)]


# This is the VB line header
class VB17_header(MyStruct):
    _pack_ = 1
    _fields_ = [
        ("FlagsAndDMALength", ctypes.c_uint32),
        ("MeasUID", ctypes.c_int32),
        ("ScanCounter", ctypes.c_uint32),
        ("TimeStamp", ctypes.c_uint32),
        ("PMUTimeStamp", ctypes.c_uint32),
        ("EvalInfoMask", ctypes.c_uint64),
        ("SamplesInScan", ctypes.c_uint16),
        ("UsedChannels", ctypes.c_uint16),
        ("Counter", LineCounter),
        ("CutOff", CutOff),
        ("CenterCol", ctypes.c_uint16),
        ("CoilSelect", ctypes.c_uint16),
        ("ReadOutOffcentre", ctypes.c_float),
        ("TimeSinceLastRF", ctypes.c_uint32),
        ("CenterLin", ctypes.c_uint16),
        ("CenterPar", ctypes.c_uint16),
        ("IceProgramPara", ctypes.c_uint16*4),
        ("FreePara", ctypes.c_uint16*4),
        ("SliceData", SliceData),
        ("ChannelId", ctypes.c_uint16),
        ("PTABPosNeg", ctypes.c_uint16)]

    def __init__(self, **kwargs):
        # initialize some fields with sane defaults
        values = {"FlagsAndDMALength": ctypes.sizeof(self), "ScanCounter": 1}
        values.update(kwargs)
        super().__init__(**values)


# VD/VE: One scan header for all channels
class Scan_header(MyStruct):
    _pack_ = 1
    _fields_ = [
        ("FlagsAndDMALength", ctypes.c_uint32),
        ("MeasUID", ctypes.c_int32),
        ("ScanCounter", ctypes.c_uint32),
        ("TimeStamp", ctypes.c_uint32),
        ("PMUTimeStamp", ctypes.c_uint32),
        ("SystemType", ctypes.c_uint16),
        ("PTABPosDelay", ctypes.c_uint16),
        ("PTABPosX", ctypes.c_int32),
        ("PTABPosY", ctypes.c_int32),
        ("PTABPosZ", ctypes.c_int32),
        ("Reserved1", ctypes.c_int32),
        ("EvalInfoMask", ctypes.c_uint64),
        ("SamplesInScan", ctypes.c_uint16),
        ("UsedChannels", ctypes.c_uint16),
        ("Counter", LineCounter),
        ("CutOff", CutOff),
        ("CenterCol", ctypes.c_uint16),
        ("CoilSelect", ctypes.c_uint16),
        ("ReadOutOffcentre", ctypes.c_float),
        ("TimeSinceLastRF", ctypes.c_uint32),
        ("CenterLin", ctypes.c_uint16),
        ("CenterPar", ctypes.c_uint16),
        ("SliceData", SliceData),
        ("IceProgramPara", ctypes.c_uint16*24),
        ("ReservedPara", ctypes.c_uint16*4),
        ("ApplicationCounter", ctypes.c_uint16),
        ("ApplicationMask", ctypes.c_uint16),
        ("CRC", ctypes.c_uint32)]

    def __init__(self, **kwargs):
        # initialize some fields with sane defaults
        values = {"FlagsAndDMALength": ctypes.sizeof(self), "ScanCounter": 1}
        values.update(kwargs)
        super().__init__(**values)


# VD/VE: One channel header per channel
class Channel_header(MyStruct):
    _pack_ = 1
    _fields_ = [
        ("TypeAndChannelLength", ctypes.c_uint32),
        ("MeasUID", ctypes.c_int32),
        ("ScanCounter", ctypes.c_uint32),
        ("Reserved1", ctypes.c_int32),
        ("SequenceTime", ctypes.c_uint32),
        ("Unused2", ctypes.c_uint32),
        ("ChannelId", ctypes.c_uint16),
        ("Unused3", ctypes.c_uint16),
        ("CRC", ctypes.c_uint32)]

    # initialize some fields with sane defaults
    def __init__(self, **kwargs):
        values = {"ScanCounter": 1}
        values.update(kwargs)
        super().__init__(**values)


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
    'noname43',
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
        infomask, 2**np.arange(64, dtype=np.uint64)).astype(bool)


def get_dma_len(mdh):
    return np.uint32(mdh.FlagsAndDMALength % (2**25))


def set_dma_len(mdh, dma_len):
    split = 2**25
    if dma_len > split-1:
        raise ValueError
    mdh.FlagsAndDMALength = mdh.FlagsAndDMALength//split + dma_len


def is_flag_set(mdh, flag):
    return bool(int(mdh.EvalInfoMask) & 1 << mask_dict[flag])


def get_flag(mdh, flag):
    return is_flag_set(mdh, flag)


def set_flag(mdh, flag, val):
    if val:
        add_flag(mdh, flag)
    else:
        remove_flag(mdh, flag)


def add_flag(mdh, flag):
    mdh.EvalInfoMask |= (1 << mask_dict[flag])


def remove_flag(mdh, flag):
    mdh.EvalInfoMask &= ~(1 << mask_dict[flag])


def get_flags(mdh):
    mask = unpack_bits(mdh.EvalInfoMask)
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
    mdh.EvalInfoMask = 0


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
