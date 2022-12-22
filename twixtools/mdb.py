import numpy as np
import os
import ctypes
import copy

import twixtools.mdh_def as mdh_def


class Mdb_base(object):
    def __init__(self, version_is_ve=True):
        self.version_is_ve = version_is_ve
        if self.version_is_ve:
            self.mdh = mdh_def.Scan_header()
            self.channel_hdr = list()
        else:
            self.mdh = mdh_def.VB17_header()
            self.channel_hdr = None

    def __str__(self):
        """Convert to string, for str()."""
        return ( f"{self.__class__.__module__}.{self.__class__.__qualname__}:\n"
                f"  ScanCounter: {self.mdh.ScanCounter}\n"
                f"  TimeStamp / PMUTimeStamp: {self.mdh.TimeStamp} / {self.mdh.PMUTimeStamp}\n"
                f"  Active Flags: {self.get_active_flags()}\n"
                f"  UsedChannels: {self.mdh.UsedChannels}\n"
                f"  SamplesInScan: {self.mdh.SamplesInScan}\n"
                f"  CutOff: {self.mdh.CutOff}\n"
                f"  ReadOutOffcentre: {self.mdh.ReadOutOffcentre:.5f}\n"
                f"  CenterCol: {self.mdh.CenterCol}\n"
                f"  CenterLin: {self.mdh.CenterLin}\n"
                f"  CenterPar: {self.mdh.CenterPar}\n"
                 "  Counter:\n"
                f"    Lin: {self.mdh.Counter.Lin}\n"
                f"    Par: {self.mdh.Counter.Par}\n"
                f"    Sli: {self.mdh.Counter.Sli}\n"
                f"    Ave: {self.mdh.Counter.Ave}\n"
                f"    Rep: {self.mdh.Counter.Rep}\n"
                f"    Set: {self.mdh.Counter.Set}\n"
                f"    Seg: {self.mdh.Counter.Seg}\n"
                f"    Ida: {self.mdh.Counter.Ida}\n"
                f"    Idb: {self.mdh.Counter.Idb}\n"
                f"    Idc: {self.mdh.Counter.Idc}\n"
                f"    Idd: {self.mdh.Counter.Idd}\n"
                f"    Ide: {self.mdh.Counter.Ide}\n"
                f"  SliceData:\n"
                f"    SlicePos:   {self.mdh.SliceData.SlicePos}\n"
                f"    Quaternion: {self.mdh.SliceData.Quaternion}"
        )

    def _get_data_len(self):
        return np.uint32(8 * self.mdh.SamplesInScan)

    def _get_block_len(self):
        if self.version_is_ve:
            return np.uint32(self.data_len + ctypes.sizeof(mdh_def.Channel_header))
        else:
            return np.uint32(self.data_len + ctypes.sizeof(mdh_def.VB17_header))

    def _get_dma_len(self):
        if self.is_flag_set('ACQEND') or self.is_flag_set('SYNCDATA'):
            return mdh_def.get_dma_len(self.mdh)

        # override value found in 'FlagsAndDMALength' which is sometimes
        # not quite correct (e.g. in case PackBit is set)
        out = self.mdh.UsedChannels * self.block_len
        if self.version_is_ve:
            out += ctypes.sizeof(mdh_def.Scan_header)
        return np.uint32(out)

    data_len = property(_get_data_len)
    block_len = property(_get_block_len)
    dma_len = property(_get_dma_len)

    def is_flag_set(self, flag):
        return mdh_def.is_flag_set(self.mdh, flag)

    def get_flag(self, flag):
        return mdh_def.get_flag(self.mdh, flag)

    def set_flag(self, flag, val):
        mdh_def.set_flag(self.mdh, flag, val)

    def add_flag(self, flag):
        mdh_def.add_flag(self.mdh, flag)

    def remove_flag(self, flag):
        mdh_def.remove_flag(self.mdh, flag)

    def get_flags(self):
        return mdh_def.get_flags(self.mdh)

    def get_active_flags(self):
        return mdh_def.get_active_flags(self.mdh)

    def set_flags(self, flags):
        mdh_def.set_flags(self.mdh, flags)

    def clear_all_flags(self):
        mdh_def.clear_all_flags(self.mdh)

    def is_image_scan(self):
        return mdh_def.is_image_scan(self.mdh)

    @property
    def cLin(self):
        return self.mdh.Counter.Lin

    @property
    def cAve(self):
        return self.mdh.Counter.Ave

    @property
    def cSlc(self):
        return self.mdh.Counter.Sli

    @property
    def cPar(self):
        return self.mdh.Counter.Par

    @property
    def cEco(self):
        return self.mdh.Counter.Eco

    @property
    def cPhs(self):
        return self.mdh.Counter.Phs

    @property
    def cRep(self):
        return self.mdh.Counter.Rep

    @property
    def cSet(self):
        return self.mdh.Counter.Set

    @property
    def cSeg(self):
        return self.mdh.Counter.Seg

    def update_CRC(self):
        if self.version_is_ve:
            # CRC is currently not used by Siemens
            self.mdh.CRC = 0
            self.channel_hdr.CRC = 0

    def set_timestamps(self, value=None):
        self.set_timestamp(value)
        self.set_pmutimestamp(value)

    def set_timestamp(self, value=None):
        if value is None:
            from time import time
            value = int(time+0.5)
        self.mdh.TimeStamp = value
        if self.version_is_ve:
            for hdr in self.channel_hdr:
                hdr.SequenceTime = value

    def set_pmutimestamp(self, value=None):
        if value is None:
            from time import time
            value = int(time+0.5)
        self.mdh.PMUTimeStamp = value

    def write_to_file(self, fid):
        pass


class Mdb_local(Mdb_base):

    def __init__(self, data=None, version_is_ve=True):
        super().__init__(version_is_ve=version_is_ve)
        if self.version_is_ve:
            self.mdh = mdh_def.Scan_header()
        else:
            self.mdh = mdh_def.VB17_header()
        if data is None:
            self.__data = (b'\x00' * 160)
            mdh_def.set_dma_len(self.mdh, 160)
        else:
            self._set_data(data)

    def _get_data(self):
        return self.__data

    def _set_data(self, value):
        if isinstance(value, bytes):
            self.__data = value
            mdh_def.set_dma_len(self.mdh, self.dma_len)
            return    
        if value.ndim > 2:
            raise ValueError
        self.__data = np.complex64(np.atleast_2d(value))
        ncha, ncol = self.__data.shape[-2:]
        self._update_hdr(ncha, ncol)

    data = property(_get_data, _set_data)

    def _update_hdr(self, ncha, ncol):
        self.mdh.SamplesInScan = ncol
        self.mdh.UsedChannels = ncha
        mdh_def.set_dma_len(self.mdh, self.dma_len)
        if self.version_is_ve:
            if len(self.channel_hdr) < ncha:
                if len(self.channel_hdr) == 0:
                    self.channel_hdr.append(mdh_def.Channel_header())
                for c in range(ncha):
                    if c >= len(self.channel_hdr):
                        self.channel_hdr.append(self.channel_hdr[0])
                    self.channel_hdr[c].ChannelId = c
            del(self.channel_hdr[ncha:])


class Mdb(Mdb_base):
    """Memory-mapped storage class for Siemens MRI raw data.

    The function `read_twix` parses the twix file and stores all
    information (such as counters as well as the position of the data
    within the twix file) in a list of Mdb objects.

    Important Attributes
    ----------
    data: This will read the data from disk and return a 2D numpy
          array (channels x columns).
    mdh: measurement data header, stored as a special numpy.dtype.
         You can get information about the order of the stored
         information (counters and such) with `print(mdh.dtype)`.

    Important Methods
    ----------
    is_image_scan(self): check for an imaging scan.

    is_flag_set(self, flag): check if a certain MDH flag is set.

    set_flag(self, flag, val): set MDH flag to value `val` (True or False).

    get_flags(self): returns a dict with bools for all MDH flags.

    get_active_flags(self): returns a list of all active MDH flags.

    convert_to_local(self): converts the Mdb to a Mdb_local object
        (no memory-mapping) that allows for array manipulation.
        This makes it possible to to manipulate data before writing
        a new twix file using `write_twix`
    """

    def __init__(self, fid=None, version_is_ve=True):
        super().__init__(version_is_ve=version_is_ve)
        self.mem_pos = None
        if fid is not None:
            self.__read_mdh(fid)

    def convert_to_local(self):
        out = Mdb_local(self.data, self.version_is_ve)
        out.mdh = copy.deepcopy(self.mdh)
        if self.version_is_ve:
            out.channel_hdr = copy.deepcopy(self.channel_hdr)
        return out

    def __read_mdh(self, fid):
        self.fid = fid
        if self.fid.closed:
            return ValueError
        self.mem_pos = self.fid.tell()
        if self.version_is_ve:
            self.mdh = mdh_def.Scan_header()
            self.fid.readinto(self.mdh)
            if self.is_flag_set('ACQEND') or self.is_flag_set('SYNCDATA'):
                return

            bytes_inter = ctypes.sizeof(mdh_def.Channel_header)
            block_sz = bytes_inter + self.data_len
            buffer = self.fid.read(self.mdh.UsedChannels * block_sz)
            self.channel_hdr = [
                mdh_def.Channel_header.from_buffer_copy(
                    buffer, c*block_sz) for c in range(self.mdh.UsedChannels)]
        else:
            self.mdh = mdh_def.VB17_header()
            self.fid.readinto(self.mdh)

    def __get_data(self):
        was_closed = self.fid.closed
        if was_closed:
            # opening/closing files adds a huge overhead, try to avoid
            self.fid = open(self.fid.name, self.fid.mode)

        if self.version_is_ve:
            bytes_initial = ctypes.sizeof(mdh_def.Scan_header)
            bytes_inter = ctypes.sizeof(mdh_def.Channel_header)
        else:
            bytes_initial = 0
            bytes_inter = ctypes.sizeof(mdh_def.VB17_header)

        self.fid.seek(self.mem_pos + bytes_initial)
        out = self.fid.read(self.dma_len - bytes_initial)

        if self.is_flag_set('ACQEND') or self.is_flag_set('SYNCDATA'):
            pass  # nothing to do here
        else:
            dt = np.dtype(
                [('skip', bytes, bytes_inter),
                 ('data', np.complex64, self.mdh.SamplesInScan)])
            out = np.frombuffer(
                out, dtype=dt, count=self.mdh.UsedChannels)['data']

        if was_closed:
            self.fid.close()

        return out

    data = property(__get_data)
