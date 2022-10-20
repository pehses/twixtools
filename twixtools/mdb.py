import numpy as np
import os

import twixtools.mdh_def as mdh_def


class Mdb_base(object):
    def __init__(self, version_is_ve=True):
        self.version_is_ve = version_is_ve
        if self.version_is_ve:
            self.mdh = np.zeros(1, dtype=mdh_def.scan_hdr_type)[0]
            self.channel_hdr = list()
        else:
            self.mdh = np.zeros(1, dtype=mdh_def.vb17_hdr_type)
            self.channel_hdr = None

    def __str__(self):
        return ( "Mdb:\n"
                f"  ScanCounter: {self.mdh['ScanCounter']}\n"
                f"  TimeStamp / PMUTimeStamp: {self.mdh['TimeStamp']} / {self.mdh['PMUTimeStamp']}\n"
                f"  Active Flags: {self.get_active_flags()}\n"
                f"  UsedChannels: {self.mdh['UsedChannels']}\n"
                f"  SamplesInScan: {self.mdh['SamplesInScan']}\n"
                f"  CutOff: {self.mdh['CutOff']}\n"
                f"  ReadOutOffcentre: {self.mdh['ReadOutOffcentre']:.5f}\n"
                f"  CenterCol: {self.mdh['CenterCol']}\n"
                f"  CenterLin: {self.mdh['CenterLin']}\n"
                f"  CenterPar: {self.mdh['CenterPar']}\n"
                 "  Counter:\n"
                f"    Lin: {self.mdh['Counter']['Lin']}\n"
                f"    Par: {self.mdh['Counter']['Par']}\n"
                f"    Sli: {self.mdh['Counter']['Sli']}\n"
                f"    Ave: {self.mdh['Counter']['Ave']}\n"
                f"    Rep: {self.mdh['Counter']['Rep']}\n"
                f"    Set: {self.mdh['Counter']['Set']}\n"
                f"    Seg: {self.mdh['Counter']['Seg']}\n"
                f"    Ida->Ide: {list(self.mdh['Counter'])[-5:]}\n"
                f"  SliceData:\n"
                f"    SlicePos:   {self.mdh['SliceData']['SlicePos']}\n"
                f"    Quaternion: {self.mdh['SliceData']['Quaternion']}"
        )

    def _get_data_len(self):
        return np.uint32(8 * self.mdh['SamplesInScan'])

    def _get_block_len(self):
        if self.version_is_ve:
            return np.uint32(self.data_len + mdh_def.channel_hdr_type.itemsize)
        else:
            return np.uint32(self.data_len + mdh_def.vb17_hdr_type.itemsize)

    def _get_dma_len(self):
        if self.is_flag_set('ACQEND') or self.is_flag_set('SYNCDATA'):
            return np.uint32(self.mdh['FlagsAndDMALength'] % (2**25))

        # override value found in 'FlagsAndDMALength' which is sometimes
        # not quite correct (e.g. in case PackBit is set)
        out = self.mdh['UsedChannels'] * self.block_len
        if self.version_is_ve:
            out += mdh_def.scan_hdr_type.itemsize
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
        return self.mdh['Counter']['Lin']

    @property
    def cAve(self):
        return self.mdh['Counter']['Ave']

    @property
    def cSlc(self):
        return self.mdh['Counter']['Sli']

    @property
    def cPar(self):
        return self.mdh['Counter']['Par']

    @property
    def cEco(self):
        return self.mdh['Counter']['Eco']

    @property
    def cPhs(self):
        return self.mdh['Counter']['Phs']

    @property
    def cRep(self):
        return self.mdh['Counter']['Rep']

    @property
    def cSet(self):
        return self.mdh['Counter']['Set']

    @property
    def cSeg(self):
        return self.mdh['Counter']['Seg']

    def update_CRC(self):
        if self.version_is_ve:
            # CRC is currently not used by Siemens
            self.mdh["CRC"] = 0
            self.channel_hdr["CRC"] = 0

    def set_timestamps(self, value=None):
        self.set_timestamp(value)
        self.set_pmutimestamp(value)

    def set_timestamp(self, value=None):
        if value is None:
            from time import time
            value = int(time+0.5)
        self.mdh["TimeStamp"] = value
        if self.version_is_ve:
            for hdr in self.channel_hdr:
                hdr["SequenceTime"] = value

    def set_pmutimestamp(self, value=None):
        if value is None:
            from time import time
            value = int(time+0.5)
        self.mdh["PMUTimeStamp"] = value

    def write_to_file(self, fid):
        pass


class Mdb_local(Mdb_base):

    def __init__(self, data=None, version_is_ve=True):
        super().__init__(version_is_ve=version_is_ve)
        if data is not None:
            self._set_data(data)

    def _get_data(self):
        return self.__data

    def _set_data(self, value):
        if value.ndim > 2:
            raise ValueError
        self.__data = np.complex64(np.atleast_2d(value))
        # ncha, ncol = self.__data.shape[-2:]
        # self._update_hdr(ncha, ncol)

    data = property(_get_data, _set_data)

    def _update_hdr(self, ncha, ncol):
        if self.mdh is None:
            self.mdh = dict()
        self.mdh["SamplesInScan"] = ncol
        self.mdh["UsedChannels"] = ncha
        if not self.is_flag_set('ACQEND') and not self.is_flag_set('SYNCDATA'):
            pass
        if self.version_is_ve:
            if self.channel_hdr is None or len(self.channel_hdr) < ncha:
                self.channel_hdr = [np.zeros(1,
                                             dtype=mdh_def.channel_hdr_type)]
            for c in range(ncha):
                if c >= len(self.channel_hdr):
                    self.channel_hdr.append(self.channel_hdr[0])
                self.channel_hdr[c]["ChannelId"] = c
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
        out.mdh = self.mdh.copy()
        out.channel_hdr = self.channel_hdr.copy()
        return out

    def __read_mdh(self, fid):
        self.fid = fid
        if self.fid.closed:
            return ValueError
        self.mem_pos = self.fid.tell()
        if self.version_is_ve:
            self.mdh = np.fromfile(
                self.fid, dtype=mdh_def.scan_hdr_type, count=1)[0]
            if not self.is_flag_set('ACQEND')\
                    and not self.is_flag_set('SYNCDATA'):
                for _ in range(self.mdh['UsedChannels']):
                    chan_hd = np.fromfile(
                        self.fid, dtype=mdh_def.channel_hdr_type, count=1)[0]
                    self.channel_hdr.append(chan_hd)
                    self.fid.seek(self.data_len, os.SEEK_CUR)
        else:
            self.mdh = np.fromfile(
                self.fid, dtype=mdh_def.vb17_hdr_type, count=1)[0]

    def __get_data(self):
        was_closed = self.fid.closed
        if was_closed:
            # opening/closing files adds a huge overhead, try to avoid
            self.fid = open(self.fid.name, self.fid.mode)
        self.fid.seek(self.mem_pos)
        if self.version_is_ve:
            self.fid.seek(mdh_def.scan_hdr_type.itemsize, os.SEEK_CUR)
            skip_bytes = mdh_def.channel_hdr_type.itemsize
        else:
            skip_bytes = mdh_def.vb17_hdr_type.itemsize

        if self.is_flag_set('ACQEND') or self.is_flag_set('SYNCDATA'):
            # channel header is in this case assumed to be part of 'data'
            dma_len_ = np.uint32(self.mdh['FlagsAndDMALength'] % (2**25))
            if not self.version_is_ve:
                self.fid.seek(mdh_def.vb17_hdr_type.itemsize, os.SEEK_CUR)
            dma_len_ -= mdh_def.scan_hdr_type.itemsize
            out = np.fromfile(self.fid, dtype='<S1', count=dma_len_)
        else:
            dt = np.dtype(
                [('skip', bytes, skip_bytes),
                 ('data', np.complex64, self.mdh['SamplesInScan'])])
            out = np.fromfile(
                self.fid, dtype=dt, count=self.mdh['UsedChannels'])['data']

        if was_closed:
            self.fid.close()
        return out

    data = property(__get_data)
