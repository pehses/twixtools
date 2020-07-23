import numpy as np
import ctypes
import os

import twixtools.mdh_def as mdh_def


def unpack_bits(infomask):
    # numpy's unpackbits does not work correctly for some reason
    infomask = infomask.view(np.uint64)[0]
    return np.bitwise_and(infomask, 2**np.arange(8*infomask.nbytes)).astype(bool)


class Mdb_base(object):
    def __init__(self, version_is_ve=True):
        self.version_is_ve = version_is_ve
        if self.version_is_ve:
            self.mdh = np.zeros(1, dtype=mdh_def.scan_hdr_type)[0]
            self.channel_hdr = list()
        else:
            self.mdh = np.zeros(1, dtype=mdh_def.vb17_hdr_type)
            self.channel_hdr = None

    def _get_data_len(self):
        return np.uint32(8 * self.mdh['ushSamplesInScan'])

    def _get_block_len(self):
        if self.version_is_ve:
            return np.uint32(self.data_len + mdh_def.channel_hdr_type.itemsize)
        else:
            return np.uint32(self.data_len + mdh_def.vb17_hdr_type.itemsize)
    
    def _get_dma_len(self):
        if self.is_flag_set('ACQEND') or self.is_flag_set('SYNCDATA'):
            dma_len = np.uint32(self.mdh['ulFlagsAndDMALength'] % (2**25))
            return dma_len

        # override value found in 'ulFlagsAndDMALength' which is sometimes
        # not quite correct (e.g. in case PackBit is set)
        out = self.mdh['ushUsedChannels'] * self.block_len
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
        return self.mdh['sLC']['ushLine']
    
    @property
    def cAcq(self):
        return self.mdh['sLC']['ushAcquisition']

    @property
    def cSlc(self):
        return self.mdh['sLC']['ushSlice']

    @property
    def cPar(self):
        return self.mdh['sLC']['ushPartition']

    @property
    def cEco(self):
        return self.mdh['sLC']['ushEcho']

    @property
    def cPhs(self):
        return self.mdh['sLC']['ushPhase']

    @property
    def cRep(self):
        return self.mdh['sLC']['ushRepetition']

    @property
    def cSet(self):
        return self.mdh['sLC']['ushSet']

    @property
    def cSeg(self):
        return self.mdh['sLC']['ushSeg']

    def update_CRC(self):
        if version_is_ve:
            # CRC is currently not used by Siemens
            self.mdh["ulCRC"] = 0
            self.channel_header["ulCRC"] = 0

    def set_timestamps(self, value=None):
        self.set_timestamp(value)
        self.set_pmutimestamp(value)

    def set_timestamp(self, value=None):
        if value is None:
            from time import time
            value = int(time+0.5)
        self.mdh["ulTimeStamp"] = value
        if self.version_is_ve:
            for hdr in self.channel_hdr:
                self.channel_hdr["ulSequenceTime"] = value

    def set_pmutimestamp(self, value=None):
        if value is None:
            from time import time
            value = int(time+0.5)
        self.mdh["ulPMUTimeStamp"] = value

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
        ncha, ncol = self.__data.shape
        # self._update_hdr(ncha, ncol)

    data = property(_get_data, _set_data)

    def _update_hdr(self, ncha, ncol):
        if self.mdh is None:
            self.mdh = dict()
        self.mdh["ushSamplesInScan"] = ncol
        self.mdh["ushUsedChannels"] = ncha
        if not self.is_flag_set('ACQEND') and not self.is_flag_set('SYNCDATA'):
            pass
            ##self.mdh["ulFlagsAndDMALength"] # set packbit to zero and dma len to expected len (self.dma_len)
        if self.version_is_ve:
            if self.channel_hdr is None or len(self.channel_hdr)<ncha:
                self.channel_hdr = [np.zeros(1, dtype=mdh_def.channel_hdr_type)]
            for c in range(ncha):
                if c>=len(self.channel_hdr):
                    self.channel_hdr.append(self.channel_hdr[0])
                self.channel_hdr[c]["ulChannelId"] = c
            del(self.channel_hdr[ncha:])


class Mdb(Mdb_base):
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
            self.mdh = np.fromfile(self.fid, dtype=mdh_def.scan_hdr_type, count=1)[0]
            if not self.is_flag_set('ACQEND') and not self.is_flag_set('SYNCDATA'):
                for c in range(self.mdh['ushUsedChannels']):
                    chan_hd = np.fromfile(self.fid, dtype=mdh_def.channel_hdr_type, count=1)[0]
                    self.channel_hdr.append(chan_hd)
                    self.fid.seek(self.data_len, os.SEEK_CUR)
        else:
            self.mdh = np.fromfile(self.fid, dtype=mdh_def.vb17_hdr_type, count=1)[0]

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
            dma_len_ = np.uint32(self.mdh['ulFlagsAndDMALength'] % (2**25))
            if not self.version_is_ve:
                self.fid.seek(mdh_def.vb17_hdr_type.itemsize, os.SEEK_CUR)
            dma_len_ -= mdh_def.scan_hdr_type.itemsize
            out = np.fromfile(self.fid, dtype='<S1', count=dma_len_)
        else:
            dt = np.dtype([('skip', bytes, skip_bytes), ('data', np.complex64, self.mdh['ushSamplesInScan'])])
            out = np.fromfile(self.fid, dtype=dt, count=self.mdh['ushUsedChannels'])['data']

        if was_closed:
            self.fid.close()
        return out

    data = property(__get_data)
