import ctypes


class SeqDataHeader(ctypes.LittleEndianStructure):
    _pack_ = 1
    _fields_ = [
        ("packet_size", ctypes.c_uint32),
        ("id", ctypes.c_char * 52),
        ("swapped", ctypes.c_uint32)
    ]


class SeqData():
    def __init__(self, data):
        self.hdr = SeqDataHeader.from_buffer_copy(data)
        self.data = data[ctypes.sizeof(self.hdr):ctypes.sizeof(self.hdr) + self.hdr.packet_size]
        self.appendix = data[ctypes.sizeof(self.hdr) + self.hdr.packet_size:]  # no idea what this is/means

    def __repr__(self):
        return f"SeqData(packet_size={self.hdr.packet_size}, id={self.hdr.id.decode('utf-8')[:50]}, swapped={self.hdr.swapped})"

    def __len__(self):
       return self.hdr.packet_size + ctypes.sizeof(SeqDataHeader) + len(self.appendix)

    def tobytes(self):
        return bytes(self.hdr) + self.data + self.appendix
