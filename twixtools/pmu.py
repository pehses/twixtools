import ctypes
import itertools
import struct
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np

pmu_magic = {
    "END": 0xFFFFFFFF,
    "ECG1": 0x0,
    "ECG2": 0x1,
    "ECG3": 0x2,
    "ECG4": 0x3,
    "PULS": 0x4,
    "RESP": 0x5,
    "EXT1": 0x6,
    "EXT2": 0x7,
    "EVENT": 0x8,
    "UNKNOWN1": 0x9,
    "UNKNOWN2": 0xA,
}

# Create reversed mapping: magic number to key
magic_pmu = {v: k for k, v in pmu_magic.items()}

Header = namedtuple(
    "Header", ["timestamp0", "timestamp", "packet_no", "duration", "unknown"]
)


class SeqDataHeader(ctypes.LittleEndianStructure):
    _pack_ = 1
    _fields_ = [
        ("packet_size", ctypes.c_uint32),
        ("id", ctypes.c_char * 52),
        ("swapped", ctypes.c_uint32),
    ]


class SeqData:
    def __init__(self, data):
        self.hdr = SeqDataHeader.from_buffer_copy(data)
        self.data = data[
            ctypes.sizeof(self.hdr) : ctypes.sizeof(self.hdr) + self.hdr.packet_size
        ]


class PMUblock:
    HEADER_FORMAT = "<IIIHH"
    HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
    MAGIC_SIZE = struct.calcsize("<I")
    PERIOD_SIZE = struct.calcsize("<I")
    DATA_SIZE = struct.calcsize("<I")

    def __init__(self, data):
        self.header = Header(
            *struct.unpack(self.HEADER_FORMAT, data[: self.HEADER_SIZE])
        )

        data_iter = itertools.islice(data, self.HEADER_SIZE, None)
        self.signal = {}
        self.trigger = {}

        while True:
            magic_bytes = bytes(itertools.islice(data_iter, self.MAGIC_SIZE))
            if not magic_bytes:
                break

            magic = struct.unpack("<I", magic_bytes)[0]
            key = magic_pmu.get(magic, "UNKNOWN")

            if key == "END":
                break

            period = struct.unpack(
                "<I", bytes(itertools.islice(data_iter, self.PERIOD_SIZE))
            )[0]
            num_pts = int(self.header.duration / period)

            block_bytes = bytes(itertools.islice(data_iter, self.DATA_SIZE * num_pts))
            block = np.frombuffer(block_bytes, dtype=np.uint32)

            if key == "EVENT":
                block = block & 0xFF
                self.trigger = {k: block == 2 for k in pmu_magic}
            else:
                block = block.astype(float)
                if key in ["ECG1", "ECG2", "ECG3", "ECG4"]:
                    block -= 2048.0  # units in mV
                else:
                    block /= 4095.0

            self.signal[key] = block

    def get_timestamp(self, key):
        num_pts = len(self.signal[key])
        return (
            self.header.timestamp
            + np.linspace(0, self.header.duration/10, num_pts, endpoint=False) / 2.5
        )

    def get_time(self, key):  # in seconds
        return self.get_timestamp(key) * 2.5e-3


class PMU:
    def __init__(self, mdbs):
        # member variables that will be populated
        self.signal = {}
        self.trigger = {}
        self.timestamp = {}
        self.pmublocks = []  # store blocks

        for mdb in mdbs:
            if not mdb.is_flag_set("SYNCDATA"):
                continue
            seqdata = SeqData(mdb.data)
            if not seqdata.hdr.id.startswith(b"PMU"):
                continue
            is_learning_phase = seqdata.hdr.id.startswith(b"PMULearnPhase")
            block = PMUblock(seqdata.data)
            self.pmublocks.append(block)
            for key in block.signal:
                pmu_key = key
                if is_learning_phase:
                    pmu_key = f"LEARN_{pmu_key}"
                if pmu_key not in self.signal:
                    self.signal[pmu_key] = []
                    self.trigger[pmu_key] = []
                    self.timestamp[pmu_key] = []
                self.signal[pmu_key].append(block.signal[key])
                self.trigger[pmu_key].append(block.trigger[key])
                self.timestamp[pmu_key].append(block.get_timestamp(key))
        for pmu_key in self.signal:
            self.signal[pmu_key] = np.concatenate(self.signal[pmu_key])
            self.trigger[pmu_key] = np.concatenate(self.trigger[pmu_key])
            self.timestamp[pmu_key] = np.concatenate(self.timestamp[pmu_key])

    def __str__(self):
        """Convert to string, for str()."""
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}:\n"
            f"  .signal: dict of pmu waveforms\n"
            f"  .trigger: dict of triggers for each channel\n"
            f"  .timestamp: dict of timestamps for each channel"
        )

    def plot(self, keys=None, show_trigger=True):
        if keys is None:
            keys = ["ECG1", "ECG2", "ECG3", "ECG4", "PULS", "RESP"]
        elif isinstance(keys, str):
            keys = [keys]

        if show_trigger:
            trig_keys = [key for key in keys[:4] if np.any(self.trigger[key])]
            if not trig_keys:
                print("No trigger signals found")
                show_trigger = False

        _, axs = plt.subplots(1 + bool(show_trigger), 1, squeeze=False, sharex=True)
        colors = {}
        for key in keys:
            axs[0, 0].plot(self.timestamp[key], self.signal[key], label=key)
            colors[key] = axs[0, 0].lines[-1].get_color()

        axs[-1, 0].set_xlabel("timestamp [2.5 us ticks from midnight]")
        axs[0, 0].set_ylabel("normalized signal")
        axs[0, 0].legend()

        # add secondary x-axis with time in seconds
        t0 = self.get_time(keys[0])[0]
        secax = axs[0, 0].secondary_xaxis(
            "top", functions=(lambda x: x * 2.5e-3 - t0, lambda x: (x + t0) / 2.5e-3)
        )
        secax.set_xlabel("time [s]")

        if show_trigger:
            color = [colors[key] for key in trig_keys]
            event = [self.timestamp[key][self.trigger[key]] for key in trig_keys]
            axs[1, 0].eventplot(event, linelengths=0.8, color=color)
            axs[1, 0].legend(trig_keys)
            axs[1, 0].set_ylabel("trigger signals")

        plt.show()

    def get_time(self, key):  # in s
        return self.timestamp[key] * 2.5e-3

    def get_signal(self, key, timestamp):
        if timestamp < self.timestamp[key][0] or timestamp > self.timestamp[key][-1]:
            return np.nan
        return np.interp(timestamp, self.timestamp[key], self.signal[key])

    def get_trigger(self, key, timestamp):
        if timestamp < self.timestamp[key][0] or timestamp > self.timestamp[key][-1]:
            return np.nan
        # Find next neighbor index
        index = np.searchsorted(self.timestamp[key], timestamp)
        # Ensure the index is within the bounds of the array
        index = np.clip(index, 0, len(timestamp[key]) - 1)
        return self.trigger[key][index]
