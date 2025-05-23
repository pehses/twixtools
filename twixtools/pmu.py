import numpy as np
import struct
import matplotlib.pyplot as plt

events_VA = {
    "Patient Table": 0x00000001,
    "ECG1": 0x00000002,
    "ECG2": 0x00000002,
    "ECG3": 0x00000002,
    "ECG4": 0x00000002,
    "PULS": 0x00000004,
    "RESP": 0x00000008,
    "EXT1": 0x00000010,
    "EXT2": 0x00000020,
    "CardPT": 0x00000040,
    "RespPT": 0x00000080,
}

events_VA61 = {
    "Patient Table": 0x00000001,
    "ECG1": 0x00000002,
    "ECG2": 0x00000002,
    "ECG3": 0x00000002,
    "ECG4": 0x00000002,
    "PULS": 0x00000004,
    "RSR1": 0x00000004,
    "RSR2": 0x00000004,
    "RESP CUSH": 0x00000004,
    "EXT1": 0x00000010,
    "EXT2": 0x00000020,
    "CardPT": 0x00000040,
    "RespPT": 0x00000080,
}


pmu_magic = {
    "END": 0x01FF0000,
    "ECG1": 0x01010000,
    "ECG2": 0x01020000,
    "ECG3": 0x01030000,
    "ECG4": 0x01040000,
    "PULS": 0x01050000,
    "RESP": 0x01060000,
    "EXT1": 0x01070000,
    "EXT2": 0x01080000
}

pmu_magic_VA = {
    "END": 0xFFFFFFFF,
    "ECG1": 0x00000000,
    "ECG2": 0x00000001,
    "ECG3": 0x00000002,
    "ECG4": 0x00000003,
    "PULS": 0x00000004,
    "RESP": 0x00000005,
    "EXT1": 0x00000006,
    "EXT2": 0x00000007,
    "EVENTS": 0x00000008,
    "CardPT": 0x00000009,
    "RespPT": 0x0000000a,
}

pmu_magic_VA61 = {
    "ECG1": 1,
    "ECG2": 2,
    "ECG3": 3,
    "ECG4": 4,
    "PULS": 5,
    "RSR1": 6,
    "RSR2": 7,
    "RESP CUSH": 8,
    "EXT1": 9,
    "EXT2": 10,
    "CardPT": 35,
    "RespPT": 36,
    "EVENTS": 49,
}

magic_pmu = dict(reversed(item) for item in pmu_magic.items())
magic_pmu_VA = dict(reversed(item) for item in pmu_magic_VA.items())
magic_pmu_VA61 = dict(reversed(item) for item in pmu_magic_VA61.items())


class PMUblock():
    def __init__(self, data, unknown_keys=set()):
        self.timestamp0, self.timestamp, self.packet_no, self.duration = struct.unpack('IIII', data[:16])
        i = 16
        self.signal = dict()
        self.trigger = dict()
        while i < len(data):
            magic, = struct.unpack('I', data[i:i+4])
            if magic in magic_pmu:
                key = magic_pmu[magic]
                if key == 'END':
                    break
            elif magic not in unknown_keys:
                print('unknown magic key: ', hex(magic))
                unknown_keys.add(magic)
            i += 4
            period, = struct.unpack('I', data[i:i+4])
            i += 4
            n_pts = int(self.duration/period)
            if magic in magic_pmu:
                block = np.frombuffer(data[i:i+4*n_pts], dtype=np.uint16).reshape((n_pts, 2)).T
                self.signal[key] = block[0].astype(float) / 4096.
                self.trigger[key] = block[1].astype(bool)
            i += 4*n_pts

    def get_timestamp(self, key, trigger=False):
        if trigger:
            n_pts = len(self.trigger[key])
        else:
            n_pts = len(self.signal[key])
        return self.timestamp + np.linspace(0, self.duration/10., n_pts, endpoint=False) / 2.5

    def get_time(self, key):  # in s
        return self.get_timestamp(key) * 2.5e-3

class PMUblockVA():
    # VA before 61
    def __init__(self, data, unknown_keys=set()):
        # Unpack all at once: 4 uint32 values (4 x 4 bytes = 16 bytes)
        self.timestamp0, self.timestamp, self.packet_no, duration_version = struct.unpack('IIII', data[:16])
        
        # Extract duration and version from the 4th uint32
        self.duration = duration_version & 0xFFFFFF       # bits 0–23
        self.version = (duration_version >> 24) & 0xFF     # bits 24–31
        
        i = 16
        self.signal = dict()
        self.trigger = dict()
        while i < len(data):
            magic, = struct.unpack('I', data[i:i+4])
            if magic in magic_pmu_VA:
                key = magic_pmu_VA[magic]
                if key == 'END':
                    break
            elif magic not in unknown_keys:
                print('unknown magic key: ', hex(magic))
                unknown_keys.add(magic)
            i += 4
            period, = struct.unpack('I', data[i:i+4])
            i += 4
            n_pts = int(self.duration/period)
            if magic in magic_pmu_VA:
                block = np.frombuffer(data[i:i+4*n_pts], dtype=np.uint32)
                if key == 'EVENTS':
                    self.ecg_method = (block >> 8) & 0xF
                    event_type = block & (~0xF00 & 0xFFFFFFFF)
                    for signal_key in events_VA:
                        self.trigger[signal_key] = event_type==events_VA[signal_key] 
                elif block.dtype == np.uint32:
                    self.signal[key] = block.astype(float) / 4095.
                else:
                    self.signal[key] = block.astype(float)
            i += 4*n_pts

    def get_timestamp(self, key, trigger=False):
        if trigger:
            n_pts = len(self.trigger[key])
        else:
            n_pts = len(self.signal[key])
        return self.timestamp + np.linspace(0, self.duration/10., n_pts, endpoint=False) / 2.5

    def get_time(self, key):  # in s
        return self.get_timestamp(key) * 2.5e-3


class PMUblockVA61():
    def __init__(self, data, unknown_keys=set()):
        # packet header: 16 bytes including 2 for the version
        version, size_header, _, size_full_packet, self.timestamp = struct.unpack(
                'HBBIQ', data[:16]
        )
        # extended packet header: 24 bytes
        self.packet_id, self.duration, _, active_signal =  struct.unpack('QIIQ', data[16:40])
        i = 40
        self.signal = dict()
        self.trigger = dict()
        while i < len(data):
            # read base signal header: 16 bytes
            # get magic key and size info
            magic, added_bytes, size_one, size_data, nb_samples = struct.unpack('BBHHH', data[i:i+8])
            if magic in magic_pmu_VA61:
                key = magic_pmu_VA61[magic]
            elif magic not in unknown_keys:
                print('unknown magic key: ', hex(magic))
                unknown_keys.add(magic)
            i+=8
            # get timestamp info
            period, delta_timestamp = struct.unpack('II', data[i:i+8])
            i+=8
            
            # read extended signal header
            ref_value, divisor, offset = struct.unpack('ddd', data[i:i+24])
            i+=24
            if magic in magic_pmu_VA61:
                block = np.frombuffer(data[i:i+size_one*nb_samples], dtype=np.float32)
                if magic == 49:
                    for signal_key in events_VA61:
                        self.trigger[signal_key] = block==events_VA61[signal_key]
                else:
                    self.signal[key] = block.astype(float) * divisor + offset
            i += size_one*nb_samples + added_bytes

    def get_timestamp(self, key, trigger=False):
        if trigger:
            n_pts = len(self.trigger[key])
        else:
            n_pts = len(self.signal[key])
        return (self.timestamp + np.linspace(0, self.duration, n_pts, endpoint=False)) / 2.5e3

    def get_time(self, key):  # in s
        return self.get_timestamp(key) * 2.5e-3


class PMU():
    def __init__(self, mdbs, syngo_version=None):
        # member variables that will be populated
        self.signal = dict()
        self.trigger = dict()
        self.timestamp = dict()
        self.timestamp_trigger = dict()
        self.pmublocks = []  # store blocks
        self.unknown_pmu_magic = set()

        for mdb in mdbs:
            if not mdb.is_flag_set('SYNCDATA'):
                continue
            seqdata = mdb.data
            if not seqdata.hdr.id.startswith(b'PMU'):
                continue
            is_learning_phase = seqdata.hdr.id.startswith(b'PMULearnPhase')
            if syngo_version is not None and syngo_version.startswith('XA'):
                if int(syngo_version[2:])>= 61:
                    block = PMUblockVA61(seqdata.data, self.unknown_pmu_magic)
                else:
                    block = PMUblockVA(seqdata.data, self.unknown_pmu_magic)
            else:
                block = PMUblock(seqdata.data, self.unknown_pmu_magic)
            self.pmublocks.append(block)
            for key in block.signal:
                pmu_key = key
                if is_learning_phase:
                    pmu_key = 'LEARN_' + pmu_key
                if pmu_key not in self.signal:
                    self.signal[pmu_key] = []
                    self.trigger[pmu_key] = []
                    self.timestamp[pmu_key] = []
                    self.timestamp_trigger[pmu_key] = []
                self.signal[pmu_key].append(block.signal[key])
                self.trigger[pmu_key].append(block.trigger[key])
                self.timestamp[pmu_key].append(block.get_timestamp(key))
                self.timestamp_trigger[pmu_key].append(block.get_timestamp(key, trigger=True))
        for pmu_key in self.signal:
            self.signal[pmu_key] = np.concatenate(self.signal[pmu_key])
            self.trigger[pmu_key] = np.concatenate(self.trigger[pmu_key])
            self.timestamp[pmu_key] = np.concatenate(self.timestamp[pmu_key])
            self.timestamp_trigger[pmu_key] = np.concatenate(self.timestamp_trigger[pmu_key])

    def __str__(self):
        """Convert to string, for str()."""
        return (f"{self.__class__.__module__}.{self.__class__.__qualname__}:\n"
                f"  .signal: dict of pmu waveforms\n"
                f"  .trigger: dict of triggers for each channel\n"
                f"  .timestamp: dict of timestamps for each channel\n"
                f"  .timestamp_trigger: dict of timestamps for the trigger of each channel")

    def plot(self, keys=None, show_trigger=True, show_learning_phase=False):

        if keys is None:
            if show_learning_phase:
                keys = [key for key in self.signal if np.ptp(self.signal[key])>0]
            else:
                keys = [key for key in self.signal if not key.startswith('LEARN_') and np.ptp(self.signal[key])>0]
        elif isinstance(keys, str):
            keys = [keys]

        if show_trigger:
            trig_keys = [key for key in keys if np.any(self.trigger[key])]
            if len(trig_keys) == 0:
                print('No trigger signals found')
                show_trigger = False

        _, axs = plt.subplots(1 + bool(show_trigger), 1, squeeze=False, sharex=True)
        colors = dict()
        for key in keys:
            # normalize signal
            min_signal, max_signal = self.signal[key].min(), self.signal[key].max()
            if min_signal == max_signal:
                normalized_signal = np.zeros_like(self.signal[key])
            else:
                normalized_signal = (self.signal[key]-min_signal)/(max_signal-min_signal)
                
            axs[0, 0].plot(self.timestamp[key], normalized_signal, label=key)
            colors[key] = axs[0, 0].lines[-1].get_color()

        axs[-1, 0].set_xlabel('timestamp [2.5 us ticks from midnight]')
        axs[0, 0].set_ylabel('normalized signal')
        axs[0, 0].legend()

        # add secondary x-axis with time in seconds
        t0 = self.get_time(keys[0])[0]
        secax = axs[0, 0].secondary_xaxis('top', functions=(lambda x: x*2.5e-3 - t0, lambda x: (x + t0) / 2.5e-3))
        secax.set_xlabel('time [s]')

        if show_trigger:
            color = [colors[key] for key in trig_keys]
            event = [self.timestamp_trigger[key][self.trigger[key]] for key in trig_keys]
            axs[1, 0].eventplot(event, linelengths=0.8, color=color)
            axs[1, 0].legend(trig_keys)
            axs[1, 0].set_ylabel('trigger signals')

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
