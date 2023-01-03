import numpy as np


class ObservableDict(dict):
    """
    Simple class that holds an argument dictionary.

    Triggers a call to .notify_change() of the observer class(es) whenever
    a key is changed (not when it is initialized)
    """

    def __init__(self, *observers, **kwds):
        self.__dict__.update(kwds)
        self.observes = observers

    def __setitem__(self, key, value):
        if key in self:
            dict.__setitem__(self, key, value)
            for obs in self.observes:
                obs.notify_change()
        else:
            dict.__setitem__(self, key, value)


def idea_version_check(f):
    prev_pos = f.tell()
    f.seek(0)
    firstInt, secondInt = np.fromfile(f, dtype=np.uint32, count=2)
    if (firstInt < 10000) and (secondInt <= 64):
        version_is_ve = True
        NScans = secondInt
    else:
        version_is_ve = False
        NScans = 1
    f.seek(prev_pos)

    return version_is_ve, NScans
