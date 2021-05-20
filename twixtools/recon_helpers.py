import numpy as np


def to_freqdomain(data, x_in_timedomain=True, axis=-1):
    if not x_in_timedomain:
        return data, False
    else:
        return np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(data, axes=[axis]), axis=axis), axes=[axis]), False


def to_timedomain(data, x_in_timedomain=False, axis=-1):
    if x_in_timedomain:
        return data, True
    else:
        return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(data, axes=[axis]), axis=axis), axes=[axis]), True

def remove_oversampling(data, x_was_in_timedomain=True):
    nx = data.shape[-1]
    data, x_in_timedomain = to_freqdomain(data, x_was_in_timedomain)
    data = data[..., nx//4:nx*3//4]
    if x_was_in_timedomain:
        data, x_in_timedomain = to_timedomain(data, x_in_timedomain)
    return data, x_in_timedomain