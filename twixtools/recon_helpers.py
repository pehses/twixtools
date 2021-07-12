import numpy as np
from scipy.integrate import cumtrapz


def to_freqdomain(data, x_in_timedomain=True, axis=-1):
    if not x_in_timedomain:
        return data, False
    else:
        return np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(data, axes=[axis]),
                               axis=axis), axes=[axis]), False


def to_timedomain(data, x_in_timedomain=False, axis=-1):
    if x_in_timedomain:
        return data, True
    else:
        return np.fft.fftshift(np.fft.fft(
            np.fft.ifftshift(data, axes=[axis]), axis=axis), axes=[axis]), True


def remove_oversampling(data, x_was_in_timedomain=True):
    nx = data.shape[-1]
    data, x_in_timedomain = to_freqdomain(data, x_was_in_timedomain)
    data = data[..., nx//4:nx*3//4]
    if x_was_in_timedomain:
        data, x_in_timedomain = to_timedomain(data, x_in_timedomain)
    return data, x_in_timedomain


def calc_regrid_traj(prot):

    meas = prot['Meas']

    if 'alRegridMode' not in meas or meas['alRegridMode'][0] < 2:
        return None

    regrid_mode = meas['alRegridMode'][0]
    ncol = meas['alRegridDestSamples'][0]
    dwelltime = meas['aflRegridADCDuration'][0] / ncol
    start = meas['alRegridDelaySamplesTime'][0]
    rampup_time = meas['alRegridRampupTime'][0]
    flattop_time = meas['alRegridFlattopTime'][0]
    rampdown_time = meas['alRegridRampdownTime'][0]
    gr_adc = np.zeros(ncol, dtype=np.single)
    time_adc = start + dwelltime * np.arange(0.5, ncol + 0.5)
    ixUp = np.where(time_adc < rampup_time)[0]
    ixFlat = np.setdiff1d(np.where(time_adc <= rampup_time + flattop_time)[0],
                          np.where(time_adc < rampup_time)[0])
    ixDn = np.setdiff1d(np.setdiff1d(np.arange(ncol), ixFlat), ixUp)
    gr_adc[ixFlat] = 1
    if regrid_mode == 2:
        # trapezoidal gradient
        gr_adc[ixUp] = time_adc[ixUp] / rampup_time
        gr_adc[ixDn] = 1 - (time_adc[ixDn] - rampup_time - flattop_time)\
            / rampdown_time
    elif regrid_mode == 4:
        gr_adc[ixUp] = np.sin(np.pi / 2 * time_adc[ixUp] / rampup_time)
        gr_adc[ixDn] = np.sin(
            np.pi / 2 * (1 + (time_adc[ixDn] - rampup_time - flattop_time)
                         / rampdown_time))
    else:
        raise Exception('regridding mode unknown')

    # make sure that gr_adc is always positive
    # (rs_traj needs to be strictly monotonic)
    gr_adc = np.maximum(gr_adc, 1e-4)
    rs_traj = (np.append(0, cumtrapz(gr_adc)) - ncol // 2) / np.sum(gr_adc)
    rs_traj -= np.mean(rs_traj[ncol//2-1:ncol//2+1])

    # scale rs_traj by kmax (only works if all slices have same FoV!!!)
    kmax = prot['MeasYaps']['sKSpace']['lBaseResolution']\
        / prot['MeasYaps']['sSliceArray']['asSlice'][0]['dReadoutFOV']
    rs_traj *= kmax

    return rs_traj


def perform_regrid(data, rs_traj, ro_shift=0):

    if rs_traj is None:
        return data

    # first correct for readout shifts
    # the nco frequency is always scaled to the max.
    # gradient amp and does account for ramp-sampling
    ncol = len(rs_traj)
    deltak = abs(rs_traj[ncol//2] - rs_traj[ncol//2+1])
    tmp = data * np.exp(1j*2*np.pi*ro_shift*(deltak*np.arange(ncol) - rs_traj))

    # interpolate
    x = np.linspace(rs_traj[0], rs_traj[-1], ncol)
    return np.asarray([np.interp(x, rs_traj, y) for y in tmp])
