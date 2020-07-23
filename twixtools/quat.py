# helper functions to convert quaternions to read/phase/slice normal vectors
# and vice vectors
# direct translation from ismrmrd.c from the ismrmrd project

import numpy as np


def quaternion_to_directions(quat):
    a, b, c, d = quat

    read_dir = 3 * [None]
    phase_dir = 3 * [None]
    slice_dir = 3 * [None]

    read_dir[0] = 1. - 2. * (b * b + c * c)
    phase_dir[0] = 2. * (a * b - c * d)
    slice_dir[0] = 2. * (a * c + b * d)

    read_dir[1] = 2. * (a * b + c * d)
    phase_dir[1] = 1. - 2. * (a * a + c * c)
    slice_dir[1] = 2. * (b * c - a * d)

    read_dir[2] = 2. * (a * c - b * d)
    phase_dir[2] = 2. * (b * c + a * d)
    slice_dir[2] = 1. - 2. * (a * a + b * b)

    return read_dir, phase_dir, slice_dir


def directions_to_quaternion(read_dir, phase_dir, slice_dir):

    r11, r21, r31 = read_dir
    r12, r22, r32 = phase_dir
    r13, r23, r33 = slice_dir

    a, b, c, d, s = 1, 0, 0, 0, 0
    trace = 0

    # verify the sign of the rotation
    if __sign_of_directions(read_dir, phase_dir, slice_dir) < 0:
        # flip 3rd column
        r13, r23, r33 = -r13, -r23, -r33

    # Compute quaternion parameters
    # http://www.cs.princeton.edu/~gewang/projects/darth/stuff/quat_faq.html#Q55
    trace = 1.0 + r11 + r22 + r33
    if trace > 0.00001:  # simplest case
        s = np.sqrt(trace) * 2
        a = (r32 - r23) / s
        b = (r13 - r31) / s
        c = (r21 - r12) / s
        d = 0.25 * s
    else:
        # trickier case...
        # determine which major diagonal element has
        # the greatest value...
        xd = 1.0 + r11 - (r22 + r33)  # 4**b**b
        yd = 1.0 + r22 - (r11 + r33)  # 4**c**c
        zd = 1.0 + r33 - (r11 + r22)  # 4**d**d
        # if r11 is the greatest
        if xd > 1.0:
            s = 2.0 * np.sqrt(xd)
            a = 0.25 * s
            b = (r21 + r12) / s
            c = (r31 + r13) / s
            d = (r32 - r23) / s
        # else if r22 is the greatest
        elif yd > 1.0:
            s = 2.0 * np.sqrt(yd)
            a = (r21 + r12) / s
            b = 0.25 * s
            c = (r32 + r23) / s
            d = (r13 - r31) / s
        # else, r33 must be the greatest
        else:
            s = 2.0 * np.sqrt(zd)
            a = (r13 + r31) / s
            b = (r23 + r32) / s
            c = 0.25 * s
            d = (r21 - r12) / s

        if a < 0.0:
            a, b, c, d = -a, -b, -c, -d

    return [a, b, c, d]


def __sign_of_directions(read_dir, phase_dir, slice_dir):
    r11, r21, r31 = read_dir
    r12, r22, r32 = phase_dir
    r13, r23, r33 = slice_dir

    # Determinant should be 1 or -1
    deti = (r11 * r22 * r33) + (r12 * r23 * r31) + (r21 * r32 * r13) -\
           (r13 * r22 * r31) - (r12 * r21 * r33) - (r11 * r23 * r32)

    if deti < 0:
        return -1
    else:
        return 1
