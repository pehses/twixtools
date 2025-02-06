#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 15:12:17 2021

@author: ehsesp
"""
import argparse
import numpy as np
import time

import twixtools
from twixtools import map_twix
from twixtools.contrib import cfl
from twixtools.recon_helpers import remove_oversampling

cfl_order = [
    'READ_DIM',
    'PHS1_DIM',
    'PHS2_DIM',
    'COIL_DIM',
    'MAPS_DIM',
    'TE_DIM',
    'COEFF_DIM',
    'COEFF2_DIM',
    'ITER_DIM',
    'CSHIFT_DIM',
    'TIME_DIM',
    'TIME2_DIM',
    'LEVEL_DIM',
    'SLICE_DIM',
    'AVG_DIM',
    'BATCH_DIM'
]

dim_map = {
    'Ide': 'MAPS_DIM',
    'Idd': 'COEFF2_DIM',
    'Idc': 'COEFF_DIM',
    'Idb': 'ITER_DIM',
    'Ida': 'LEVEL_DIM',
    'Seg': 'BATCH_DIM',
    'Set': 'CSHIFT_DIM',
    'Rep': 'TIME_DIM',
    'Phs': 'TIME2_DIM',
    'Eco': 'TE_DIM',
    'Par': 'PHS2_DIM',
    'Sli': 'SLICE_DIM',
    'Ave': 'AVG_DIM',
    'Lin': 'PHS1_DIM',
    'Cha': 'COIL_DIM',
    'Col': 'READ_DIM'
}


def read_to_cfl(twixmap, data_type, remove_os, no_average):
    if data_type.lower() == 'noise':
        # noise data is treated differently:
        #   -- concatenate all noise ADC into one column
        sig = list()
        for mdb in twixmap.mdb_list:
            if remove_os:
                sig.append(remove_oversampling(mdb.data)[0])
            else:
                sig.append(mdb.data)
        sig = np.moveaxis(np.asarray(sig), 1, 0)
        sig = sig.reshape([sig.shape[0], -1])

        # move all additional dims to READ_DIM
        sz = np.ones(1, twixmap.size.dtype)[0]
        sz['Cha'], sz['Col'] = sig.shape
        sig = sig.reshape(sz)
    else:
        twixmap.flags['remove_os'] = remove_os
        twixmap.flags['average']['Seg'] = not no_average
        if no_average:
            twixmap.flags['average']['Ave'] = False
        if data_type.lower() == 'image':
            # image data is zero-filled to the correct matrix size
            twixmap.flags['zf_missing_lines'] = True
        else:
            # all other data is cropped at the min and max line/partition indices
            twixmap.flags['skip_empty_lead'] = True
        sig = twixmap[:]

    return np.moveaxis(sig, np.arange(twixmap.ndim), [cfl_order.index(v) for v in dim_map.values()])


def convert_to_cfl(twix_filename, out_filename, meas_no, data_type, remove_os, apply_rawdatacorr, no_average):

    print('\n ---- Parsing twix file ----\n')
    t_parse_start = time.time()
    twixmap = map_twix(twix_filename)

    if apply_rawdatacorr:
        twixmap.setCoilInfoTo('nova_ptx')  # wip, these should be read from protocol

    n_meas = len(twixmap)
    if meas_no+1 > n_meas:
        raise IndexError('meas_no is out of range')
    elif n_meas > 1:
        print(f'\n ---- Selecting measurement {slice(meas_no).indices(n_meas)[1]+1} out of {n_meas} ----')
    twixmap = twixmap[meas_no]

    if data_type == 'all':
        for item in twixmap:
            if isinstance(twixmap[item], twixtools.twix_array):
                fname = out_filename + '_' + item
                print(f'\n --- Now reading {item} data and writing to cfl file {fname} ---')
                sig = read_to_cfl(twixmap[item], item, remove_os, no_average)
                cfl.writecfl(fname, sig)
    else:
        print(f'\n --- Now reading {data_type} data and writing to cfl file {out_filename}')
        sig = read_to_cfl(twixmap[data_type], data_type, remove_os, no_average)
        cfl.writecfl(out_filename, sig)
    t_convert_stop = time.time()
    t = t_convert_stop - t_parse_start
    print(f'\n--- Completed in {t//60} min and {t%60} s ---')


def main():
    parser = argparse.ArgumentParser(description='convert Siemens twix file to bart cfl file')
    parser.add_argument('infile', type=str, help='Siemens twix input file')
    parser.add_argument('outfile', type=str, help='bart cfl output filename (without extension)')
    parser.add_argument('--meas_no', type=int, default=-1, help='measurement number in case of multi-raid file (default: -1)')
    parser.add_argument('--remove_os', action='store_true', help='remove oversampling (factor 2)')
    parser.add_argument('--rawdatacorr', action='store_true', help='apply rawdata correction factors (wip)')
    parser.add_argument('--no-average', dest='no_average', action='store_true', help='No averaging along BATCH or AVG dim.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--type', dest='data_type', type=str, default='image', help='raw data type (default: image)')
    group.add_argument('--read_all', action='store_true', help='reads all known data types, outfile name is used as prefix')
    args = parser.parse_args()

    if args.read_all:
        args.data_type = 'all'

    convert_to_cfl(args.infile, args.outfile, args.meas_no, args.data_type, args.remove_os, args.rawdatacorr, args.no_average)


if __name__ == "__main__":
    main()
