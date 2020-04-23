# twixzip</span>.py Installation & Usage

## Purpose

twixzip is a Python based command line tool for Siemens MRI raw data compression. Following compression methods can be selected via the command line:

* Oversampling removal
* Lossy floating point compression using the zfp library
* Single coil compression (scc) based on svd
* Geometric coil compression (gcc) based on svd
* Optionally FID navigators can be removed

Before applying the selected compression method(s), lossless compression (gzip) is applied to the header and meta data information which is then added to a hdf5 file. All additional meta information necessary for decompression (e.g. coil compression matrices) are also stored in the hdf5 file.

## Installation

Navigate to the twixtools folder in an open terminal and install twixzip</span>.py with pip:

    pip install .

Installation through python setup</span>.py install is currently not possible.

## Requirements

The tool works under Python 3.7 with the following packages installed:

* numpy &ge; 1.17.3
* pyzfp &ge; 0.3.1
* h5py &ge; 2.10.0

The pyzfp and h5py libraries can be installed via pip:

    pip install pyzfp
    pip install h5py

## Usage

Executing the command twixzip</span>.py in an open terminal gives an overview of all possible arguments. Optional arguments are:

    -h:  help  
    -d:  decompress data

Input and output directories & filenames are required arguments that can be selected via:

    -i infile:  input file  
    -o outfile: output file

In the compression mode the input file should be an MRI raw data file, in the decompression mode (`-d`) it should be the hdf5 file containing the compressed data. The output file is then an hdf5 file (compression mode) or an MRI raw data file (decompression mode).

Compression methods can be selected via:

    --remove_fidnav:            removes FID navigators  
    --remove_os:                removes oversampling  
    --gcc -n NCC:               geometric coil compression with NCC virtual coils
    --scc -n NCC:               single coil compression with NCC virtual coils
    --zfp --zfp_tol ZFP_TOL:    floating point compression with ZFP_TOL tolerance
    --zfp --zfp_prec ZFP_PREC:  floating point compression with ZFP_PREC precision (not recommended)

The optional argument `--testmode` can be used to automatically decompress the data after compression. The created decompressed MRI raw data filename contains the selected compression method.
