# twixtools

## twix file reader (read_twix)

### Purpose

twixtools provide reading and writing capability of Siemens MRI raw data files (.dat). In addition, it also includes the compression utility twixzip (see further below for a description).


### Requirements

The tool works under Python 3.7 with the following package installed:

* numpy &ge; 1.17.3


### Installation

Navigate to the twixtools folder in an open terminal and install twixtools with pip:

    pip install .

Installation through `python setup.py install` is currently not possible.


### Brief Tutorial

The raw data file can be parsed using the read_twix function:

```python
import twixtools
multi_twix = twixtools.read_twix(filename)
```

The function returns a list of individual measurements (of length >=1). The last measurement usually corresponds to the imaging scan, earlier measurements often include calibration data. Each measurement contains a python dict() with the following entries:

* 'mdb': measurement data divided into blocks (return type: list)
* 'hdr_str': dict of protocol header strings (divided into different protocol types)
* 'hdr': dict of partially parsed protocol header strings (each dict element contains another dict with protocol information)
* ('init': required for twix file writing, otherwise of little importance)


Each invididual 'mdb' in the list of mdbs consists of a data and a header (line counters and such) part, which can be accessed as follows:

```python
mdb = multi_twix[-1]['mdb'][0] # first mdb element of last measurement
mdb.data # data of first mdb (may or may not be imaging data)
mdb.mdh # full miniheader information stored as a numpy dtype object
 ```

Different data types can be distinguished by returning a list of active flags, or by directly checking whether the data is assumed to be from an imaging scan (and not from a calibration scan such as a phase correction scan or a noise measurement):

```python
mdb.get_active_flags() # get all active MDH flags
mdb.is_image_scan() # check if this an image scan (True or False)
```

Line Counters can be accessed as follows:
```python
mdb.cLin   # returns line number
mdb.cPar   # returns partition number
mdb.c<tab> # with line completion enabled, this should give you a list of all counters
```

The full minidata header (mdh) information is stored in a `mdb.mdh` special numpy dtype object. You can print a list of its entry names by printing `mdb.mdh.dtype.names`.


### Example code
```python
import numpy as np
import twixtools

# read all image data from file
def read_image_data(filename):
    out = list()
    for mdb in twixtools.read_twix(filename)[-1]['mdb']:
        if mdb.is_image_scan():
            out.append(mdb.data)
    return np.asarray(out)  # 3D numpy array [acquisition_counter, n_channel, n_column]


# read image data from list of mdbs and sort into 3d k-space (+ coil dim.)
def import_kspace(mdb_list)
    image_mdbs = []
    for mdb in mdb_list:
        if mdb.is_image_scan():
            image_mdbs.append(mdb)

    n_line = 1 + max([mdb.cLin for mdb in image_mdbs])
    n_part = 1 + max([mdb.cPar for mdb in image_mdbs])
    n_channel, n_column = image_mdbs[0].data.shape

    out = np.zeros([n_part, n_line, n_channel, n_column], dtype=np.complex64)
    for mdb in image_mdbs:
        # '+=' takes care of averaging, but careful in case of other counters (e.g. echoes)
        out[mdb.cPar, mdb.cLin] += mdb.data

    return out  # 4D numpy array [n_part, n_line, n_channel, n_column]
```


## twixzip Compression Utility

### Purpose

twixzip is a Python based command line tool for Siemens MRI raw data compression. Following compression methods can be selected via the command line:

* Oversampling removal
* Lossy floating point compression using the zfp library
* Single coil compression (scc) based on singular value decomposition (SVD)
* Geometric coil compression (gcc) based on SVD
* Optionally FID navigators can be removed

Before applying the selected compression method(s), lossless compression (gzip) is applied to the header and meta data information which is then added to a hdf5 file. All additional meta information necessary for decompression (e.g. coil compression matrices) are also stored in the hdf5 file.


### Additional Requirements

* pyzfp &ge; 0.3.1
* pytables &ge; 3.6.1

The pyzfp and pytables libraries can be installed via pip:

    pip install pyzfp
    pip install tables


### Usage

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
    --scc -n NCC:               single coil compression (SCC) - keep NCC virtual coils
    --scc -t CC_TOL:            SCC - number of coils is calculated with a tolerance for the singular values
    --scc_bart -n NCC:          SCC using BART  
    --gcc -n NCC:               geometric coil compression (GCC) - keep NCC virtual coils
    --gcc -t CC_TOL:            GCC - number of coils is calculated with a tolerance for the singular values
    --gcc_bart -n NCC:          GCC using the Berkeley Advanved Reconstruction Toolbox (BART) [1]         
    --zfp --zfp_tol ZFP_TOL:    floating point compression with ZFP_TOL tolerance
    --zfp --zfp_prec ZFP_PREC:  floating point compression with ZFP_PREC precision (not recommended)

The optional argument `--testmode` can be used to automatically decompress the data after compression. The created decompressed MRI raw data filename contains the selected compression method. The option `--profile` can be used to profile the compression code.

[1] BART Toolbox for Computational Magnetic Resonance Imaging, DOI: 10.5281/zenodo.592960