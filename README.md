# twixtools

[![Python package](https://github.com/pehses/twixtools/actions/workflows/test_action.yml/badge.svg)](https://github.com/pehses/twixtools/actions/workflows/test_action.yml)

### Purpose

twixtools provide reading and limited writing capability of Siemens MRI raw data files (.dat).


### Installation

Navigate to the twixtools folder in an open terminal and install twixtools with pip:

    pip install .

Installation through `python setup.py install` is currently not possible.


### Demo code

A jupyter notebook that demonstrates the basic functionality of the  `read_twix`, `map_twix`, and `write_twix` tools can be found in `demo/recon_example.ipynb`.


## read_twix: "low-level" access to twix data

The raw data file can be parsed using the read_twix function:

```python
import twixtools
multi_twix = twixtools.read_twix(filename)
```

The function returns a list of individual measurements (of length >=1). The last measurement usually corresponds to the imaging scan, earlier measurements often include calibration data. Each measurement contains a python dict() with the following entries:

* 'mdb': measurement data divided into blocks (return type: list)
* 'hdr': dict of parsed protocol header strings (each dict element contains another dict with protocol information)
* 'hdr_str': dict of original protocol header strings (divided into different protocol types)
  - note that this is the protocol information that is used for twix file writing (by `write_twix`), so make sure to make necessary adjustments here and not in ['hdr']
* 'pmu': physiological (PMU) data (if available and parse_pmu is set to True)
* ('raidfile_hdr': required for twix file writing, otherwise of little importance)


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

## map_twix: "high level" access to twix data
`map_twix` is a high-level function that takes the data obtained from `read_twix` (in the form of `Mdb` objects), and maps it to multi-dimensional "k-space" arrays. These `twix_array` objects are generated for different data types (image/noise adjust/phase-correction/... scan) and can be accessed with `numpy.ndarray` array-slicing syntax.


Optional flags control additional feature and also have an impact on size and shape of the multidimensional arrays. The following flags are currently available (stored in the `flags` dict within each `twix_array` object):
  * `average`: dict of bools that determines which dimensions should be averaged.
  * `squeeze_ave_dims`: bool that determines whether averaged dimensions should be removed/squeezed from the array's shape.
  * `remove_os`: oversampling removal. Reduces the number of columns by a factor of two.
  * `regrid`: bool that controls ramp-sampling regridding (if applicable)
  * `skip_empty_lead`: skips to first line & partition that is found in mdb list (e.g. if first line counter is 10, the output array starts at line counter 10).
  * `zf_missing_lines`: zero-fill k-space to include lines and partitions that are higher than the maximum counter found in the mdb list, but are still within the k-space matrix according to the twix header.

If available, physiological (PMU) data is stored in the returned dict under the 'pmu' key.

For example code, please look at the `demo/recon_example.ipynb` jupyter file.

## Acknowledgements

The protocol header parsing code originates from William Clarke's excellent pymapvbvd project (https://github.com/wexeee/pymapvbvd).
