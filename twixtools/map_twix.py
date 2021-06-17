import copy
import numpy as np
import twixtools
from twixtools.recon_helpers import (
    remove_oversampling, calc_regrid_traj, perform_regrid
)


# define categories in which the twix data should be sorted based on MDH flags
# that must or must not be set (True/False)
# only one 'any' per category allowed (however, it is possible to add other
# appropriate functions (even synonyms for any))
twix_category = {
    'image':         {'RTFEEDBACK': False, 'HPFEEDBACK': False,
                      'REFPHASESTABSCAN': False, 'PHASESTABSCAN': False,
                      'PHASCOR': False, 'NOISEADJSCAN': False,
                      'noname60': False},
    'noise':         {'NOISEADJSCAN': True},
    'phasecorr':     {'PHASCOR': True, 'PATREFSCAN': False, 'noname60': False},
    'phasestab':     {'PHASESTABSCAN': True, 'REFPHASESTABSCAN': False,
                      'noname60': False,
                      any: {'PATREFSCAN': False, 'PATREFANDIMASCAN': True}},
    'refphasestab0': {'REFPHASESTABSCAN': True, 'PHASESTABSCAN': False,
                      'noname60': False,
                      any: {'PATREFSCAN': False, 'PATREFANDIMASCAN': True}},
    'refphasestab1': {'REFPHASESTABSCAN': True, 'PHASESTABSCAN': True,
                      'noname60': False,
                      any: {'PATREFSCAN': False, 'PATREFANDIMASCAN': True}},
    'refscan':       {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True},
                      'PHASCOR': False, 'PHASESTABSCAN': False,
                      'REFPHASESTABSCAN': False, 'RTFEEDBACK': False,
                      'HPFEEDBACK': False, 'noname60': False},
    'ref_pc':        {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True},
                      'PHASCOR': True},
    'ref_ps':        {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True},
                      'REFPHASESTABSCAN': False, 'PHASESTABSCAN': True},
    'ref_ps_ref0':   {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True},
                      'REFPHASESTABSCAN': True, 'PHASESTABSCAN': False},
    'ref_ps_ref1':   {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True},
                      'REFPHASESTABSCAN': True, 'PHASESTABSCAN': True},
    'rt_feedback':   {any: {'RTFEEDBACK': True, 'HPFEEDBACK': True},
                      'MDH_VOP': False},
    'vop':           {'MDH_VOP': True},
    'fidnav':        {'noname60': True}  # why we include the 'noname60' checks
}


def map_twix(input):
    """ creates a list of measurements (or a single dict if input was dict)
    with data for each measurement mapped to a twix_array object.

    Parameter
    ----------
    input: string, int, list, or dict
        If the filename or its measurement id are passed as a string or int,
        respectively, the corresponding twix file is first parsed using
        `read_twix`. Alternatively, it is possible to directly pass a scan list
        (as returned by `read_twix`) to `map_twix` or to pass only a dict that
        includes header information and mdb list of a single twix scan.

    Returns:
    ----------
    out: dict of twix_array objects
        A twix_array object is created for each data category (as defined by
        `twix_category`) that is encountered in the input.
        The twix_array object includes the header information (twix_array.hdr)
        as well as access to the underlying data via array slicing of a virtual
        'k-space'-like array that is designed to closely mimick a
        `numpy.ndarray` object (and indeed returns a `numpy.ndarray`).

    Examples:
    ----------
    Read the data and then select only the twix_array object that contains
    image data:
    >>> twix = map_twix(filename)
    >>> im_array = twix['image']

    Now set a few optional flags that control additional features and determine
    the shape of the output array:
    >>> im_array.flags['remove_os'] = True  # activate automatic os removal
    >>> im_array.flags['regrid'] = True  # activate ramp sampling regridding
    >>> im_array.flags['average']['Rep'] = True  # average all repetitions
    >>> im_array.flags['squeeze_ave_dims'] = True  # reduces the dimensionality

    Print all available flags and their values:
    >>> print(im_array.flags)

    Print the shape of the data and the names of the active dimensions:
    >>> print(im_array.shape)
    >>> print(im_array.dims)

    And finally read the data:
    >>> im_data = im_array[:]

    Alternatively, we can for example only select the data for the first
    receiver channel:
    >>> im_data0 = im_array[...,0,:]

    All standard array slicing operations should be supported.
    """

    if isinstance(input, list):
        # assume list of measurements
        twix = input
    elif isinstance(input, dict):
        # assume measurement dict
        # return twix_array of the input (no measurement list)
        twix = [input]
    else:
        # assume that this is the filename or a meas id
        twix = twixtools.read_twix(input)

    out = list()
    for meas in twix:

        if not isinstance(meas, dict):
            continue  # first "meas" may store first 10240 bytes of file

        # append new dict to output list
        out.append(dict())

        # sort mdbs into categories
        for mdb in meas['mdb']:
            if mdb.is_flag_set('SYNCDATA'):  # ignore syncdata
                continue
            if mdb.is_flag_set('ACQEND'):
                break

            for category, rqmts in twix_category.items():
                include_in_cat = True
                for flag in rqmts.keys():
                    if isinstance(flag, str):  # check whether flag is set
                        if mdb.is_flag_set(flag) != rqmts[flag]:
                            include_in_cat = False
                            break
                    else:  # assume a function call (probably any())
                        checks = list()
                        for flag2 in rqmts[flag].keys():
                            checks.append(
                                mdb.is_flag_set(flag2) == rqmts[flag][flag2])
                        if not flag(checks):
                            include_in_cat = False
                            break
                if include_in_cat:
                    if category not in out[-1]:
                        out[-1][category] = list()
                    out[-1][category].append(mdb)

        # convert each categories' mdb list to twix_array
        for category in out[-1].keys():
            out[-1][category] = twix_array(out[-1][category],
                                           meas['hdr'].copy())

        # set few default flag(s) for 'image' category
        if 'image' in out[-1]:
            out[-1]['image'].flags['zf_missing_lines'] = True

        # include hdr in dict
        out[-1]['hdr'] = meas['hdr'].copy()
        out[-1]['hdr_str'] = meas['hdr_str'].copy()

    # go back to dict if input was dict
    if isinstance(input, dict):
        out = out[0]

    return out


class twix_array():
    """Memory-mapped storage class for Siemens MRI raw data.

    The twix array class constructs a virtual multi-dimensional array from a
    list of mdb objects, that tries to closely resemble a numpy.ndarray with
    standard array slicing operations. The selected array is then read from the
    twix file (all reading operations are handled by the Mdb class) and
    returned in the form of a multi-dimensional numpy.ndarray.

    Note that additional flags can change the shape of the virtual array.

    Important Attributes
    ----------
    ndim: int
        number of output dimensions. May change depending on `flags`.
    shape: tuple
        shape of the output array. May change depending on `flags`.
    dims: list
        List of names of output dimensions. May change depending on `flags`.
    non_singleton_dims: list
        Returns list of non-singleton dimensions.
    dim_order: tuple
        List of the standard dimension order (immutable).
    hdr: dict
        twix header information
    flags: dict
        Dict of optional flags. The following flags are currently supported:
        - 'average': dict of bools that determines which dimensions should
            be averaged.
        - 'squeeze_ave_dims': bool that determines whether averaged
            dimensions should be removed/squeezed from the array's shape.
        - 'remove_os': oversampling removal. Reduces the number of columns
            by a factor of two.
        - 'regrid': bool that controls ramp-sampling regridding (if applicable)
        - 'skip_empty_lead': skips to first line & partition that is found
            in mdb list (e.g. if first line counter is 10, the output array
            starts at line counter 10).
        - 'zf_missing_lines': zero-fill k-space to include lines and partitions
           that are higher than the maximum counter found in the mdb list, but
           are still within the k-space matrix according to the twix header.
    """

    def __init__(self, mdb_list, hdr=None, flags=None):

        self.mdb_list = mdb_list.copy()
        self.hdr = None
        if hdr is not None:
            self.hdr = copy.deepcopy(hdr)

        self.rs_traj = calc_regrid_traj(self.hdr)

        # delete 'ACQEND' and 'SYNCDATA' flags if present
        twixtools.del_from_mdb_list(
            self.mdb_list,
            lambda b: b.is_flag_set('ACQEND') or b.is_flag_set('SYNCDATA'))

        self._dim_order = (
            "Ide", "Idd", "Idc", "Idb", "Ida", "Seg", "Set", "Rep",
            "Phs", "Eco", "Par", "Sli", "Ave", "Lin", "Cha", "Col"
        )

        # dtype that includes all dims:
        self.dt_dims = np.dtype([(name, "<u2") for name in self.dim_order])

        # dtype that only includes counters (no cha & col)
        self.dt_counters = np.dtype([(n, "<u2") for n in self.dim_order[:-2]])

        self.key_map = {
            'Ide': 'ushIde', 'Idd': 'ushIdd', 'Idc': 'ushIdc',
            'Idb': 'ushIdb', 'Ida': 'ushIda', 'Seg': 'ushSeg',
            'Set': 'ushSet', 'Rep': 'ushRepetition', 'Phs': 'ushPhase',
            'Eco': 'ushEcho', 'Par': 'ushPartition', 'Sli': 'ushSlice',
            'Ave': 'ushAcquisition', 'Lin': 'ushLine'
        }

        self.sorted_mdh_keys = [self.key_map[d] for d in self.dim_order[:-2]]

        # determine k-space shape by finding max index
        shp = np.ones(len(self.dt_dims), dtype=self.dt_dims[1])
        self._first_ix = 1024 * np.ones(len(self.dt_dims)-2,
                                        dtype=self.dt_dims[1])

        for mdb in self.mdb_list:
            sLC = mdb.mdh['sLC']
            sLC = np.asarray(sLC[self.sorted_mdh_keys].tolist(),
                             dtype=sLC[0].dtype)
            req_shape = 1 + sLC
            # add channels & columns
            req_shape = np.concatenate([req_shape,
                                        [mdb.mdh['ushUsedChannels'],
                                         mdb.mdh['ushSamplesInScan']]])
            shp = np.maximum(shp, req_shape)
            self._first_ix = np.minimum(self._first_ix, sLC)

        self.base_size = np.ones(1, dtype=self.dt_dims)[0]
        for key, item in enumerate(shp):
            # complicated, can we do this conversion better? (proper casting?)
            self.base_size[key] = item

        self._flags = {'average': {item: False for item in self.dim_order},
                       'remove_os': False,
                       'regrid': False,
                       'squeeze_ave_dims': False,
                       'skip_empty_lead': False,
                       'zf_missing_lines': False}

        # 'Ave' should be averaged by default, Idx indices should be ignored:
        for dim in ['Ide', 'Idd', 'Idc', 'Idb', 'Ida', 'Ave']:
            self._flags['average'][dim] = True

        # set flags that were passed in constructor call
        if flags is not None:
            for key, item in flags.items():
                try:
                    self.flags[key] = item.copy()
                except Exception:
                    self.flags[key] = item

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        self._flags = self._flags.copy()
        return twix_array(self.mdb_list, self.hdr, self.flags)

    @property
    def dim_order(self):
        return self._dim_order

    @property
    def dims(self):
        if not self.flags['squeeze_ave_dims']:
            return self.dim_order
        else:
            return [n for n in self.dim_order if not self.flags['average'][n]]

    @property
    def non_singleton_dims(self):
        return [dim for dim in self.dim_order if self.size[dim] > 1]

    @property
    def ndim(self):
        return len(self.dims)

    @property
    def flags(self):
        # wip: although the flags dict itself is write-protected, its entries
        #      are currently not and can be overwritten by garbage!
        return self._flags

    @property
    def size(self):
        # self.size returns the shape of the data as a dtype with named
        # elements for easier access
        # averaged dims will be kept even if 'squeeze_ave_dims' is set to True
        sz = self.base_size.copy()
        if not self.flags['average']['Col'] and self.flags['remove_os']:
            sz[self.dim_order.index('Col')] //= 2

        if self.hdr is not None and self.flags['zf_missing_lines']:
            if not self.flags['average']['Lin']:
                hdr_lin = \
                    self.hdr['MeasYaps']['sKSpace']['lPhaseEncodingLines']
                sz['Lin'] = max(sz['Lin'], hdr_lin)
            if not self.flags['average']['Par']\
                    and self.hdr['MeasYaps']['sKSpace']['ucDimension'] > 2:
                hdr_par = self.hdr['MeasYaps']['sKSpace']['lPartitions']
                sz['Par'] = max(sz['Par'], hdr_par)

        if self.flags['skip_empty_lead']:
            if not self.flags['average']['Lin']:
                sz['Lin'] -= self._first_ix[self.dim_order.index('Lin')]
            if not self.flags['average']['Par']:
                sz['Par'] -= self._first_ix[self.dim_order.index('Par')]

        for dim in range(len(sz)):
            if self.flags['average'][self.dim_order[dim]]:
                sz[dim] = 1
        return sz

    @property
    def shape(self):
        # self.shape is the more numpy compatible version of self.size by
        # returning a tuple
        # 'squeeze_ave_dims': averaged dims are removed from shape
        if not self.flags['squeeze_ave_dims']:
            return self.size.item()
        else:
            return [sz for sz, name in zip(self.size.item(),
                    self.size.dtype.names) if not self.flags['average'][name]]

    def __getitem__(self, index):
        # implement array slicing here
        # returns numpy.ndarray
        if not isinstance(index, tuple):
            index = (index,)  # make sure to pass along tuple
        if len(index) > self.ndim:
            raise IndexError(
                "too many indices for array: array is %d-dimensional, "
                "but %d were indexed" % (self.ndim, len(index)))
        ellipsis_in_index = False
        selection = list()
        remove_dim = list()
        for key, item in enumerate(index):
            if ellipsis_in_index:
                key += self.ndim - len(index)
            if item is Ellipsis:
                if ellipsis_in_index:
                    raise IndexError(
                        "an index can only have a single ellipsis ('...')")
                ellipsis_in_index = True
                # fill selection with slice(None)
                for _ in range(self.ndim - len(index) + 1):
                    selection.append(slice(None))
            elif isinstance(item, slice):
                if item == slice(None):
                    selection.append(item)
                    continue
                if item.start is not None and item.start > self.shape[key]:
                    raise IndexError(
                        "index %d is out of bounds for axis %d with size %d"
                        % (item.start, key, self.shape[key]))
                else:
                    ix = item.indices(self.shape[key])
                    selection.append(range(ix[0], ix[1], ix[2]))
            else:
                if isinstance(item, int):
                    item = [item]
                    remove_dim.append(key)
                for k, i in enumerate(item):
                    if (i < -int(self.shape[key])) or (i >= self.shape[key]):
                        raise IndexError("index %d is out of bounds for axis "
                                         "%d with size %d"
                                         % (i, key, self.shape[key]))
                    # make sure to only keep positive indices
                    if i < 0:
                        item[k] = self.shape[key] + i
                selection.append(item)

        target_sz = list(self.shape)

        # to follow the python convention, single indices
        # will reduce the output's dimensionality
        for key, item in enumerate(selection):
            if item != slice(None):
                target_sz[key] = len(item)

        out = np.zeros(target_sz, dtype='complex64')
        # make sure that cha & col dim exist
        if self.flags['squeeze_ave_dims']:
            if self.flags['average']['Cha']:
                out = out[..., np.newaxis]
            if self.flags['average']['Col']:
                out = out[..., np.newaxis]
            if out.ndim < 3:
                out = out[np.newaxis]

        # 'vectorize' the output array for now
        out = out.reshape([-1, out.shape[-2], out.shape[-1]])

        # average counter to scale the data properly later
        ave_counter = np.zeros(np.prod(out.shape[:-2]), dtype=np.uint16)

        # now that we have our selection, we can read the data
        # for this, we simply go through all mdb's and fill them in if selected
        # this is not very efficient for large files, but fool-proof
        for mdb in self.mdb_list:

            counters = mdb.mdh['sLC'][self.sorted_mdh_keys]\
                .astype(self.dt_counters)

            if self.flags['skip_empty_lead']:
                lpos, ppos = self.dim_order.index('Lin'),\
                    self.dim_order.index('Par')
                counters[lpos] -= self._first_ix[lpos]
                counters[ppos] -= self._first_ix[ppos]

            # check if we have to read this mdb
            do_not_read = False
            for key, sel in enumerate(selection):
                if sel == slice(None):
                    # all data selected, no counter check required for this dim
                    continue
                if key >= self.ndim-2:
                    # skip col & cha
                    continue
                if self.flags['average'][self.dims[key]]:
                    # averaged dims are completely read
                    continue
                if counters[key] not in sel:
                    do_not_read = True
                    break

            if do_not_read:
                # go to next mdb
                continue

            # read data
            data = mdb.data

            # average cha if requested
            if self.flags['average']['Cha']:
                data = data.mean(-2, keepdims=True)

            # average col or apply oversampling removal if requested
            if self.flags['average']['Col']:
                data = data.mean(-1, keepdims=True)
            else:
                if self.flags['regrid'] and self.rs_traj is not None\
                        and not mdb.is_flag_set('SKIP_REGRIDDING'):
                    data = perform_regrid(
                        data, self.rs_traj, mdb.mdh["fReadOutOffcentre"])

                if self.flags['remove_os']:
                    data, _ = remove_oversampling(data)

            # reflect data if mdh flag is set
            if mdb.is_flag_set('REFLECT'):
                data = data[..., ::-1]

            ix = int(0)
            requests = 1
            for key, dim in enumerate(self.dims[:-2]):
                if self.flags['average'][dim]:
                    pass  # nothing to add
                elif key >= len(selection) or selection[key] == slice(None):
                    ix += int(counters[dim] * np.prod(target_sz[key+1:-2]))
                else:
                    ix += int(selection[key].index(counters[dim])
                              * np.prod(target_sz[key+1:-2]))

            # only keep selected channels & columns
            if len(selection) > self.ndim-2:
                # select channels
                data = data[selection[-2]]
            if len(selection) > self.ndim-1:
                # select columns
                data = data[:, selection[-1]]

            out[ix] += requests * data

            # increment average counter for ix
            ave_counter[ix] += 1

        # scale data according to ave_counter:
        ave_counter = np.maximum(ave_counter, 1)
        out /= ave_counter[..., np.newaxis, np.newaxis]

        # to follow the numpy convention,
        # single indices will reduce the output's dimensionality
        target_sz = [target_sz[key] for key in range(len(target_sz))
                     if key not in remove_dim]

        return out.reshape(target_sz)
