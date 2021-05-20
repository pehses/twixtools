import copy
import numpy as np
from numpy.core.fromnumeric import shape
import twixtools
from twixtools.recon_helpers import remove_oversampling


# define categories in which the twix data should be sorted based on MDH flags that must or must not be set (True/False)
# only one 'any' per category allowed (however, it is possible to add other appropriate functions (even synonyms for any))
twix_category = dict()
twix_category['image']           = {'RTFEEDBACK': False, 'HPFEEDBACK': False, 'REFPHASESTABSCAN': False, 'PHASESTABSCAN': False,
                                    'PHASCOR': False, 'NOISEADJSCAN': False, 'noname60': False}
twix_category['noise']           = {'NOISEADJSCAN': True}
twix_category['phasecorr']       = {'PHASCOR': True, 'PATREFSCAN': False, 'noname60': False}
twix_category['phasestab']       = {'PHASESTABSCAN': True, 'REFPHASESTABSCAN': False, any: {'PATREFSCAN': False, 'PATREFANDIMASCAN': True}, 'noname60': False} 
twix_category['phasestab_ref0']  = {'REFPHASESTABSCAN': True, 'PHASESTABSCAN': False, any: {'PATREFSCAN': False, 'PATREFANDIMASCAN': True}, 'noname60': False} 
twix_category['phasestab_ref1']  = {'REFPHASESTABSCAN': True, 'PHASESTABSCAN': True, any: {'PATREFSCAN': False, 'PATREFANDIMASCAN': True}, 'noname60': False} 
twix_category['refscan']         = {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True}, 'PHASCOR': False, 'PHASESTABSCAN': False, 
                                   'REFPHASESTABSCAN': False, 'RTFEEDBACK': False, 'HPFEEDBACK': False, 'noname60': False}
twix_category['refscan_pc']      = {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True}, 'PHASCOR': True}
twix_category['refscan_ps']      = {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True}, 'REFPHASESTABSCAN': False, 'PHASESTABSCAN': True}
twix_category['refscan_ps_ref0'] = {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True}, 'REFPHASESTABSCAN': True, 'PHASESTABSCAN': False}
twix_category['refscan_ps_ref1'] = {any: {'PATREFSCAN': True, 'PATREFANDIMASCAN': True}, 'REFPHASESTABSCAN': True, 'PHASESTABSCAN': True}
twix_category['rt_feedback']     = {any: {'RTFEEDBACK': True, 'HPFEEDBACK': True}, 'MDH_VOP': False}
twix_category['vop']             = {'MDH_VOP': True}
twix_category['fidnav']          = {'noname60': True}  # this is the real reason for including the 'noname60' check above


def map_twix(input):
    # creates a list of measurements
    # with data for each measurement mapped to a twix_array object

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
            pass  # first "meas" may store first 10240 bytes of file
        
        # append new dict to output list
        out.append(dict())
        
        # sort mdbs into categories
        for mdb in meas['mdb']:
            if mdb.is_flag_set('SYNCDATA'): # ignore syncdata
                continue
            if mdb.is_flag_set('ACQEND'):
                break

            for category, requirements in twix_category.items():
                include_in_cat = True
                for flag in requirements.keys():
                    if isinstance(flag, str):   # simple check whether flag is set
                        if mdb.is_flag_set(flag) != requirements[flag]:
                            include_in_cat = False
                            break
                    else:  # assume that this is a function call (probably any())
                        checks = list()
                        for flag2 in requirements[flag].keys():
                            checks.append(mdb.is_flag_set(flag2)==requirements[flag][flag2])
                        if not flag(checks):
                            include_in_cat = False
                            break
                if include_in_cat:
                    if category not in out[-1]:
                        out[-1][category] = list()
                    out[-1][category].append(mdb)

        # convert each categories' mdb list to twix_array
        for category in out[-1].keys():
            out[-1][category] = twix_array(out[-1][category])

        # include hdr in dict
        out[-1]['hdr'] = meas['hdr'].copy()
        out[-1]['hdr_str'] = meas['hdr_str'].copy()

        # go back to dict if input was dict
        if isinstance(input, dict):
            out = out[0]

    return out


class twix_array():

    def __init__(self, mdb_list):
        
        self.mdb_list = mdb_list.copy()
        
        # delete 'ACQEND' and 'SYNCDATA' flags if present
        twixtools.del_from_mdb_list(self.mdb_list, lambda mdb: mdb.is_flag_set('ACQEND') or mdb.is_flag_set('SYNCDATA'))

        self.dim_order = ["Ide", "Idd", "Idc", "Idb", "Ida", "Seg", "Set", "Rep", "Phs", "Eco", "Par", "Sli", "Ave", "Lin", "Cha", "Col"]
        self.dt_dims = np.dtype([(name, "<u2") for name in self.dim_order])  # dtype that includes all dims
        self.dt_counters = np.dtype([(name, "<u2") for name in self.dim_order[:-2]])  # dtype that only includes counters (no cha & col)

        self.key_map = {'Ide': 'ushIde', 'Idd': 'ushIdd', 'Idc': 'ushIdc', 'Idb': 'ushIdb', 'Ida': 'ushIda',
        'Seg': 'ushSeg', 'Set': 'ushSet', 'Rep': 'ushRepetition', 'Phs': 'ushPhase', 'Eco': 'ushEcho',
        'Par': 'ushPartition', 'Sli': 'ushSlice', 'Ave': 'ushAcquisition', 'Lin': 'ushLine'}
        self.sorted_mdh_keys = [self.key_map[dim] for dim in self.dims[:-2]]
        
        # determine k-space shape by finding max index
        shp = np.ones(len(self.dt_dims), dtype=self.dt_dims[1])
        for mdb in self.mdb_list:

            sLC = mdb.mdh['sLC']
            sLC = np.asarray(sLC[self.sorted_mdh_keys].tolist(), dtype = sLC[0].dtype)
            req_shape = 1 + sLC
            # add channels & columns
            req_shape = np.concatenate([req_shape, [mdb.mdh['ushUsedChannels'], mdb.mdh['ushSamplesInScan']]])
            shp = np.maximum(shp, req_shape)

        self.base_shape = np.ones(1, dtype=self.dt_dims)[0]
        for key, item in enumerate(shp): # complicated, can we do this converston better (proper casting?)
            self.base_shape[key]=item 
        # todo: coil-compression ('cc', 'ncc')
        self._flags = {'average_dim': np.zeros(self.ndim, dtype=bool), 'remove_os': False, 'regrid': False, 'cc': False, 'ncc': -1}

    @property
    def dims(self):
        return list(self.dt_dims.names)

    @property
    def non_singleton_dims(self):
        return self.dims[self.shape>1]

    @property
    def ndim(self):
        return len(self.dims)

    @property
    def flags(self):
        # wip: although the flags dict itself is write-protected, its entries are currently not and can be overwritten by garbage!
        return self._flags

    @property
    def shape(self):
        shp = self.base_shape.copy()
        if self.flags['remove_os']:
            shp[-1] //= 2
        for dim in range(len(shp)):
            if self.flags['average_dim'][dim]:
                shp[dim] = 1
        return shp

    def __getitem__(self, index):
        # implement array slicing here
        # returns numpy.ndarray
        if not isinstance(index, tuple):
            index = (index,) # make sure to pass along tuple
        if len(index) > self.ndim:
            raise IndexError("too many indices for array: array is %d-dimensional, but %d were indexed"%(self.ndim, len(index)))
        ellipsis_in_index = False
        selection = list()
        for key, item in enumerate(index):
            if ellipsis_in_index:
                key += self.ndim - len(index)
            if item is Ellipsis:
                if ellipsis_in_index:
                    raise IndexError("an index can only have a single ellipsis ('...')")
                ellipsis_in_index = True
                # fill selection with slice(None)
                for k in range(self.ndim - len(index) + 1):
                    selection.append(slice(None))
            elif isinstance(item, slice):
                if item == slice(None):
                    selection.append(item)
                    continue
                if item.start is not None and item.start > self.shape[key]:
                    raise IndexError("index %d is out of bounds for axis %d with size %d"%(item.start, key, self.shape[key]))
                else:
                    ix = item.indices(self.shape[key])
                    selection.append(range(ix[0], ix[1], ix[2]))
            else:
                if isinstance(item, int):
                    item = [item] # create list
                for i in item:
                    if (i < -int(self.shape[key])) or (i>=self.shape[key]):
                        raise IndexError("index %d is out of bounds for axis %d with size %d"%(i, key, self.shape[key]))
                selection.append(item)
        
        target_sz = list(self.shape.copy().tolist())
        for key, item in enumerate(selection):
            if item != slice(None):
                target_sz[key] = len(item)

        out = np.zeros(target_sz, dtype='complex64')
        out = out.reshape([-1, out.shape[-2], out.shape[-1]]) # for now 'vectorize' it

        # average counter to scale the data properly later
        ave_counter = np.zeros(np.prod(out.shape[:-2]), dtype=np.uint16)

        # now that we have our selection and allocated memory, we can read in the data
        # for this, we simply go through all mdb's and fill them in if selected
        # this is not very efficient for large files, but fool-proof
        for mdb in self.mdb_list:

            counters = mdb.mdh['sLC'][self.sorted_mdh_keys].astype(self.dt_counters)
            
            # check if we have to read this mdb
            do_not_read = False
            for key, sel in enumerate(selection):
                if sel == slice(None):
                    # all data selected, no counter check required for this dim
                    continue
                if key>=self.ndim-2:
                    # skip col & cha
                    continue
                if self.flags['average_dim'][key]:
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
            if self.flags['average_dim'][-2]:
                data = data.mean(-2, keepdims=True)

            # average col or apply oversampling removal if requested
            if self.flags['average_dim'][-1]:
                data = data.mean(-1, keepdims=True)
            else:
                if self.flags['regrid']:
                    pass  # wip
                if self.flags['remove_os']:  
                    data,_ = remove_oversampling(data)

            ix = int(0)
            requests = 1
            for dim in range(self.ndim-2):
                if self.flags['average_dim'][dim]:
                    pass  # nothing to add
                elif dim >= len(selection) or selection[dim] == slice(None):
                    ix += int(counters[dim] * np.prod(target_sz[dim+1:-2]))
                else:
                    # wip
                    # test if this is working correctly + need to handle repeated indices (requests)
                    print(selection[dim], counters[dim], selection[dim].index(counters[dim]))
                    ix += int(selection[dim].index(counters[dim]) * np.prod(target_sz[dim+1:-2]))
            
            # only keep selected channels & columns
            if len(selection) > self.ndim-2:
                # select channels
                data = data[selection[-2]]
            if len(selection) > self.ndim-1:
                # select columns
                data = data[:,selection[-1]]

            out[ix] += requests * data

            # increment average counter for ix
            ave_counter[ix] += ix

        # scale data according to ave_counter:
        ave_counter = np.maximum(ave_counter, 1)
        out /= ave_counter[..., np.newaxis, np.newaxis]

        return out.reshape(target_sz)
