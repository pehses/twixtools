import copy
from twixtools.mdh_def import is_flag_set
import numpy as np
from numpy.core.fromnumeric import shape
from scipy.integrate import cumtrapz
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
    # creates a list of measurements (or a single dict if input was dict)
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
            continue  # first "meas" may store first 10240 bytes of file
        
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
            out[-1][category] = twix_array(out[-1][category], meas['hdr'].copy())

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

    def __init__(self, mdb_list, hdr=None, flags=None):
        
        self.mdb_list = mdb_list.copy()
        self.hdr = None
        if hdr is not None:
            self.hdr = copy.deepcopy(hdr)

        self.rs_traj = calc_regrid_traj(self.hdr)
        
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

        self.base_size = np.ones(1, dtype=self.dt_dims)[0]
        for key, item in enumerate(shp): # complicated, can we do this converston better (proper casting?)
            self.base_size[key]=item 
        # todo: coil-compression ('cc', 'ncc')
        #       'skip_missing_lines': False, 'cc': -1}
        self._flags = {'average': {item:False for item in self.dims}, 'remove_os': False, 'regrid': False, 'zf_missing_lines': False}
        
        # averages should be averaged by default:
        self._flags['average']['Ave'] = True

        # set flags that were passed in constructor call        
        if flags is not None:
            for key,item in flags.items():
                try:
                    self._flags[key] = item.copy()
                except:
                    self._flags[key] = item


    def copy(self):
        return self.__copy__()

    def __copy__(self):
        self._flags = self._flags.copy()
        return twix_array(self.mdb_list, self.hdr, self.flags)

    @property
    def dims(self):
        return list(self.dt_dims.names)

    @property
    def non_singleton_dims(self):
        return [dim for dim in self.dims if self.size[dim]>1]

    @property
    def ndim(self):
        return len(self.dims)

    @property
    def flags(self):
        # wip: although the flags dict itself is write-protected, its entries are currently not and can be overwritten by garbage!
        return self._flags

    @property
    def size(self):
        # self.size returns the shape of the data as a dtype with named elements for easier access
        sz = self.base_size.copy()
        if self.flags['remove_os']:
            sz[-1] //= 2
        if self.hdr is not None and self.flags['zf_missing_lines']:
            hdr_lin = self.hdr['MeasYaps']['sKSpace']['lPhaseEncodingLines']
            sz['Lin'] = max(sz['Lin'], hdr_lin)
            if self.hdr['MeasYaps']['sKSpace']['ucDimension'] > 2:
                hdr_par = self.hdr['MeasYaps']['sKSpace']['lPartitions']
                sz['Par'] = max(sz['Par'], hdr_par)
        for dim in range(len(sz)):
            if self.flags['average'][self.dims[dim]]:
                sz[dim] = 1
        return sz


    @property
    def shape(self):
        # self.shape is the more numpy compatible version of self.size by returning a tuple
        return self.size.item()


    def __getitem__(self, index):
        # implement array slicing here
        # returns numpy.ndarray
        if not isinstance(index, tuple):
            index = (index,) # make sure to pass along tuple
        if len(index) > self.ndim:
            raise IndexError("too many indices for array: array is %d-dimensional, but %d were indexed"%(self.ndim, len(index)))
        ellipsis_in_index = False
        selection = list()
        remove_dim = list()
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
                    item = [item]
                    remove_dim.append(key)
                for i in item:
                    if (i < -int(self.shape[key])) or (i>=self.shape[key]):
                        raise IndexError("index %d is out of bounds for axis %d with size %d"%(i, key, self.shape[key]))
                selection.append(item)
        
        target_sz = list(self.shape)
         
         # to follow the python convention, single indices will reduce the output's dimensionality
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
                if self.flags['regrid'] and self.rs_traj is not None and not mdb.is_flag_set('SKIP_REGRIDDING'):
                    # first correct for readout shifts
                    # the nco frequency is always scaled to the max.
                    # gradient amp and does account for ramp-sampling
                    ncol = mdb.mdh["ushSamplesInScan"]
                    ro_shift = mdb.mdh["fReadOutOffcentre"]
                    deltak = max(abs(np.diff(self.rs_traj)))
                    adcphase = deltak * ro_shift * np.arange(ncol)
                    fovphase = ro_shift * self.rs_traj
                    data *= np.exp(1j*2*np.pi*(adcphase - fovphase))                    
                    
                    # interpolate
                    x = np.linspace(self.rs_traj[0], self.rs_traj[-1], ncol)
                    data = np.asarray([np.interp(x, self.rs_traj, y) for y in data])

                if self.flags['remove_os']:
                    data,_ = remove_oversampling(data)

            # reflect data if mdh flag is set
            if mdb.is_flag_set('REFLECT'):
                data = data[...,::-1]

            ix = int(0)
            requests = 1
            for dim in range(self.ndim-2):
                if self.flags['average'][self.dims[dim]]:
                    pass  # nothing to add
                elif dim >= len(selection) or selection[dim] == slice(None):
                    ix += int(counters[dim] * np.prod(target_sz[dim+1:-2]))
                else:
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
            ave_counter[ix] += 1

        # scale data according to ave_counter:
        ave_counter = np.maximum(ave_counter, 1)
        out /= ave_counter[..., np.newaxis, np.newaxis]

        # to follow the numpy convention, single indices will reduce the output's dimensionality
        target_sz = [target_sz[key] for key in range(len(target_sz)) if key not in remove_dim]

        return out.reshape(target_sz)



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
        gr_adc[ixDn] = 1 - (time_adc[ixDn] - rampup_time - flattop_time) / rampdown_time
    elif regrid_mode == 4:
        gr_adc[ixUp] = np.sin(np.pi / 2 * time_adc[ixUp] / rampup_time)
        gr_adc[ixDn] = np.sin(np.pi / 2 * (1 + (time_adc[ixDn] - rampup_time - flattop_time) / rampdown_time))
    else:
        raise Exception('regridding mode unknown')

    # make sure that gr_adc is always positive (rs_traj needs to be strictly monotonic)
    gr_adc = np.maximum(gr_adc, 1e-4)
    rs_traj = (np.append(0, cumtrapz(gr_adc)) - ncol // 2) / np.sum(gr_adc)
    rs_traj -= np.mean(rs_traj[ncol//2-1 : ncol//2+1])
 
    # scale rs_traj by kmax (only works if all slices have same FoV!!!)
    kmax = prot['MeasYaps']['sKSpace']['lBaseResolution'] / prot['MeasYaps']['sSliceArray']['asSlice'][0]['dReadoutFOV']
    rs_traj *= kmax

    return rs_traj
