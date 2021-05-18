import twixtools

def map_twix(input):
    # creates a list of measurements
    # with data for each measurement mapped to a twix_array object

    return_list = True
    if isinstance(input, list):
        # assume list of measurements
        twix = input
    elif isinstance(input, dict):
        # assume measurement dict
        twix = [input]       # put dict into list for code compatibility
        return_list = False  # the user probably does not expect a list, return twix_array instead
    else:
        # assume that this is the filename or a meas id
        twix = twixtools.read_twix(input)

    out = list()
    for meas in twix:
        if not isinstance(meas, dict):
            pass  # first "meas" may store first 10240 bytes of file
        out.append(twix_array(meas))
    
    if not return_list:
        out = out[-1]

    return out


class twix_array(dict):
    # WIP: not ready yet
    def __init__(self, mdb_list):
        super(twix_array, self).__init__({key:list() for key in ['ACQEND', 'RTFEEDBACK', 'HPFEEDBACK', 'SYNCDATA', 'REFPHASESTABSCAN', 'PHASESTABSCAN', 'PHASCOR', 'NOISEADJSCAN', 'noname60', 'PATREFSCAN', 'IMASCAN']})
        self.mdb_list = copy.deepcopy(mdb_list)
        self.sLC = dict()
        self.mdb_shape = dict()
        self.parse()
        

    def parse(self):
        self.sLC = {key:list() for key in self.keys()}
        self.mdb_shape = {key:list() for key in self.keys()}
        for mdb in self.mdb_list:
            if mdb.is_flag_set('ACQEND') or mdb.is_flag_set('SYNCDATA'):
                continue
            for cat in list(self.keys())[:-1]:
                if mdb.is_flag_set(cat):
                    self.get(cat).append(mdb)
                    self.sLC[cat].append(mdb.mdh['sLC'])
                    self.mdb_shape[cat].append([mdb.mdh['ushUsedChannels'], mdb.mdh['ushSamplesInScan']])
            if mdb.is_image_scan():
                self['IMASCAN'].append(mdb)
                self.sLC['IMASCAN'].append(mdb.mdh['sLC'])
                self.mdb_shape['IMASCAN'].append([mdb.mdh['ushUsedChannels'], mdb.mdh['ushSamplesInScan']])
        
        for cat in list(self.keys()):
            if len(self.get(cat))==0:
                # remove empty categories
                del(self[cat])
                del(self.sLC[cat])
                del(self.mdb_shape[cat])
            else:
                # add category to self
                pass


    # def __getitem__(self, index):
    #     pass
