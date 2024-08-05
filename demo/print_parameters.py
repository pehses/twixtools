import twixtools

dat_file = 'xxx.dat'

twixobj = twixtools.read_twix(dat_file)

hdr_twix = twixobj[-1]['hdr']

MB = hdr_twix['Phoenix']['sSliceAcceleration']['lMultiBandFactor']         # multi-band
N_slices = hdr_twix['Phoenix']['sSliceArray']['lSize']                     # slices
N_segments = hdr_twix['Phoenix']['sFastImaging']['lSegments']              # segments
N_EchoTrainLength = int(hdr_twix['Meas']['EchoTrainLength'])               # echo train length
N_Accel_PE = hdr_twix['Phoenix']['sPat']['lAccelFactPE']                   # phase-encoding acceleration

dwelltime_in_s = 1e-9 * hdr_twix['MeasYaps']['sRXSPEC']['alDwellTime'][0]  # dwell time
base_res = hdr_twix['MeasYaps']['sKSpace']['lBaseResolution']              # base resolution
ADC_dur = dwelltime_in_s * base_res * 2                                    # ADC duration (assume oversampling factor of 2)
bandwidth_per_pixel = 1 / ADC_dur                                          # bandwidth per pixel