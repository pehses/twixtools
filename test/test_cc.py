import os
import numpy as np
import tempfile
import unittest
from twixtools import twixzip, read_twix
from twixtools.mdh_def import is_flag_set
from twixtools.twixzip import suppress_stdout_stderr


infile = 'example_data/gre.dat'

class test_scc(unittest.TestCase):

    def test(self):

        with suppress_stdout_stderr():
            twix_orig = read_twix(infile)[-1]
        
        nc = twix_orig['mdb'][1].mdh['ushUsedChannels']

        with tempfile.NamedTemporaryFile(suffix='.dat') as out_python:
            with tempfile.NamedTemporaryFile(suffix='.h5') as out_h5:
                twixzip.compress_twix(infile=infile, outfile=out_h5.name, cc_mode='scc', ncc=nc)
                twixzip.reconstruct_twix(infile=out_h5.name, outfile=out_python.name)
            
            with suppress_stdout_stderr():
                twix_orig = read_twix(infile)[-1]
                twix_python = read_twix(out_python.name)[-1]

            for mdb_orig, mdb_python in zip(twix_orig['mdb'], twix_python['mdb']):
                if mdb_orig.is_flag_set('ACQEND'):
                    continue
                elif mdb_orig.is_flag_set('SYNCDATA'):
                    continue

                self.assertTrue(mdb_orig.mdh == mdb_python.mdh, 'scc: mdhs do not match')
                self.assertTrue(np.allclose(mdb_orig.data, mdb_python.data), 'scc: mdb data not within tolerance')


class test_gcc(unittest.TestCase):

    def test(self):

        with suppress_stdout_stderr():
            twix_orig = read_twix(infile)[-1]
        
        nc = twix_orig['mdb'][1].mdh['ushUsedChannels']

        with tempfile.NamedTemporaryFile(suffix='.dat') as out_python:
            with tempfile.NamedTemporaryFile(suffix='.h5') as out_h5:
                twixzip.compress_twix(infile=infile, outfile=out_h5.name, cc_mode='gcc', ncc=nc)
                twixzip.reconstruct_twix(infile=out_h5.name, outfile=out_python.name)
            
            with suppress_stdout_stderr():
                twix_orig = read_twix(infile)[-1]
                twix_python = read_twix(out_python.name)[-1]

            for mdb_orig, mdb_python in zip(twix_orig['mdb'], twix_python['mdb']):
                if mdb_orig.is_flag_set('ACQEND'):
                    continue
                elif mdb_orig.is_flag_set('SYNCDATA'):
                    continue

                self.assertTrue(mdb_orig.mdh == mdb_python.mdh, 'gcc: mdhs do not match')
                self.assertTrue(np.allclose(mdb_orig.data, mdb_python.data), 'gcc: mdb data not within tolerance')



if __name__ == '__main__':
    unittest.main()
