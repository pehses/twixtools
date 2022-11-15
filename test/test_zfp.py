import os
import numpy as np
import tempfile
import unittest
from twixtools import twixzip, read_twix
from twixtools.twixzip import suppress_stdout_stderr


infile = 'example_data/gre.dat'


class test_zfp(unittest.TestCase):

    def test(self):
        sz_orig = os.path.getsize(infile)

        zfp_tol = 1e-5
        with tempfile.NamedTemporaryFile(suffix='.dat') as out_dat:
            with tempfile.NamedTemporaryFile(suffix='.h5') as out_h5:
                twixzip.compress_twix(infile=infile, outfile=out_h5.name,
                                      zfp=True, zfp_tol=zfp_tol)
                twixzip.reconstruct_twix(infile=out_h5.name,
                                         outfile=out_dat.name)

            self.assertEqual(sz_orig, os.path.getsize(out_dat.name),
                             'zfp: file size not equal to original')

            with suppress_stdout_stderr():
                twix_orig = read_twix(infile, keep_syncdata_and_acqend=True)[-1]
                twix_new = read_twix(out_dat.name, keep_syncdata_and_acqend=True)[-1]

            self.assertTrue(
                (np.all(twix_orig['hdr_str'] == twix_new['hdr_str'])),
                'zfp: headers do not match')

            for mdb_orig, mdb_new in zip(twix_orig['mdb'], twix_new['mdb']):
                if mdb_orig.is_flag_set('ACQEND'):
                    continue
                elif mdb_orig.is_flag_set('SYNCDATA'):
                    continue

                self.assertTrue(mdb_orig.mdh == mdb_new.mdh,
                                'zfp: mdhs do not match')
                self.assertTrue(
                    np.allclose(mdb_orig.data, mdb_new.data, atol=zfp_tol),
                    'zfp: mdb data not within zfp tolerance')


if __name__ == '__main__':
    unittest.main()
