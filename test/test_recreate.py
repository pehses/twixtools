import os
import numpy as np
import tempfile
import unittest
from twixtools import twixzip, read_twix
from twixtools.mdh_def import is_flag_set


def md5(fname):
    import hashlib
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.digest()


class Lossless_Test(unittest.TestCase):

    def test(self):
        infile = 'test/singlechannel.dat'
        md5_orig = md5(infile)

        with tempfile.NamedTemporaryFile(suffix='.h5') as out_h5, tempfile.NamedTemporaryFile(suffix='.dat') as out_dat:
            twixzip.compress_twix(infile=infile, outfile=out_h5.name)
            twixzip.reconstruct_twix(infile=out_h5.name, outfile=out_dat.name)
            md5_new = md5(out_dat.name)

        self.assertEqual(md5_orig, md5_new, 'lossless compression: md5 hash does not match with original')


class ZFP_Test(unittest.TestCase):

    def test(self):
        infile = 'test/singlechannel.dat'
        sz_orig = os.path.getsize(infile)

        zfp_tol = 1e-5
        with tempfile.NamedTemporaryFile(suffix='.h5') as out_h5, tempfile.NamedTemporaryFile(suffix='.dat') as out_dat:
            twixzip.compress_twix(infile=infile, outfile=out_h5.name, zfp=True, zfp_tol=zfp_tol)
            twixzip.reconstruct_twix(infile=out_h5.name, outfile=out_dat.name)
            
            self.assertEqual(sz_orig, os.path.getsize(out_dat.name), 'zfp: file size not equal to original')

            twix_orig = read_twix(infile)[-1]
            twix_new = read_twix(out_dat.name)[-1]
             
            self.assertTrue((np.all(twix_orig['hdr_str']==twix_new['hdr_str'])), 'zfp: headers do not match')

            for mdb_orig, mdb_new in zip(twix_orig['mdb'], twix_new['mdb']):
                if mdb_orig.is_flag_set('ACQEND'):
                    continue
                elif mdb_orig.is_flag_set('SYNCDATA'):
                    continue
                
                self.assertTrue(mdb_orig.mdh == mdb_new.mdh, 'zfp: mdhs do not match')
                self.assertTrue(np.allclose(mdb_orig.data, mdb_new.data, atol=zfp_tol), 'zfp: mdb data not within zfp tolerance')


class remove_os_Test(unittest.TestCase):

    def test(self):
        infile = 'test/singlechannel.dat'
        sz_orig = os.path.getsize(infile)

        with tempfile.NamedTemporaryFile(suffix='.h5') as out_h5, tempfile.NamedTemporaryFile(suffix='.dat') as out_dat:
            twixzip.compress_twix(infile=infile, outfile=out_h5.name, remove_os=True)
            twixzip.reconstruct_twix(infile=out_h5.name, outfile=out_dat.name)
            sz_new = os.path.getsize(out_dat.name)
        
        self.assertEqual(sz_orig, sz_new, 'remove_os: file size not equal to original')
