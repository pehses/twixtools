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
            sz_new = os.path.getsize(out_dat.name)
        
            all_close = True
            for mdb_orig, mdb_new in zip(read_twix(infile)[-1]['mdb'], read_twix(out_dat.name)[-1]['mdb']):
                if mdb_orig.is_flag_set('ACQEND') or mdb_orig.is_flag_set('SYNCDATA'):
                    continue
                if not np.allclose(mdb_orig.data, mdb_new.data, atol=zfp_tol):
                    all_close = False
                    break

        self.assertEqual(sz_orig, sz_new, 'zfp: file size not equal to original')
        self.assertTrue(all_close, 'zfp: not np.allclose()')


class remove_os_Test(unittest.TestCase):

    def test(self):
        infile = 'test/singlechannel.dat'
        sz_orig = os.path.getsize(infile)

        with tempfile.NamedTemporaryFile(suffix='.h5') as out_h5, tempfile.NamedTemporaryFile(suffix='.dat') as out_dat:
            twixzip.compress_twix(infile=infile, outfile=out_h5.name, remove_os=True)
            twixzip.reconstruct_twix(infile=out_h5.name, outfile=out_dat.name)
            sz_new = os.path.getsize(out_dat.name)
        
        self.assertEqual(sz_orig, sz_new, 'remove_os: file size not equal to original')
