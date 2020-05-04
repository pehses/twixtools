from twixtools import twixzip
import tempfile
import hashlib
import os
import unittest


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


class ZFP_OSrm_Test(unittest.TestCase):

    def test(self):
        infile = 'test/singlechannel.dat'
        sz_orig = os.path.getsize(infile)

        with tempfile.NamedTemporaryFile(suffix='.h5') as out_h5, tempfile.NamedTemporaryFile(suffix='.dat') as out_dat:
            twixzip.compress_twix(infile=infile, outfile=out_h5.name, zfp=True, remove_os=True)
            twixzip.reconstruct_twix(infile=out_h5.name, outfile=out_dat.name)
            sz_new = os.path.getsize(out_dat.name)

        self.assertEqual(sz_orig, sz_new, 'zfp & removeOS: file size not equal to original')
