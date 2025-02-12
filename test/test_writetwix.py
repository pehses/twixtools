import unittest
import twixtools
import os
import hashlib
import tempfile

infile = 'example_data/gre.dat'


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.digest()


class test_writetwix(unittest.TestCase):

    def test(self):
        orig_size = os.path.getsize(infile)
        orig_md5 = md5(infile)
        twix = twixtools.read_twix(infile, keep_acqend=True)

        try:
            # Create a NamedTemporaryFile but close and delete it immediately
            with tempfile.NamedTemporaryFile(suffix='.dat', delete=True) as temp_file:
                out_name = temp_file.name

            twixtools.write_twix(twix, out_name)
            saved_size = os.path.getsize(out_name)
            saved_md5 = md5(out_name)
            self.assertEqual(orig_size, saved_size, "File size does not match after writing to file.")
            self.assertEqual(orig_md5, saved_md5, "MD5 checksums do not match after writing to file.")
        finally:
            os.remove(out_name)


if __name__ == '__main__':
    unittest.main()
