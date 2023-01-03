import unittest
import twixtools


infile = 'example_data/gre.dat'


class test_maptwix(unittest.TestCase):

    def test(self):

        twix = twixtools.map_twix(infile)
        _ = twix[-1]['image'][:]

        self.assertEqual(1, 1, 'should never happen')


if __name__ == '__main__':
    unittest.main()
