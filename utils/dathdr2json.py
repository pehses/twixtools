#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import twixtools
import json

parser = argparse.ArgumentParser(description="Output the header of a .dat file")
parser.add_argument("datfile", type=str, help="Path to .dat file")
args = parser.parse_args()

twix = twixtools.read_twix(args.datfile, parse_data=False, verbose=False)
print(json.dumps(twix[-1]['hdr'], indent=4))
