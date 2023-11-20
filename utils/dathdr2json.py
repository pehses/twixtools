#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import twixtools
import json

parser = argparse.ArgumentParser(description="Output the header of a .dat file")
parser.add_argument("datfile", type=str, help="Path to .dat file")
parser.add_argument("--outfile", type=str, help="Output file / prefix")
args = parser.parse_args()

twix = twixtools.read_twix(args.datfile, parse_data=False, verbose=False)

for i,twi in enumerate(twix):
    js = twix[-1]['hdr']
    raw_keys = []
    for key in js:
        if "_raw" in key:
            raw_keys.append(key)
    for k in raw_keys:
        del js[k]

    if args.outfile is None:
        print(json.dumps(twix[-1]['hdr'], indent=4))
    else:
        ofname = args.outfile
        if len(twix) > 1:
            ofname += f"_{i}"
        with open(ofname, "w") as f:
            print(json.dumps(twix[-1]['hdr'], indent=4), file=f)
