#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import twixtools
import json
import sys


def main():
    parser = argparse.ArgumentParser(description="Output the header of a .dat file")
    parser.add_argument("datfile", type=str, help="Path to .dat file")
    parser.add_argument("--out", type=str, help="Output file / prefix instead of stdout")
    parser.add_argument("--keys", type=str, help="Comma-separated list of hdrs to obtain; h to list available keys")
    args = parser.parse_args()

    twix = twixtools.read_twix(args.datfile, parse_data=False, verbose=False)

    if args.keys:
        selected_keys = [k.strip() for k in args.keys.split(',')]

    if args.keys and 'h' == selected_keys[0]:
        print(f"Allowed key types: {','.join([k for k in twix[-1]['hdr'].keys() if not '_raw' in k])}")
        sys.exit(0)

    for i, twix_i in enumerate(twix):

        if args.keys:
            for key in selected_keys:
                if not key in twix_i['hdr'].keys():
                    print(f"Requested key {key} not found in twix header no. {i}", file=sys.stderr)
                    sys.exit(1)

        del_keys = []

        for key in twix_i['hdr']:
            if "_raw" in key:
                del_keys.append(key)
            elif args.keys and not key in selected_keys:
                del_keys.append(key)

        for k in del_keys:
            del twix_i['hdr'][k]

        if args.out is None:
            print(json.dumps(twix_i['hdr'], indent=4))
        else:
            ofname = args.out
            if len(twix) > 1:
                ofname += f"_{i}"
            with open(ofname, "w") as f:
                print(json.dumps(twix_i['hdr'], indent=4), file=f)


if __name__ == "__main__":
    main()
