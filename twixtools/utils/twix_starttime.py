#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import twixtools
import datetime

def main():
    parser = argparse.ArgumentParser(description="Output the acquisition start time of the first scan in the datfile in iso format.")
    parser.add_argument("datfile", type=str, help="Path to .dat file")
    args = parser.parse_args()

    twix = twixtools.read_twix(args.datfile, parse_data=True, verbose=False, parse_pmu=False)

    # timestamp is in 2.5ms since midnight.
    ts = twix[0]['mdb'][0].mdh.TimeStamp
    seconds = ts // 400
    us = (ts - 400 * seconds) * 2500

    t = datetime.datetime.fromtimestamp(seconds).time()
    t.replace(microsecond = us)

    print(t.isoformat(timespec="milliseconds"))

if __name__ == "__main__":
    main()
