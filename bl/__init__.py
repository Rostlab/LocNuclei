# -*- coding: utf8 -*-
from __future__ import print_function
import argparse


def main():
    # parse command line options
    parser = argparse.ArgumentParser()
    parser.add_argument('param', help='parm-help')
    parser.add_argument('-v', '--verbose', help='Toggles verbose mode on', action='store_true')
    args = parser.parse_args()


def error(*objs):
    print("ERROR: ", *objs, file=sys.stderr)


if __name__ == "__main__":
    main()