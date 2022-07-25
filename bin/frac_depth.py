#!/usr/bin/env python

__doc__ = """Compute fraction of bases for MN908947.3 at or above a desired depth."""

import sys, csv


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
