#!/usr/bin/env python


"""This script is used to get the midpoint of region"""

import re
import sys
from FEICP import open_file

try:
    bed_file = sys.argv[1]
except:
    sys.exit(__doc__)

mid_size = 200

with open_file(bed_file) as fh:
    for line in fh:
        if re.search('^#', line):
            continue
        line = line.strip()
        lineL = line.split('\t')
        feature_len = int(lineL[2]) - int(lineL[1])
        if feature_len > mid_size:
            mid_point = (int(lineL[1]) + int(lineL[2])) / 2
            lineL[1] = int(mid_point - mid_size / 2)
            lineL[2] = int(mid_point + mid_size / 2)
        print('\t'.join([str(i) for i in lineL]))
