#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def read_continuum_file(filename):
    regions = []
    x = y = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Detect coordinates
            if line.startswith('(') and line.endswith(')'):
                coords = line.strip('()').split(',')
                x = int(coords[0])
                y = int(coords[1])
                print(f"[INFO] Using spaxel (x={x}, y={y})")
            
            else:
                parts = line.split()
                if len(parts) == 2:
                    lmin, lmax = map(float, parts)
                    regions.append((lmin, lmax))

    return np.array(regions), x, y