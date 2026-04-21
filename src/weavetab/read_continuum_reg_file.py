#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of weavetab.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

import numpy as np

def read_continuum_file(filename, redshift):
    regions = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) == 2:
                lmin, lmax = map(float, parts)
                regions.append((lmin, lmax))

    regions_restframe = np.array(regions)
    regions_obs = (1+redshift) * regions_restframe
    return regions_obs