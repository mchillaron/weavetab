#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of pypistrello.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
#

from pathlib import Path

import argparse
import os
import re

def extract_cube_spectra(cube_file, output_dir):
    pass


def main():
    parser = argparse.ArgumentParser(description="Extract flux-calibrated spectra and their spatial coordinates from WEAVE FITS data cubes.")
    parser.add_argument('-F', '--cube-file', type=str, help="Name of the WEAVE FITS data cube file.")
    parser.add_argument('-n', '--number-simulations', type=int, default=0, help="Number of simulations to perform.")
    args = parser.parse_args()

    cube_file = args.cube_file
    cube_path = Path('.') / cube_file
    # get the absolute path of the cube file
    cube_path = cube_path.resolve()
    print(f"Processing cube file: {cube_path}")
    number_simulations = args.number_simulations

    extract_cube_spectra(cube_path, number_simulations)


if __name__ == "__main__":
    main()