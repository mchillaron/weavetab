# weavetab
![Python](https://img.shields.io/badge/Python-3.9+-3776AB?logo=python&logoColor=white)

**weavetab** is a Python command-line tool designed to extract flux-calibrated spectra and their spatial coordinates from **WEAVE FITS data cubes**. The tool can convert original cubes into tabular spectral 
data or new flux-calibrated cubes, and can generate simulated cubes too.

It is intended for working with astronomical data products from the **WEAVE** spectrograph (WHT), allowing users to easily manipulate spectral cubes and perform simulations.

---

# Features

- Extract spectra from **WEAVE** FITS data cubes and calibrate them in flux.
- Save the wavelength axis to a CSV file.
- Convert cubes into FITS tables or new cubes with:
  - **Flux-calibrated** spectra
  - Inverse variance converted to **sigma**
- Generate **simulated cubes** based on the extracted data.
- Optionally process only a selected **region of interest** (x,y) instead of the whole cube.

---

# Installation

To install the package in your environment, first clone the repository from GitHub by running in your terminal:

```bash
git clone https://github.com/mchillaron/weavetab.git
```

Then, navigate into the megaradrpsimul/ folder:
```bash
cd megaradrpsimul/
```
And install the package in editable mode:
```bash
pip install -e .
```
The code will be now ready to use!

---

# Command Line Options

| Argument | Description |
|--------|-------------|
| `-F`, `--cube-file` | Name of the WEAVE FITS data cube file (required). |
| `-reg`, `--region` | Region of interest in FITS pixel coordinates: `x1 x2 y1 y2`. |
| `-n`, `--number-simulations` | Number of simulated cubes to generate (default: `0`). |
| `-s`, `--save-table` | Save a FITS table containing spatial coordinates and spectra per pixel (default: `False`). |

---

# Examples

### Extract spectra from a full cube

```bash
weavetab -F stackcube.fits
```

### Extract spectra from a specific region

```bash
weavetab -F stackcube.fits -reg 100 200 100 200
```

### Generate simulated cubes

```bash
weavetab -F stackcube.fits -n 5
```

This will produce:

- 1 reconstructed cube
- 5 simulated cubes

### Save spectra as a FITS table

```bash
weavetab -F stackcube.fits -s True
```

---

# Output

The program creates an output directory named:

```
<cube_filename>_spectra_cubes
```

Inside this directory you will find:

- `wavelength_<arm>.csv` — wavelength axis
- reconstructed FITS cube
- simulated cubes directory (if requested)
- optional FITS table with spectra (if requested)

If the directory already exists, the program will warn before overwriting files.

---

# Author

Maintainer: Maria Chillaron <mchill01@ucm.es>