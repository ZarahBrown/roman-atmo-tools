# Roman CGI Bandpass Tools

This repository provides tools to compute band-averaged fluxes using the Nancy Grace Roman Space Telescope Coronagraph Instrument (CGI) filters, and to convert flux densities between frequency-based units (Fν) and wavelength-based units (Fλ).

---

## Contents

- `cgibands.py`: Compute mean flux within CGI-defined or user-defined wavelength bands.
- `fnu2flam.py`: Convert flux values from Fν (e.g., Jy or erg/cm²/s/Hz) to Fλ (e.g., W/m²/μm).

---

## `cgibands.py` Features

- Use official Roman CGI transmission curves to compute weighted average fluxes.
- Optionally provide your own custom wavelength bands.
- Plot the spectrum and optionally overlay CGI bandpass curves.
- Support for computing planet-to-star flux ratios (Fp/Fs) if stellar spectrum and radii are provided.

### Inputs

- `wavelength`, `flux`: Your spectrum (wavelengths + fluxes).
- `flux_unit`, `wavelength_unit`: Specify units (e.g., 'W/m2/um', 'erg/cm2/s/A', etc.).
- Optional: `stellar_flux`, `stellar_wavelength`, `stellar_radius`, `planet_radius`.
- Optional: Custom band definitions as a Python dictionary.
- Flags: `plot`, `verbose`, `show_color_filter_curves`, `log_y`.

---

## `fnu2flam.py` Features

- Convert flux from frequency-based units to wavelength-based units.
- Supports units like `Jy`, `W/m²/Hz`, `erg/cm²/s/Hz`, etc.
- Output units can be specified (e.g., `W/m²/μm`, `erg/cm²/s/Å`).

### Example Usage

```python
from fnu2flam import convert_fnu_to_flambda

wavelength_um = np.linspace(0.5, 2.5, 100)
flux_jy = np.ones_like(wavelength_um) * 1.0

w_out, f_out = convert_fnu_to_flambda(wavelength_um, flux_jy,
                                      flux_unit='Jy',
                                      wavelength_unit='um',
                                      output_unit='W / (m2 um)')
```

---

## Filter Data

To use the official CGI filter curves:
Download from [roman.ipac.caltech.edu](https://roman.ipac.caltech.edu/page/additional-coronagraph-instrument-parameters-model-and-data-html#Color_Filter_Curves)

---

## Demo Notebook

A Jupyter notebook (`cgibands_demo.ipynb`) is available with examples:
- Using built-in CGI filters
- Using user-defined bands
- Reading in stellar and planetary spectra
- Converting flux units

---

## Requirements

- Python ≥ 3.8
- `numpy`, `matplotlib`, `astropy`, `pandas`

---

## License

MIT License. See `LICENSE.txt`.

---

## Author

Zarah Brown – Lunar and Planetary Laboratory, University of Arizona
