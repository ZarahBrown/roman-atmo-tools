import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.units import Unit
import os
import warnings
from astropy.utils.exceptions import AstropyWarning
from astropy.constants import R_jup, R_sun

DEFAULT_FILTER_DIR = 'C:/Research/Projects/Roman/Data/Filters' # Replace this with your directory for the CGI color contrast files

warnings.filterwarnings('ignore', category=AstropyWarning)

# ----------------------
# CGI Band Filter Utilities
# ----------------------
def read_transmission_csv(filepath):
    # Read the Color Filter Curve CSV files
    data = pd.read_csv(filepath, comment='#', encoding='latin1')
    data.columns = ['wavelength_nm', 'transmission_percent']
    data['transmission'] = data['transmission_percent'] / 100.0
    return data[['wavelength_nm', 'transmission']]

def get_transmission_bands(folder=DEFAULT_FILTER_DIR):
    files = {
        'B1': 'transmission_ID-01_1F_v0.csv',
        'B2': 'transmission_ID-02_2F_v0.csv',
        'B3': 'transmission_ID-03_3F_v0.csv',
        'B4': 'transmission_ID-04_4F_v0.csv'
    }

    bands = {}
    for band, fname in files.items():
        path = os.path.join(folder, fname)
        if not os.path.exists(path):
            print(f"Warning: Missing file for {band}: {fname}")
            continue
        bands[band] = read_transmission_csv(path)

    return bands

def fwhm_clip(wavelengths, transmission):
    half_max = 0.5 * transmission.max()
    return transmission > half_max

# ----------------------
# Ensure Flux_lamda
# ----------------------

def is_fnu_unit(unit_string):
    """
    Determine if the given unit string is in frequency units (Fν).
    """
    try:
        unit = Unit(unit_string)
        return unit.is_equivalent(u.Jy) or unit.is_equivalent(u.erg / (u.cm**2 * u.s * u.Hz))
    except Exception:
        return False

# ----------------------
# Main compute function
# ----------------------
def compute_cgibands(wavelength, flux_array, flux_unit,
                     filter_folder=None, plot=False,
                     wavelength_unit='nm', custom_bands=None,
                     show_color_filter_curves=False, log_y=False,
                     stellar_wavelength=None, stellar_flux=None,
                     stellar_radius=None,  # in R_sun
                     planet_radius=None,   # in R_jup
                     stellar_wavelength_unit=None, stellar_flux_unit=None):


    wavelength = np.array(wavelength)
    flux_array = np.array(flux_array)

    # Convert wavelength to nm
    wavelength_nm = (wavelength * Unit(wavelength_unit)).to(u.nm).value
    if stellar_wavelength is not None and stellar_wavelength_unit is not None:
        stellar_wavelength_nm = (np.array(stellar_wavelength) * Unit(stellar_wavelength_unit)).to('nm').value
    else:
        stellar_wavelength_nm = None

    # Check that the flux is NOT F_nu
    if is_fnu_unit(flux_unit):
        raise ValueError(f"The unit '{flux_unit}' appears to be an F_nu unit (per frequency). "
                         f"Please convert to F_lambda (per wavelength) before using this function.")
    if stellar_flux is not None:
        if is_fnu_unit(stellar_flux_unit):
            raise ValueError("Stellar flux must be in F_lambda units (e.g., erg/cm2/s/A, W/m2/nm). F_nu units not yet supported.")
        else:
            stellar_flux = np.array(stellar_flux)
    
    flux_unit_label = f"Planet Flux [{flux_unit}]"
    fp_over_fs = None  # default

    if custom_bands is None:
        bands = get_transmission_bands(folder=filter_folder or DEFAULT_FILTER_DIR)
    else:
        # Create flat-top filters for custom bands
        bands = {}
        for k, (start, end) in custom_bands.items():
            trans = np.zeros_like(wavelength_nm)
            mask = (wavelength_nm >= start) & (wavelength_nm <= end)
            trans[mask] = 1.0
            bands[k] = pd.DataFrame({'wavelength_nm': wavelength_nm, 'transmission': trans})

    results = {}

    if stellar_flux is not None:
        if stellar_radius is None or planet_radius is None:
            raise ValueError("Please provide stellar_radius in R_sun and planet_radius in R_jup in order to calculate Fp/Fs.")
        if stellar_wavelength is None:
            raise ValueError("If stellar_flux is provided, stellar_wavelength must also be provided.")
        # Convert stellar wavelengths to nm
        stellar_wavelength_nm = (np.array(stellar_wavelength) * Unit(stellar_wavelength_unit)).to(u.nm).value
        # Interpolate stellar flux onto the same wavelength_nm grid
        if stellar_flux is not None and stellar_wavelength_nm is not None:
            stellar_flux_interp = np.interp(wavelength_nm, stellar_wavelength_nm, stellar_flux)
    else:
        stellar_flux_interp = None
    if stellar_radius is not None and planet_radius is not None:
        # Convert to meters
        R_s = stellar_radius * R_sun
        R_p = planet_radius * R_jup
        area_ratio = (R_p / R_s)**2
    else:
        area_ratio = None


    for band_label, trans_df in bands.items():
        interp_trans = np.interp(wavelength_nm, trans_df['wavelength_nm'], trans_df['transmission'])
        fwhm_mask = fwhm_clip(wavelength_nm, interp_trans)

        flux_vals = flux_array[fwhm_mask]
        trans_vals = interp_trans[fwhm_mask]

        numerator = np.sum(flux_vals * trans_vals)
        denominator = np.sum(trans_vals)
        mean_flux = numerator / denominator if denominator > 0 else 0.0

        stellar_mean = None
        if stellar_flux_interp is not None:
            stellar_flux_vals = stellar_flux_interp[fwhm_mask]
            stellar_numerator = np.sum(stellar_flux_vals * trans_vals)
            stellar_mean = stellar_numerator / denominator if denominator > 0 else 0.0

            # Apply geometric scaling if radii are provided
            if area_ratio is not None and stellar_mean > 0:
                fp_over_fs = (mean_flux / stellar_mean) * area_ratio.value
        else:
            fp_over_fs = None

        fwhm_mask_native = fwhm_clip(trans_df['wavelength_nm'].values, trans_df['transmission'].values)
        fwhm_range = trans_df['wavelength_nm'].values[fwhm_mask_native].min(), trans_df['wavelength_nm'].values[fwhm_mask_native].max()

        results[band_label] = {
            'mean': mean_flux,
            'stellar_mean': stellar_mean,
            'fp_over_fs': fp_over_fs,
            'formatted': f"{mean_flux:.2e}",
            'fwhm_range_nm': fwhm_range,
            'color': {
                'B1': 'green',
                'B2': 'orange',
                'B3': 'red',
                'B4': 'purple'
            }.get(band_label, 'gray')
        }

        print(f"\n--- Band {band_label[1]} ---")
        print(f"  Wavelength Range   : {fwhm_range[0]:.1f} – {fwhm_range[1]:.1f} nm")
        print(f"  Mean Planetary Flux: {mean_flux:.2e} [{flux_unit}]")

        if stellar_flux_interp is not None:
            print(f"  Mean Stellar Flux  : {stellar_mean:.2e} [{stellar_flux_unit}]")
            if area_ratio is not None:
                print(f"  (R_p / R_s)^2      : {area_ratio.value:.2e}")
                print(f"  Flux Ratio (Fp/Fs) : {fp_over_fs:.2e}")
            else:
                print("  [!] Provide stellar_radius (R_sun) and planet_radius (R_jup) to compute Fp/Fs")


    if plot:
        _plot_band_fluxes(wavelength_nm, flux_array, flux_unit, results,
                        bands=bands if show_color_filter_curves else None,
                        stellar_wavelength=stellar_wavelength,
                        stellar_flux_interp=stellar_flux_interp,stellar_flux_unit=stellar_flux_unit,
                        show_color_filter_curves=show_color_filter_curves,
                        log_y=log_y)

# ----------------------
# Plotting function
# ----------------------
def _plot_band_fluxes(wavelength_nm, flux_array, flux_unit, results,
                      bands=None,
                      stellar_wavelength=None, stellar_flux_interp=None, stellar_flux_unit=None,
                      flux_unit_label="Flux", show_color_filter_curves=False,
                      log_y=False):

    # Plot planetary flux
    fig, ax = plt.subplots(figsize=(10, 5))
    x_min, x_max = 400, 1000
    zoom_mask = (wavelength_nm >= x_min) & (wavelength_nm <= x_max)
    zoomed_flux = flux_array[zoom_mask]
    max_flux = np.max(zoomed_flux)
    min_flux = np.min(zoomed_flux[zoomed_flux > 0]) if log_y else 0

    ax.plot(wavelength_nm, flux_array, color='black', label="Planet Spectrum")

    for band_label, result in results.items():
        fwhm_min, fwhm_max = result['fwhm_range_nm']
        color = result['color']
        val = result['mean']
        ax.hlines(val, fwhm_min, fwhm_max, color=color, linewidth=2,
                  label=f"{band_label} mean")

    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel(f"{flux_unit}")
    ax.set_title("Planet Flux")
    if log_y:
        ax.set_yscale("log")
        ax.set_ylim(min_flux, max_flux * 1.1)
    else:
        ax.set_ylim(0, max_flux * 1.1)

    if bands and show_color_filter_curves:
        ax_trans = ax.twinx()
        ax_trans.set_ylabel("In-Band Transmission", color='gray')
        ax_trans.set_ylim(0, 1.05)
        ax_trans.tick_params(axis='y', labelcolor='gray')
        for band_label, trans_df in bands.items():
            ax_trans.plot(trans_df['wavelength_nm'],
                          trans_df['transmission'],
                          linestyle=':', color='gray', alpha=0.7)
            center = (trans_df['wavelength_nm'].min() + trans_df['wavelength_nm'].max()) / 2
            ax_trans.text(center, 1.01, band_label, ha='center', color='gray')

    plt.tight_layout()
    plt.show()

    # Plot stellar flux if provided
    if stellar_flux_interp is not None:
        fig, ax = plt.subplots(figsize=(10, 5))
        zoomed_flux = stellar_flux_interp[zoom_mask]
        max_flux = np.max(zoomed_flux)
        min_flux = np.min(zoomed_flux[zoomed_flux > 0]) if log_y else 0

        ax.plot(wavelength_nm, stellar_flux_interp, color='black', label="Stellar Spectrum")

        for band_label, result in results.items():
            fwhm_min, fwhm_max = result['fwhm_range_nm']
            color = result['color']
            val = result['stellar_mean']
            if val is not None:
                ax.hlines(val, fwhm_min, fwhm_max, color=color, linewidth=2,
                          label=f"{band_label} stellar mean")

        ax.set_xlim(x_min, x_max)
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel(f"Flux {stellar_flux_unit}")
        ax.set_title("Stellar Flux")
        if log_y:
            ax.set_yscale("log")
            ax.set_ylim(min_flux, max_flux * 1.1)
        else:
            ax.set_ylim(0, max_flux * 1.1)

        if bands and show_color_filter_curves:
            ax_trans = ax.twinx()
            ax_trans.set_ylabel("Transmission", color='gray')
            ax_trans.set_ylim(0, 1.05)
            ax_trans.tick_params(axis='y', labelcolor='gray')
            for band_label, trans_df in bands.items():
                ax_trans.plot(trans_df['wavelength_nm'],
                              trans_df['transmission'],
                              linestyle=':', color='gray', alpha=0.7)
                center = (trans_df['wavelength_nm'].min() + trans_df['wavelength_nm'].max()) / 2
                ax_trans.text(center, 1.01, band_label, ha='center', color='gray')

        plt.tight_layout()
        plt.show()
