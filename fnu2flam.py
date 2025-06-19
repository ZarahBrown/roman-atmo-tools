import numpy as np
import astropy.units as u
from astropy.constants import c

def convert_fnu_to_flambda(wavelength, flux_fnu, 
                           flux_unit='Jy', 
                           wavelength_unit='um', 
                           output_unit='W / (m2 um)'):
    """
    Convert F_nu flux values to F_lambda.

    Parameters
    ----------
    wavelength : array-like
        Wavelength values (same length as flux_fnu).
    flux_fnu : array-like
        Flux density values in F_nu units.
    flux_unit : str
        Input flux unit (e.g., 'Jy', 'erg / (cm2 s Hz)').
    wavelength_unit : str
        Unit of wavelength (e.g., 'um', 'nm', 'Angstrom').
    output_unit : str
        Desired output unit for F_lambda.

    Returns
    -------
    wavelength_out : Quantity array
        Wavelength array in the desired unit (same as wavelength_unit).
    flux_flambda : Quantity array
        Converted flux values in F_lambda units.
    """
    try:
        # Parse input units
        wave = wavelength * u.Unit(wavelength_unit)
        flux_fnu_quantity = flux_fnu * u.Unit(flux_unit)

        # Convert wavelength to frequency
        freq = wave.to(u.Hz, equivalencies=u.spectral())

        # Convert F_nu to F_lambda (same shape)
        flambda = flux_fnu_quantity * (freq**2) / c

        # Convert to desired output unit
        flambda_out = flambda.to(u.Unit(output_unit), equivalencies=u.spectral_density(wave))

        return wave, flambda_out

    except Exception as e:
        raise ValueError(f"Conversion failed: {e}")
