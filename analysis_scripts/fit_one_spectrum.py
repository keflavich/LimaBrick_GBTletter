"""
This is dumb because it doesn't account for turbulence
"""
import pyspeckit
import pyradex
from pyradex import synthspec
from pyspeckit.spectrum.models.model import SpectralModel

radex_pars = dict(temperature=20, column=1e13,
                  abundance=10**-8.5,
                  collider_densities={'H2':1e4})

R = pyradex.Radex(species='oh2co-h2', **radex_pars)
h2co_spec = synthspec.SyntheticSpectrum.from_RADEX(wcs,
