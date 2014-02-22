# Create H2CO modeling data by running RADEX
execfile('create_radex_grid.py')

# Create the optical depth ratio and density cubes
execfile('tau_ratio_cube.py')

# Do some continuum comparisons
execfile('comparison_maps.py')
execfile('continuum_cband_law_comparison.py')

# Integral images of Brick
# (must do this before bgsub_image_figures)
execfile('integrate_brick.py')
execfile('integrate_fitgaussian.py')

# Try background subtraction
execfile('create_background_cube.py')
execfile('bgsub_plotting.py') # must come before bgsub_image_figures
execfile('bgsub_image_figures.py')

# Examine the ratios [should probably be in plots]
# requires bgsub_plotting
execfile('ratio_dens_histogram.py')

# Extract some spectra and plot them
execfile('extract_spectra.py')
execfile('spectral_grid.py')
