from leap.lib.tools import generic
from leap.lib.units.angles import *

settings = generic.Class()

settings.nside = 2048
settings.frequency = 1200 # 600, 860, 1200
settings.show = False
settings.save = True



########## ud grade settings ##############
settings.polarization = True
settings.in_filename = "./maps/py_sm/2048/2048_MJySr_cmb_freef_spinn_synch_therm_%dp0_2048.fits" %settings.frequency
settings.out_filename = "512_MJySr_cmb_freef_spinn_synch_therm_%dp0_512.fits" %settings.frequency
settings.out_nside = 512
settings.smooth = False
settings.convert_gal_to_eq = False
settings.convert_eq_to_gal = False
settings.map_type = "signal" #or "hitmap" if hitmap. This is for the nside conversion.


settings.fwhm = from_arcmin(8.0)
settings.coord = "GALACTIC"
settings.plotting = False
