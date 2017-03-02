#! /usr/bin/python

import numpy as np

from leap.lib.geometry.coordinates import *
from leap.lib.units.angles import *
#from leap.lib.units import times
from leap.lib.leap_app import leap_app

class Patch(object):

    def __init__(self, name, center, xsize_deg, ysize_deg, obs_time_hours, delta_label_deg):
        self.name = name
        self.center = center
        self.lon = from_degrees(center[0])
        self.lat = from_degrees(center[1])
        self.xsize = from_degrees(xsize_deg)
        self.ysize = from_degrees(ysize_deg)
        self.xsize_deg = xsize_deg
        self.ysize_deg = ysize_deg
        self.obs_time = obs_time_hours / 3600.0
        self.obs_time_hours = obs_time_hours
        self.delta_label = delta_label_deg

class ScanStrategy(leap_app.App):

    def run(self):
        patches = [Patch("patch0", [-106.41, 10.631, 0], 5.0, 5.0, 96.0, 2.0)] #, Patch("patch0-peter",  [249.5, 14.25, 0], 5.0, 5.0, 146.0, 2.0),
                   #Patch("patch1a", [-119.424, -22.331, 0], 2.0, 2.0, 48.0, 1.0), Patch("patch1b", [-46.996, 31.454, 0], 2.0, 2.0, 48.0, 1.0), 
                   #Patch("patch2", [-42.234, 11.614, 0], 1.0, 10.0, 96.0, 1.0)]
        for patch in patches:
            self.plot_scanning_patch(patch)

    def plot_scanning_patch(self, patch):
        lst = 0
        lat = from_degrees(-73)
        scan_speed = from_degrees(0.5)
        el_step = from_arcmin(1.0)
        num_steps = int(round(patch.ysize / el_step))
        dt = 0.1
        i = 0
        azs = np.empty((0))
        els = np.empty((0))
        lsts = np.empty((0))
        ra_center, dec_center = gal_to_eq(patch.lon, patch.lat)
        num_points = patch.xsize/scan_speed/dt
        while i < num_steps:
            az_center, el_center = eq_to_hor(ra_center, dec_center, lat, lst)
            if i == 0:
                el = el_center-patch.ysize/2.0
            az_from_throw = (np.linspace(az_center-patch.xsize/2.0, az_center + patch.xsize/2.0, num_points) %(2*np.pi))
            el_from_throw = np.linspace(el, el, num_points)
            lst_from_throw = np.linspace(lst, lst+patch.xsize/scan_speed, num_points)
            el += el_step
            lst += patch.xsize/scan_speed
            azs = np.concatenate((azs, az_from_throw))
            els = np.concatenate((els, el_from_throw))
            lsts = np.concatenate((lsts, lst_from_throw))
            i += 1
        lats = np.linspace(lat, lat, els.size)

        final_ras, final_decs = hor_to_eq(azs, els, lats, lsts)
        pl.title("%s with (l, b) = %s deg, or (ra, dec) = %s deg" %(patch.name, str(patch.center), str([ra_center, dec_center])))
        pl.plot(to_degrees(azs), to_degrees(els), '.')
        #pl.plot(to_degrees(final_ras), to_degrees(final_decs), '.')
        #gal_lats, gal_lons = eq_to_gal(*hor_to_eq(azs, els, lats, lsts))
        #pl.plot(to_degrees(gal_lons), to_degrees(gal_lats), 'o')
        pl.show()


if __name__ == "__main__":
    scan_strategy = ScanStrategy()
    scan_strategy.run()
    scan_strategy.end()
