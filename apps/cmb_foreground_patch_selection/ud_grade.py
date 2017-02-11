#! /usr/bin/env python

import os
import pylab as pl
import healpy
from leap.lib.leap_app import leap_app
from leap.lib.geometry.healpy_map_conversions import map_gal_to_eq, map_eq_to_gal


class ConversionApp(leap_app.App):

    def run(self):
        if self.settings.polarization:
            in_map = pl.array(healpy.read_map(self.settings.in_filename, field=(0, 1, 2)))
            range_ = [0, 1, 2]
        else:
            in_map = pl.array([healpy.read_map(self.settings.in_filename)])
            range_ = [0]
        if self.settings.plotting:
            healpy.mollview(in_map[0], sub=211)
        if self.settings.smooth:
            in_map_copy = in_map.copy()
            in_map_copy[pl.isnan(in_map)] = 0.0
            in_map_copy[in_map == healpy.UNSEEN] = 0.0
            for i in range_:
                in_map_copy[i] = healpy.smoothing(in_map_copy[i], fwhm=self.settings.fwhm)
            in_map = in_map_copy.copy()
            in_map[in_map_copy == 0] = pl.nan
        in_nside = healpy.npix2nside(in_map[0].size)
        if in_nside != self.settings.out_nside:
            out_map = []
            for i in range_:
                if self.settings.map_type == "hitmap":
                    out_map.append(healpy.ud_grade(in_map[i], self.settings.out_nside, 
                                   order_in="RING", order_out="RING", power=-2))
                elif self.settings.map_type == "signal":
                    out_map.append(healpy.ud_grade(in_map[i], self.settings.out_nside, 
                                   order_in="RING", order_out="RING"))
                elif self.settings.map_type == "variance":
                    out_map.append(healpy.ud_grade(in_map[i], self.settings.out_nside, 
                                   order_in="RING", order_out="RING", power=2))
            out_map = pl.array(out_map)
        else:
            out_map = in_map.copy()
        if self.settings.convert_gal_to_eq:
            for i in range_:
                out_map[i] = map_gal_to_eq(out_map[i])
        if self.settings.convert_eq_to_gal:
            for i in range_:
                out_map[i] = map_eq_to_gal(out_map[i])
        self.create_output()
        out_full_filename = os.path.join(self.out_path, self.settings.out_filename)
        if not self.settings.polarization:
            healpy.write_map(out_full_filename, out_map[0], coord=self.settings.coord)
        else:
            healpy.write_map(out_full_filename, out_map, coord=self.settings.coord)
        if self.settings.plotting:
            out_map = healpy.read_map(out_full_filename)
            healpy.mollview(out_map, sub=212, title='smoothed and converted map')
            pl.show()


if __name__ == "__main__":
    app = ConversionApp()
    app.run()
    app.end()
