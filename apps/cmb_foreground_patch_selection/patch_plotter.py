#! /usr/bin/python

import pylab as pl
import numpy as np
import healpy as hp
import matplotlib as mpl
import os

from leap.lib.leap_app import leap_app
from leap.lib.units.angles import *
from leap.lib.geometry import coordinates


fontsize = 22
mpl.rcParams["font.size"] = fontsize
mpl.rcParams["xtick.labelsize"] = fontsize
mpl.rcParams["ytick.labelsize"] = fontsize
mpl.rcParams["legend.fontsize"] = fontsize
mpl.rcParams["axes.titlesize"] = fontsize


class Patch(object):

    def __init__(self, name, center, xsize_deg, ysize_deg, obs_time_hours, delta_label_deg):
        self.name = name
        self.center = center
        self.xsize = from_degrees(xsize_deg)
        self.ysize = from_degrees(ysize_deg)
        self.xsize_deg = xsize_deg
        self.ysize_deg = ysize_deg
        self.obs_time = obs_time_hours / 3600.0
        self.obs_time_hours = obs_time_hours
        self.delta_label = delta_label_deg
        self.ra, self.dec = coordinates.gal_to_eq(from_degrees(center[0]), from_degrees(center[1]))

    def plot_contour(self):
        lon_min = self.center[0] - self.xsize_deg/2.0
        lon_max = self.center[0] + self.xsize_deg/2.0
        lat_min = self.center[1] - self.ysize_deg/2.0
        lat_max = self.center[1] + self.ysize_deg/2.0
        lons = np.array([lon_min, lon_max, lon_max, lon_min, lon_min])
        lats = np.array([lat_min, lat_min, lat_max, lat_max, lat_min])
        hp.projplot(self.upsample(lons), self.upsample(lats), '-', color='green', lw=2, lonlat=True, coord="G")

    def upsample(self, x, delta=100):
        return np.interp(np.linspace(0, x.size-1, delta), np.arange(x.size), x)

    def __str__(self):
        return "%s: RA = %s, Dec = %s" %(self.name, to_hours_minutes_seconds_string(self.ra), to_degrees_minutes_seconds_string(self.dec))
        
class PatchPlotter(leap_app.App):

    def plot_selected_patch(self, patch):
        print patch
        delta = 2.0
        phis = np.round(np.arange(patch.center[0]-patch.xsize_deg/2.0, patch.center[0]+patch.xsize_deg/2.0, patch.delta_label))
        thetas = np.round(np.arange(patch.center[1]-patch.ysize_deg/2.0, patch.center[1]+patch.ysize_deg/2.0, patch.delta_label))
        reso = 1.7
        #reso = 1.0
        xsize = to_arcmin(from_degrees(patch.xsize_deg))/reso
        ysize = to_arcmin(from_degrees(patch.ysize_deg))/reso
        # Remove mean of Pol over patch for better SNR estimate
        imap = hp.gnomview(self.maps[0], rot=patch.center, coord="G", reso=reso, xsize=xsize, ysize=ysize, return_projected_map=True, fig=5)
        qmap = hp.gnomview(self.maps[1], rot=patch.center, coord="G", reso=reso, xsize=xsize, ysize=ysize, return_projected_map=True, fig=5)
        umap = hp.gnomview(self.maps[2], rot=patch.center, coord="G", reso=reso, xsize=xsize, ysize=ysize, return_projected_map=True, fig=5)
        imean = imap.flatten().mean()
        qmean = qmap.flatten().mean()
        umean = umap.flatten().mean()
        pol_map = np.sqrt((self.maps[1]-qmean)**2 + (self.maps[2]-umean)**2)
        print 'In patch', patch.name, 'the means in MJy/sr are Qmean = %f and Umean = %f' %(qmean, umean)
        pl.close(5)
        pl.suptitle("%s Gal Coord Center: %s" %(patch.name, str(patch.center)), y=0.999)
        # I
        hp.gnomview(self.maps[0]-imean, 
                             rot=patch.center, coord="G", min=-6.0, max=6.0, cbar=False, reso=reso, xsize=xsize, ysize=ysize, sub=(2, 2, 1), notext=True, cmap="seismic", title="")
        self.make_cbar("Mean removed I (MJy/sr)", "%1.1f")
        self.make_coord_labels_patches(phis, thetas, patch, coord="G", lw=1)
        # P
        hp.gnomview(pol_map, rot=patch.center, coord="G", min=0, max=1.0, cbar=False, reso=reso, xsize=xsize, ysize=ysize, sub=(2, 2, 2), notext=True, cmap="seismic", title="")
        self.make_cbar("Differential P (MJy/sr)", "%1.1f")
        self.make_coord_labels_patches(phis, thetas, patch, coord="G", lw=1)
        # p
        hp.gnomview(self.maps[3]/abs(self.maps[0])*100, 
                             rot=patch.center, coord="G", min=0, max=30, cbar=False, reso=reso, xsize=xsize, ysize=ysize, sub=(2, 2, 3), notext=True, cmap="seismic", title="")
        self.make_cbar("P/I (percent) ", "%1.1f")
        self.make_coord_labels_patches(phis, thetas, patch, coord="G", lw=1)
        # Pol SNR
        sn_map, pix_area, patch_area = self.get_snr_map(patch.xsize_deg, patch.ysize_deg, patch.obs_time_hours, pol_map, pix_area=from_arcmin(reso)*from_arcmin(reso))
        #sn_map, pix_area, patch_area = self.get_snr_map(patch.xsize_deg, patch.ysize_deg, patch.obs_time_hours, pol_map)
        sn_gnom_map = hp.gnomview(sn_map, rot=patch.center, coord="G", min=0, max=10.0, cbar=False, reso=reso, xsize=xsize, ysize=ysize, unit="unitless", sub=(2, 2, 4), notext=True, title="", cmap="seismic", return_projected_map=True)
        self.make_coord_labels_patches(phis, thetas, patch, coord="G", lw=1)
        self.make_cbar("Polarized SNR", "%1.0f") #\n(mean removed in Q/U)
        self.save_and_show(os.path.join(self.out_path, "%s_SNR_%dGHz_%d_%2.0fdeg_%2.0fh.png" %(patch.name, self.settings.frequency, self.settings.nside, to_degrees(to_degrees(patch_area)), patch.obs_time_hours)), 
                           figsize=(16, 10), save=self.settings.save, show=self.settings.show) 
        print "For patch %s, there are %d pixels with SNR > 3" %(patch.name, (sn_gnom_map > 3).sum())

    def make_cbar(self, label, format_, ticks=None):
        fig = pl.gcf()
        sub = pl.gca()
        im = sub.get_images()[0]
        cbaxes, kwd = mpl.colorbar.make_axes(sub, orientation='vertical', shrink=0.8) #, anchor=(0.5, 1.0)) #, , location="bottom")
        cbar = fig.colorbar(cax=cbaxes, mappable=im, orientation='vertical', format=format_, ticks=ticks)
        cbar.set_label(label, labelpad=2, fontsize=(fontsize-4))
        cbar.ax.tick_params(labelsize=fontsize-4, width=1.0, length=4.0)

    def get_nan_pixels_for_orientation_check(self):
        import astropy.coordinates
        from leap.resources.constants import astro_locations
        disc_cartesian = astropy.coordinates.spherical_to_cartesian(1, astro_locations.rcw38_gal[1], astro_locations.rcw38_gal[0])
        rcw38_disc_pixels = hp.query_disc(self.settings.nside, disc_cartesian, from_degrees(1.0))
        invalid_pixels = rcw38_disc_pixels
        self.maps[0][invalid_pixels] = np.nan
        gal_pixels = hp.query_strip(self.settings.nside, np.pi/2.0-from_degrees(0.2), np.pi/2.0+from_degrees(0.2))
        self.maps[0][gal_pixels] = np.nan

    def load_pysm(self):
        nside = self.settings.nside
        map_path = os.path.join("./maps/py_sm/%d" %nside, "%d_MJySr_cmb_freef_spinn_synch_therm_%dp0_%d.fits"%(nside, self.settings.frequency, nside))
        self.maps = list(hp.read_map(map_path, field=(0, 1, 2)))
        self.maps.append(np.sqrt(self.maps[1]**2 + self.maps[2]**2))
        self.maps = np.array(self.maps)

    def get_snr_map(self, patch_xsize_deg, patch_ysize_deg, obs_time_hours, map_, pix_area=None):
        if pix_area is None:
            pix_area = hp.nside2pixarea(hp.npix2nside(self.maps[0].size))
        patch_area = from_degrees(patch_xsize_deg)*from_degrees(patch_ysize_deg)
        sensitivity = self.get_sensitivity(patch_area, pix_area, obs_time_hours)
        sn_map = map_/sensitivity
        return sn_map, pix_area, patch_area

    def plot_mollzoom_to_search_for_best_patch(self, patch_xsize_deg, patch_ysize_deg, obs_time_hours):
        # Intensity
        #hp.mollzoom(np.log10(self.maps[0]), title="%d GHz Intensity" %self.settings.frequency, rot=(0, 180, 180), 
        #                     coord=["G", "C"], cmap='seismic', min=-1, max=3.0, unit="log10(MJy/sr)")
        self.get_nan_pixels_for_orientation_check()
        hp.mollzoom(np.log10(self.maps[0]), title="%d GHz Intensity" %self.settings.frequency, 
                             cmap='seismic', min=-1, max=3.0, unit="log10(MJy/sr)", coord=["G"])
        hp.graticule()
        hp.set_g_clim(-1, 3.0)
        self.plot_blast_contour(coord="C")
        # SNR
        pol_map = np.sqrt((self.maps[1]-self.maps[1].mean())**2 + (self.maps[2]-self.maps[2].mean())**2)
        sn_map, pix_area, patch_area = self.get_snr_map(patch_xsize_deg, patch_ysize_deg, obs_time_hours, pol_map)
        #hp.mollzoom(sn_map, title="SNR at %d GHz, patch area %3.0f sq deg\npix size = %1.1f\', obs time %2.0f hours" %(self.settings.frequency, to_degrees(to_degrees(patch_area)), 
        #            to_arcmin(np.sqrt(pix_area)), obs_time_hours), coord=["G", "C"], cmap='seismic', min=0, max=10.0, unit="unitless", rot=(0, 180, 180))
        hp.mollzoom(sn_map, title="SNR at %d GHz, patch area %3.0f sq deg\npix size = %1.1f\', obs time %2.0f hours" %(self.settings.frequency, to_degrees(to_degrees(patch_area)), 
                    to_arcmin(np.sqrt(pix_area)), obs_time_hours), coord=["G"], cmap='seismic', min=0, max=10.0, unit="unitless")
        hp.graticule()
        hp.set_g_clim(0, 10.0)
        self.plot_blast_contour(coord="C")
        # Pol fraction
        hp.mollzoom((self.maps[3]/abs(self.maps[0]))*100, title="%d GHz Pol Fraction" %self.settings.frequency, 
                             coord=["G"], cmap='seismic', min=0.0, max=20.0, unit="Percent")
        hp.graticule()
        hp.set_g_clim(0, 20.0)
        self.plot_blast_contour(coord="C")
        pl.show()

    def plot_I_P_p_fullsky(self, patches, plot_patch=False):
        pl.suptitle("BLAST Coverage (Equatorial Coordinates)")
        lon_range = [-115, 65] # [-135, 120]
        lat_range = [-70, -10] # [-90, 0]
        phis = pl.arange(-120, 60, 30)
        thetas = pl.arange(-60, 0, 30)
        hp.cartview(np.log10(self.maps[0]), title="%d GHz Intensity" %self.settings.frequency, 
                             rot=(0, 180, 180), coord=["G", "C"], cmap='seismic', min=-1, max=3.0, cbar=True, lonra=lon_range, latra=lat_range, unit="log10(MJy/sr)", sub=(2, 2, 1))
        self.plot_blast_contour()
        if plot_patch:
            for patch in patches:
                patch.plot_contour()
        self.make_coord_labels(phis, thetas, 30, coord="C")
        hp.cartview(np.log10(self.maps[3]), title="%d GHz Pol Power" %self.settings.frequency, 
                             rot=(0, 180, 180,), coord=["G", "C"], cmap='seismic', min=-1, max=0.2, cbar=True, lonra=lon_range, latra=lat_range, unit="log10(MJy/sr)", sub=(2, 2, 2))
        self.plot_blast_contour()
        if plot_patch:
            for patch in patches:
                patch.plot_contour()
        self.make_coord_labels(phis, thetas, 30)
        hp.cartview((self.maps[3]/abs(self.maps[0]))*100, title="%d GHz Pol Fraction" %self.settings.frequency, 
                             rot=(0, 180, 180,), coord=["G", "C"], cmap='seismic', min=0.0, max=20.0, cbar=True, lonra=lon_range, latra=lat_range, unit="Percent", sub=(2, 2, 3))
        self.plot_blast_contour()
        if plot_patch:
            for patch in patches:
                patch.plot_contour()
        self.make_coord_labels(phis, thetas, 30)
        self.save_and_show(os.path.join(self.out_path, "I_P_p_%dGHz_%d.png" %(self.settings.frequency, self.settings.nside)), 
                           figsize=(16, 10), save=self.settings.save, show=self.settings.show)
        pl.clf()
        pl.close()

    def plot_sn_fullsky(self, patch):
        phis = pl.arange(-120, 60, 30)
        thetas = pl.arange(-60, 0, 30)
        sn_map, pix_area, patch_area = self.get_snr_map(patch.xsize_deg, patch.ysize_deg, patch.obs_time_hours, self.maps[3]-self.maps[3].mean())
        lon_range = [-115, 65] # [-135, 120]
        lat_range = [-70, -10] # [-90, 0]
        hp.cartview(sn_map, title="%s SNR at %d GHz, patch area %3.0f sq deg\npix size = %1.1f\', obs time %2.0f hours" %(patch.name, self.settings.frequency, to_degrees(to_degrees(patch_area)), 
                    to_arcmin(np.sqrt(pix_area)), patch.obs_time_hours), rot=(0, 180, 180), coord=["G", "C"], cmap='seismic', min=0, max=10.0, 
                    cbar=True, lonra=lon_range, latra=lat_range, unit="unitless")
        self.plot_blast_contour()
        patch.plot_contour()
        self.make_coord_labels(phis, thetas, 30)
        self.save_and_show(os.path.join(self.out_path, "%s_full_SNR_%dGHz_%d_%2.0fdeg_%2.0fh.png" %(patch.name, self.settings.frequency, self.settings.nside, to_degrees(to_degrees(patch_area)), patch.obs_time_hours)), 
                           figsize=(16, 10), save=self.settings.save, show=self.settings.show)
        pl.clf()
        pl.close()

    def get_sensitivity(self, patch_area, pixel_area, obs_time_hours):
        # Nf = {1200: 0.3856, 860: 0.2873, 600: 0.1547} # in MJy/sr sqrt(s)
        Nf = {1200: 0.318, 860: 0.237, 600: 0.127} # in MJy/sr sqrt(s), less pessimistic
        scanning_time = obs_time_hours*3600.0 # 4 days in seconds
        sensitivity = Nf[self.settings.frequency] * np.sqrt(patch_area / pixel_area) * 1.0/np.sqrt(scanning_time)
        return sensitivity

    def plot_blast_contour(self, coord="C"):
        lon_10h = np.array([910,776,1050,2566,2690,2421,910])*(360.0/4000)
        lat_10h = np.array([800,644,339,345,655,794,800])*(180.0/2000)-90
        hp.projplot(self.upsample(lon_10h), self.upsample(lat_10h), '-', color='yellow', lw=2, lonlat=True, coord=coord)

    def upsample(self, x, delta=100):
        return np.interp(np.linspace(0, x.size-1, delta), np.arange(x.size), x)
        
    def make_coord_labels(self, phis, thetas, delta_grid, coord="C", lw=2):
        hp.graticule(coord=coord, dpar=delta_grid, dmer=delta_grid, color="white", local=False, linestyle=":", lw=lw, verbose=False)
        xlabels = ["%d$\degree$" %label for label in (180+phis)]
        ylabels = ["%d$\degree$" %label for label in (thetas)]
        for i in range(phis.size):
            hp.projtext(phis[i], -70, xlabels[i], lonlat=True, fontsize=20, ha='center', va="top", direct=True)
        for i in range(thetas.size):
            hp.projtext(60, thetas[i], ylabels[i], lonlat=True, fontsize=20, va='center', ha="right", direct=True)

    def make_coord_labels_patches(self, phis, thetas, patch, coord="C", lw=2):
        hp.graticule(coord=coord, dpar=patch.delta_label/2.0, dmer=patch.delta_label/2.0, color="white", local=False, linestyle=":", lw=lw, verbose=False)
        xlabels = ["%d$\degree$" %label for label in (phis)]
        ylabels = ["%d$\degree$" %label for label in (thetas)]
        for i in range(phis.size):
            hp.projtext(phis[i], patch.center[1]-patch.ysize_deg/2.0, xlabels[i], lonlat=True, fontsize=fontsize-4, ha='center', va="top", direct=False)
        for i in range(thetas.size):
            hp.projtext(patch.center[0]+patch.xsize_deg/2.0+0.2, thetas[i], ylabels[i], lonlat=True, fontsize=fontsize-4, va='center', ha="right", direct=False)

    def run(self):
        self.load_pysm()
        patches = [Patch("patch0", [-106.41, 10.631, 0], 5.0, 5.0, 96.0, 2.0)] #, 
                         #Patch("patch0b", [-110.5, 14.25, 0], 5.0, 5.0, 96.0, 2.0)] #,
                         #Patch("patch1a", [-119.424, -22.331, 0], 2.0, 2.0, 48.0, 1.0), Patch("patch1b", [-46.996, 31.454, 0], 2.0, 2.0, 48.0, 1.0)] #,
                         #Patch("patch2", [-42.234, 11.614, 0], 1.0, 10.0, 96.0, 1.0)]
        #self.plot_I_P_p_fullsky(patches, plot_patch=True)
        self.patches = patches
        if False:
            for patch in patches:
                self.plot_sn_fullsky(patch)
                #self.plot_mollzoom_to_search_for_best_patch(patch.xsize_deg, patch.ysize_deg, patch.obs_time_hours)
        if True:
            for patch in patches:
                self.plot_selected_patch(patch)

if __name__ == "__main__":
    patch_plotter = PatchPlotter()
    patch_plotter.run()
    patch_plotter.end() 
