from astropy.io import fits as pf
import glob
import numpy as np
from matplotlib import pyplot as pl

flist = glob.glob('*.fits')
print("Found %d fits files" % len(flist))

pl.ioff()

for f in flist:
    ff = pf.open(f)
    wave = ff[1].data['WAVELENGTH']
    flux = ff[1].data['FLUX']
    ff.close()
    name = f.split('.fits')[0]
    wlmin = np.nanmin(wave)
    wlmax = np.nanmax(wave)
    wllen = len(wave)
    print("%20s %10.2f %10.2f %d" % (f, wlmin, wlmax, wllen))
    fig = pl.figure()
    pl.yscale('log')
    pl.plot(wave, flux, "+")
    pl.xlim([wlmin-100., np.min([12000, wlmax])+100.])
    if wlmax > 12100:
        seen = flux[np.where(wave <= 12100.)]
        flmin = np.nanmin(seen)
        flmax = np.nanmax(seen)
        pl.ylim((flmin, flmax))
    pl.title("%s: %.2f - %.2f A, %d pts" % (name, wlmin, wlmax, wllen))
    pl.xlabel('WAVELENGTH')
    pl.ylabel('FLUX')
    pl.savefig(name+'.png')
    pl.close(fig)

