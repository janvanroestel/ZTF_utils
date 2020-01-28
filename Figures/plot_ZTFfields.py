import numpy as np

from matplotlib import rc
rc('text', usetex=False)
import matplotlib.pyplot as plt



import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import healpy as hp
from astropy.io import fits
from astropy.io.votable import parse
import sys
from healpy.visufunc import mollview,gnomview,cartview,orthview
import argparse

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

# do not use latex


# plot arguments
parser = argparse.ArgumentParser(description='Plot ZTF and TESS coverage')
parser.add_argument('--fieldlists', nargs='+', default=[],
                    help='list of filenames with fieldIDs')
parser.add_argument('--TESSsector','-t', type=int, default=None,
                    help='sectors to show')
parser.add_argument('-c', type=str, default='C',
                    help='coordinate system to plot')
parser.add_argument('-i', type=bool, default=False,
                    help='show interactive plot')
parser.add_argument('-p', nargs='+', default=[],
                    help='list of dates for which to plot opposition lines')
parser.add_argument('-f', type=bool, default=False,
                    help='print FIDs on figure')
parser.add_argument('-m', type=str, default=None,
                    help='show moon')
parser.add_argument('-b', type=str, default='Gaia',
                    help='background hist')
args = parser.parse_args()

print(args)

datadir = "../Data/"

# load field grid
colnames = "ID Ra Dec Ebv GalLong GalLat EclLong EclLat Entry".split()
ZTFgrid = np.genfromtxt(datadir+'ZTF_Fields.txt',dtype=float,comments='#',names=colnames)
ZTFfields = dict(zip(ZTFgrid['ID'], ZTFgrid))

# init figure
fig1 = plt.figure(figsize=(14,8))

if args.b=='Gaia':
    # background setup
    NSIDE = 256
    NPIX = hp.nside2npix(NSIDE)
    coordsys = ['C',args.c]
    nest = True

    # colormap
    cm = plt.cm.get_cmap('viridis') # colorscale
    cm.set_under('w')
    cm.set_bad('w')

    # load the data
    hdulist = fits.open(datadir+'Gaia_hp8_densitymap.fits')
    hist = hdulist[1].data['srcdens'][np.argsort(hdulist[1].data['hpx8'])]

    # plot the data in healpy
    norm ='log' 
    hp.mollview(hist,fig=1,norm=norm,unit='Stars per sq. arcmin.',cbar=True,nest=nest,title='',
        coord=coordsys,notext=False,cmap=cm,flip='astro',nlocs=4)



# borders
lw = 3
pi = np.pi
dtor = pi/180.
theta = np.arange(0,181)*dtor
hp.projplot(theta, theta*0-pi,'-k',
                               lw=lw,direct=True)
hp.projplot(theta, theta*0+0.9999*pi,'-k',
                               lw=lw,direct=True)
phi = np.arange(-180,180)*dtor
hp.projplot(phi*0+1.e-10, phi,'-k',
                               lw=lw,direct=True)
hp.projplot(phi*0+pi-1.e-10, phi,'-k',
                               lw=lw,direct=True)

# galaxy
for gallat in [15,0,-15]:
    theta = np.arange(0., 360, 0.036)
    phi = gallat*np.ones_like(theta)
    hp.projplot(theta, phi, 'w-', coord=['G'],lonlat=True,lw=2)

# ecliptic
for ecllat in zip([0,-30,30],[2,1,1]):
    theta = np.arange(0., 360, 0.036)
    phi = gallat*np.ones_like(theta)
    hp.projplot(theta, phi, 'w-', coord=['E'],lonlat=True,lw=2,ls=':')

# graticule
hp.graticule(ls='-',alpha=0.1,lw=0.5)

# NWES
fig = plt.gcf()
ax = plt.gca()
plt.text(0.0,0.5,r'E',ha='right',transform=ax.transAxes,weight='bold')
plt.text(1.0,0.5,r'W',ha='left',transform=ax.transAxes,weight='bold')
plt.text(0.5,0.992,r'N',va='bottom',ha='center',transform=ax.transAxes,weight='bold')
plt.text(0.5,0.0,r'S',va='top',ha='center',transform=ax.transAxes,weight='bold')


for p in args.p:
    # plot opposition 
    from astropy.time import Time
    from astropy.coordinates import solar_system_ephemeris, EarthLocation
    from astropy.coordinates import get_sun# get_body_barycentric, get_body, get_moon

    # calculate RA of opposition
    t = Time("%s 00:00" %p)
    with solar_system_ephemeris.set('builtin'):
        opp_ra = (get_sun(t).ra.deg+180)%360

    # plot line
    for RA in [opp_ra,]:
        phi = np.arange(-90., 90, 0.036)
        hp.projplot(opp_ra, phi, 'w-', coord=['C'],lonlat=True,lw=2,ls=':')

if args.m is not None:
    from astropy.time import Time
    from astropy.coordinates import solar_system_ephemeris, EarthLocation
    from astropy.coordinates import get_moon # get_body_barycentric, get_body, get_moon

    # calculate RA of opposition
    t = Time("%s 00:00" %args.m)
    with solar_system_ephemeris.set('builtin'):
        moonloc = get_moon(t)

    hp.projplot(moonloc.ra, moonloc.dec, 'wo',ms=10, coord=['C'],lonlat=True)



# list of fields to plot
fieldlists = args.fieldlists
if len(fieldlists)>0:
    fields = [np.genfromtxt(f) for f in fieldlists]
    colours1 = ['C3','C1','C6','C7','C8',]
    for fields,c in zip(fields,colours1):
        for line in fields:
            try:
                alpha = 1.0
                temp = (datadir+'/field_edges/%d.dat' %(line))
                edges = np.loadtxt(temp, unpack=True)
                hp.projplot(edges,c,coord=coordsys,lonlat=True,lw=2,alpha=alpha)
                if args.f:
                    hp.projtext(ZTFfields[line][1],ZTFfields[line][2],
                    s='%d'%line,text='%d'%line,
                    coord=coordsys,
                    lonlat=True,
                    fontsize=8,
                    color=c,zorder=5,
                    ha='center',va='center')
            except Exception as e:
                print("Could not plot field %s \n" %line, e)



# plot background grid
colours2 = ['k',]
allfields = np.arange(245,882,1)
for fields,c in zip([allfields,],colours2):
    for line in fields:
        try:
            temp = (datadir+'/field_edges/%d.dat' %(line))
            edges = np.loadtxt(temp, unpack=True)
            hp.projplot(edges,c,coord=coordsys,lonlat=True,lw=1,alpha=0.1)
            if args.f:
                color = 'w'
                if args.m is not None:
                    from astropy.time import Time
                    from astropy.coordinates import solar_system_ephemeris, EarthLocation
                    from astropy.coordinates import get_moon # get_body_barycentric, get_body, get_moon

                    # calculate RA of opposition
                    t = Time("%s 00:00" %args.m)
                    with solar_system_ephemeris.set('builtin'):
                        moonloc = get_moon(t)

                    c1 = SkyCoord(ZTFfields[line][1]*u.degree,
                                  ZTFfields[line][2]*u.degree, frame='icrs')
                    sep = moonloc.separation(c1)
                    if sep.deg<35:
                        color = 'k'
                hp.projtext(ZTFfields[line][1],ZTFfields[line][2],
                s='%d'%line,text='%d'%line,
                coord=coordsys,
                lonlat=True,
                fontsize=8,
                color=color,zorder=5,
                ha='center',va='center')
        except Exception as e:
            print(e)

plt.show()

