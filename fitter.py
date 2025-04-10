import os
import json
import glob
import time
import utils  # Import of utils.py
import getVpec # Import of getVpec.py
import astropy
import requests
import warnings
import datetime
import numpy as np
import time as systime
import snpy as snoopyfit
import sncosmo as salt3fit
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.table import Table, conf
from astroquery.sdss import SDSS
from collections import OrderedDict
from astropy.time import Time as astrotime
from astropy.coordinates import SkyCoord, Galactic
from astro_ghost.ghostHelperFunctions import getTransientHosts, getGHOST

class sneObj:
    # Construction Functions ----------------------------------------------------------------------------------------- #
    def __init__(self, source: str, algo: str, path: str):
        source = source.lower()
        if source == 'empty' or len(path) == 0:
            print('[+++] Creating empty SNe object...')
            pass
        elif source == 'class' and os.path.exists(path):
            print('[+++] Creating SNe object from class file...')
            self.load_class(path)
        elif (source in ['csp-91bg', 'atlas-91bg', 'ztf-91bg', 'combined_91bg',
                         'csp-norm', 'atlas-norm', 'ztf-norm', 'combined_norm']
              and os.path.exists(path)):
            print(f"[+++] Creating '{source}' SNe object using '{path}'...")

            # Extra details
            self.z_cmb = np.nan
            self.origin = source
            self.originalname = path.split('/')[-1].split('.')[0]
            self.objname = (self.originalname.replace('CSP', '').
                            replace('ATLAS', '').
                            replace('ZTF', ''))
            self.path = f"classes/{source}/{source}_{self.originalname}_{algo}_class.txt"

            # Load and clean data

            var_tbl = self.make_objTbl(source, path)
            var_tbl = self.clean_objTbl(var_tbl)
            self.coords = [np.average(var_tbl['ra']), np.average(var_tbl['dec'])]
            # if 'z' in var_tbl.colnames:
            #     self.z = np.average(var_tbl['z'][~np.isnan(var_tbl['z'])])
            # self.get_details()
            self.check_tarlist(source)

            # # Fix offset peak light
            # if source == 'combined_91bg':
            #     var_tbl = self.adjust_time_scale(var_tbl)

            # Adjust ztf time scale
            if source == 'ztf-91bg' or source == 'combined_91bg' or source == 'combined_norm':
                var_tbl = var_tbl[(var_tbl['mjd'] > self.discdate - 30) & (var_tbl['mjd'] < self.discdate + 120)]

            # Set arrays
            self.mjd = np.array(var_tbl['mjd'])
            self.mag = np.array(var_tbl['mag'])
            self.dmag = np.array(var_tbl['dmag'])
            self.flux = np.array(var_tbl['flux'])
            self.dflux = np.array(var_tbl['dflux'])
            self.filter = np.array(var_tbl['filter'])
            self.zp = np.array(var_tbl['zp'])
            self.params = {} # Initalized for later
            self.covariance = [] # Initalized for later

            # Save class file
            self.save_class()
        else:
            print(f"[!!!] Unknown source/algo/path '{source}'/'{algo}'/'{path}'...")
        return
    def __str__(self):
        return (f"{self.objname}, {self.origin} @ {self.path}\n"
                f"{str(self.coords)} | ({self.z}, {str(round(self.z_cmb, 2))})")
    def make_objTbl(self, source: str, path: str) -> Table:
        # Load data and make intial table
        if source == 'atlas-91bg' or source == 'atlas-norm':
            # Load data
            with open(path, 'r') as f: hdr = f.readline()[1:-1].split(',')
            data = np.genfromtxt(path, dtype='str', delimiter=',', skip_header=1)

            # Make table
            var_table = Table()
            for h in hdr:
                try: var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError: var_table[h] = data[:, hdr.index(h)]

            # Add parity with CSP & ZTF
            var_table.remove_columns(['err', 'chi/N', 'x', 'y', 'maj', 'min', 'phi', 'apfit', 'mag5sig', 'Sky', 'Obs'])
            for h_old, h_new in zip(['JD', 'm', 'dm', 'uJy', 'duJy', 'F', 'RA', 'Dec'],
                                    ['mjd', 'mag', 'dmag', 'flux', 'dflux', 'filter', 'ra', 'dec']):
                var_table[h_old].name = h_new

            # Recreate Zero-Points using ZP = m + 2.5log_10(flux)
            var_table['zp'] = np.full(var_table['mag'].shape, 23.9)
        elif source == 'csp-91bg' or source == 'csp-norm':
            var_table = Table(names=['filter', 'zp', 'z', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux'],
                              dtype=[str, float, float, float, float, float, float, float, float, float])
            with open(path, 'r') as f:
                csp_objname, csp_z, csp_ra, csp_dec = f.readline()[2:-1].split(' ')
                for l in f.readlines():
                    l = l.split(' ')
                    # Filter line
                    if len(l) == 2:
                        csp_filter = str(l[1][:-1])
                        csp_zp = float(utils.get_constants()['csp_zpts_'+csp_filter])
                    else:
                        csp_mjd, csp_mag, csp_dmag = float(l[-3])+53000, float(l[-2]), float(l[-1])
                        csp_flux = 10 ** ((csp_mag - csp_zp) / -2.5)
                        csp_dflux = np.abs(csp_flux) * np.log(10) * ((1 / 2.5) * csp_dmag)

                        var_table.add_row([csp_filter, csp_zp, csp_z, csp_ra, csp_dec, csp_mjd,
                                           csp_mag, csp_dmag, csp_flux, csp_dflux])
            self.z = var_table['z'][0]  # Special case for CSP since CSP-I didn't report all SNe redshifts to TNS
        elif source == 'ztf-91bg' or source == 'ztf-norm':
            # # Load ZTF data from file
            # ztf_data = np.genfromtxt(path, delimiter=',', dtype=str)
            # ztf_hdr, ztf_data = list(ztf_data[0, :]), ztf_data[1:, :]
            # ztf_ra = ztf_data[:, ztf_hdr.index('ra')].astype(float)
            # ztf_dec = ztf_data[:, ztf_hdr.index('dec')].astype(float)
            # ztf_filter = ztf_data[:, ztf_hdr.index('filtercode')]
            # ztf_mjd = ztf_data[:, ztf_hdr.index('mjd')].astype(float)
            # ztf_mag = ztf_data[:, ztf_hdr.index('mag')].astype(float)
            # ztf_dmag = ztf_data[:, ztf_hdr.index('magerr')].astype(float)
            # ztf_zp = ztf_data[:, ztf_hdr.index('magzp')].astype(float)
            #
            # # Calculate ZTF flux values
            # ztf_flux = 10 ** ((ztf_mag.astype(float) - ztf_zp.astype(float)) / (-2.5))
            # ztf_dflux = np.abs(ztf_flux) * (1 / 2.5) * np.log(10) * ztf_dmag.astype(float)
            #
            # # Make table of combined data
            # var_table = Table(names=['ra', 'dec', 'filter', 'mjd', 'mag', 'dmag', 'flux', 'dflux', 'zp'],
            #                   data=[ztf_ra, ztf_dec, ztf_filter, ztf_mjd, ztf_mag, ztf_dmag, ztf_flux, ztf_dflux,
            #                         ztf_zp])
            # Load hdr & data
            with open(path, 'r') as f:
                for i in range(3): f.readline()
                ztf_ra = float(f.readline().split(' ')[-2])
                ztf_dec = float(f.readline().split(' ')[-2])
            data = np.genfromtxt(path, delimiter=' ', dtype=str, skip_header=54)
            hdr, data = list(data[0]), data[1:]
            for i in range(len(hdr)): hdr[i] = hdr[i][:-1]

            # Make table
            var_table = Table()
            for h in hdr:
                try:
                    var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError:
                    var_table[h] = data[:, hdr.index(h)]

            # Add RA & DEC
            var_table['ra'] = np.full(len(var_table), ztf_ra)
            var_table['dec'] = np.full(len(var_table), ztf_dec)

            # Fix time, JD to MJD
            var_table['jd'] = np.array(var_table['jd']).astype(float) - 2400000.5
            var_table['forcediffimflux'][var_table['forcediffimflux'] == 'null'] = 'nan'
            var_table['forcediffimfluxunc'][var_table['forcediffimfluxunc'] == 'null'] = 'nan'
            var_table['mag'] = (-2.5 * np.log10(np.array(var_table['forcediffimflux']).astype(float))) + var_table['zpdiff']
            var_table['dmag'] = np.abs(-1.08573620476 * (np.array(var_table['forcediffimfluxunc']).astype(float)
                                                         / np.array(var_table['forcediffimflux']).astype(float)))

            # # ATLAS Pairty flux
            # new_flux = var_table['forcediffimflux'].astype(float) * 10 ** (-0.4 * (23.9 - var_table['zpdiff'].astype(float)))
            # new_dflux = var_table['forcediffimfluxunc'].astype(float) * 10 ** (-0.4 * (23.9 - var_table['zpdiff'].astype(float)))
            # # new_dflux = np.abs((new_flux * -0.4 * (23.9-var_table['zpdiff'].astype(float)) * var_table['forcediffimfluxunc'].astype(float)) / (var_table['forcediffimflux'].astype(float)))
            # var_table['forcediffimflux'] = new_flux
            # var_table['forcediffimfluxunc'] = new_dflux

            # Add parity with CSP & ATLAS
            var_table.remove_columns(['index', 'field', 'ccdid', 'qid', 'pid', 'infobitssci', 'sciinpseeing',
                                      'scibckgnd', 'scisigpix', 'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms',
                                      'clrcoeff', 'clrcoeffunc', 'ncalmatches', 'exptime', 'adpctdif1', 'adpctdif2',
                                      'diffmaglim', 'programid', 'rfid', 'forcediffimsnr', 'forcediffimchisq',
                                      'forcediffimfluxap', 'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr',
                                      'dnearestrefsrc', 'nearestrefmag', 'nearestrefmagunc', 'nearestrefchi',
                                      'nearestrefsharp', 'refjdstart', 'refjdend', 'procstatu'])
            for h_old, h_new in zip(['zpdiff', 'jd', 'forcediffimflux', 'forcediffimfluxunc'],
                                    ['zp', 'mjd', 'flux', 'dflux']):
                var_table[h_old].name = h_new
        elif source == 'combined_91bg' or source == 'combined_norm':
            # Load data
            with open(path, 'r') as f:
                f.readline() # Skip header
                hdr = f.readline().rstrip('\n').split(', ')
            data = np.genfromtxt(path, dtype='str', delimiter=', ', skip_header=2)

            # Make table
            var_table = Table()
            for h in hdr:
                try:
                    var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError:
                    var_table[h] = data[:, hdr.index(h)]

            print(f"[~~~] Making class using combined dataset... {np.unique(np.array(var_table['source']))}")
        else:
            raise ValueError(f'[!!!] Unknown source, {source}! Must be (csp/atlas/ztf/combined)-(91bg/norm)...')
        return var_table
    def new_make_objTbl(self, source: str, path: str) -> Table:
        # Load data and make intial table
        if source == 'atlas-91bg' or source == 'atlas-norm':
            # Load data from file
            atlas_data = np.genfromtxt(path, delimiter=', ', dtype=str)
            atlas_hdr, atlas_data = list(atlas_data[0, :]), atlas_data[1:, :]
            atlas_ra = atlas_data[:, atlas_hdr.index('RA')]
            atlas_dec = atlas_data[:, atlas_hdr.index('Dec')]
            atlas_filter = atlas_data[:, atlas_hdr.index('F')]
            atlas_mjd = atlas_data[:, atlas_hdr.index('MJD')]
            atlas_mag = atlas_data[:, atlas_hdr.index('m')]
            atlas_dmag = atlas_data[:, atlas_hdr.index('dm')]
            atlas_flux = atlas_data[:, atlas_hdr.index('uJy')]
            atlas_dflux = atlas_data[:, atlas_hdr.index('duJy')]

            # Calculate zeropoints
            atlas_zp = -2.5 * np.log10(atlas_flux.astype(float)) + 23.9

            # Make ATLAS table
            var_table = Table(names=['ra', 'dec', 'filter', 'mjd', 'mag', 'dmag', 'flux', 'dflux', 'zp'],
                              data=[atlas_ra, atlas_dec, atlas_filter, atlas_mjd, atlas_mag, atlas_dmag, atlas_flux,
                                    atlas_dflux, atlas_zp],
                              dtype=[float, float, str, float, float, float, float, float, float])
            # # Load data
            # with open(path, 'r') as f: hdr = f.readline()[1:-1].split(',')
            # data = np.genfromtxt(path, dtype='str', delimiter=',', skip_header=1)
            #
            # # Make table
            # var_table = Table()
            # for h in hdr:
            #     try: var_table[h] = data[:, hdr.index(h)].astype(float)
            #     except ValueError: var_table[h] = data[:, hdr.index(h)]
            #
            # # Add parity with CSP & ZTF
            # var_table.remove_columns(['err', 'chi/N', 'x', 'y', 'maj', 'min', 'phi', 'apfit', 'mag5sig', 'Sky', 'Obs'])
            # for h_old, h_new in zip(['JD', 'm', 'dm', 'uJy', 'duJy', 'F', 'RA', 'Dec'],
            #                         ['mjd', 'mag', 'dmag', 'flux', 'dflux', 'filter', 'ra', 'dec']):
            #     var_table[h_old].name = h_new
            #
        elif source == 'csp-91bg' or source == 'csp-norm':
            var_table = Table(names=['filter', 'zp', 'z', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux'],
                              dtype=[str, float, float, float, float, float, float, float, float, float])
            with open(path, 'r') as f:
                csp_objname, csp_z, csp_ra, csp_dec = f.readline()[2:-1].split(' ')
                for l in f.readlines():
                    l = l.split(' ')
                    # Filter line
                    if len(l) == 2:
                        csp_filter = str(l[1][:-1])
                        csp_zp = float(utils.get_constants()['csp_zpts_'+csp_filter])
                    else:
                        csp_mjd, csp_mag, csp_dmag = float(l[-3])+53000, float(l[-2]), float(l[-1])
                        csp_flux = 10 ** ((csp_mag - csp_zp) / -2.5)
                        csp_dflux = np.abs(csp_flux) * np.log(10) * ((1 / 2.5) * csp_dmag)

                        var_table.add_row([csp_filter, csp_zp, csp_z, csp_ra, csp_dec, csp_mjd,
                                           csp_mag, csp_dmag, csp_flux, csp_dflux])
            self.z = var_table['z'][0]  # Special case for CSP since CSP-I didn't report all SNe redshifts to TNS
        elif source == 'ztf-91bg' or source == 'ztf-norm':
            # Load data from file
            ztf_data = np.genfromtxt(path, delimiter=', ', dtype=str)

            ztf_hdr, ztf_data = list(ztf_data[0, :]), ztf_data[1:, :]
            ztf_ra = ztf_data[:, ztf_hdr.index('ra')]
            ztf_dec = ztf_data[:, ztf_hdr.index('dec')]
            ztf_filter = ztf_data[:, ztf_hdr.index('filtercode')]
            ztf_mjd = ztf_data[:, ztf_hdr.index('mjd')]
            ztf_mag = ztf_data[:, ztf_hdr.index('mag')]
            ztf_dmag = ztf_data[:, ztf_hdr.index('magerr')]
            ztf_zp = ztf_data[:, ztf_hdr.index('magzp')]



            # Calculate flux values
            ztf_flux = 10 ** ((ztf_mag.astype(float) - ztf_zp.astype(float)) / (-2.5))
            ztf_dflux = np.abs(ztf_flux) * (1 / 2.5) * np.log(10) * ztf_dmag.astype(float)

            # Make ZTF table
            var_table = Table(names=['ra', 'dec', 'filter', 'mjd', 'mag', 'dmag', 'flux', 'dflux', 'zp'],
                              data=[ztf_ra, ztf_dec, ztf_filter, ztf_mjd, ztf_mag, ztf_dmag, ztf_flux, ztf_dflux, ztf_zp],
                              dtype=[float, float, str, float, float, float, float, float, float])
            # # Load hdr & data
            # with open(path, 'r') as f:
            #     for i in range(3): f.readline()
            #     ztf_ra = float(f.readline().split(' ')[-2])
            #     ztf_dec = float(f.readline().split(' ')[-2])
            # data = np.genfromtxt(path, delimiter=' ', dtype=str, skip_header=54)
            # hdr, data = list(data[0]), data[1:]
            # for i in range(len(hdr)): hdr[i] = hdr[i][:-1]
            #
            # # Make table
            # var_table = Table()
            # for h in hdr:
            #     try:
            #         var_table[h] = data[:, hdr.index(h)].astype(float)
            #     except ValueError:
            #         var_table[h] = data[:, hdr.index(h)]
            #
            # # Add RA & DEC
            # var_table['ra'] = np.full(len(var_table), ztf_ra)
            # var_table['dec'] = np.full(len(var_table), ztf_dec)
            #
            # # Fix time, JD to MJD
            # var_table['jd'] = np.array(var_table['jd']).astype(float) - 2400000.5
            # var_table['forcediffimflux'][var_table['forcediffimflux'] == 'null'] = 'nan'
            # var_table['forcediffimfluxunc'][var_table['forcediffimfluxunc'] == 'null'] = 'nan'
            # var_table['mag'] = (-2.5 * np.log10(np.array(var_table['forcediffimflux']).astype(float))) + var_table['zpdiff']
            # var_table['dmag'] = np.abs(-1.08573620476 * (np.array(var_table['forcediffimfluxunc']).astype(float)
            #                                              / np.array(var_table['forcediffimflux']).astype(float)))
            #
            # # Add parity with CSP & ATLAS
            # var_table.remove_columns(['index', 'field', 'ccdid', 'qid', 'pid', 'infobitssci', 'sciinpseeing',
            #                           'scibckgnd', 'scisigpix', 'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms',
            #                           'clrcoeff', 'clrcoeffunc', 'ncalmatches', 'exptime', 'adpctdif1', 'adpctdif2',
            #                           'diffmaglim', 'programid', 'rfid', 'forcediffimsnr', 'forcediffimchisq',
            #                           'forcediffimfluxap', 'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr',
            #                           'dnearestrefsrc', 'nearestrefmag', 'nearestrefmagunc', 'nearestrefchi',
            #                           'nearestrefsharp', 'refjdstart', 'refjdend', 'procstatu'])
            # for h_old, h_new in zip(['zpdiff', 'jd', 'forcediffimflux', 'forcediffimfluxunc'],
            #                         ['zp', 'mjd', 'flux', 'dflux']):
            #     var_table[h_old].name = h_new
        elif source == 'combined_91bg' or source == 'combined_norm':
            # Load data
            with open(path, 'r') as f:
                f.readline() # Skip header
                hdr = f.readline().rstrip('\n').split(', ')
            data = np.genfromtxt(path, dtype='str', delimiter=', ', skip_header=2)

            # Make table
            var_table = Table()
            for h in hdr:
                try:
                    var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError:
                    var_table[h] = data[:, hdr.index(h)]

            print(f"[~~~] Making class using combined dataset... {np.unique(np.array(var_table['source']))}")
        else:
            raise ValueError(f'[!!!] Unknown source, {source}! Must be (csp/atlas/ztf/combined)-(91bg/norm)...')
        return var_table
    def clean_objTbl(self, tbl: Table) -> Table:
        # Remove '>' from ATLAS mags
        for i, n_mag in enumerate(tbl['mag']):
            if str(n_mag)[0] == '>': tbl['mag'][i] = float(n_mag[1:])

        # Remove nulls
        for col in ['mag', 'dmag', 'flux', 'dflux']:
            tbl = tbl[~np.isnan(np.array(tbl[col]).astype(float))]

        # Remove non-detections (negative/zero magnitudes/fluxes)
        tbl = tbl[tbl['mag'].astype(float) > 0]
        tbl = tbl[tbl['flux'].astype(float) > 0]

        # Remove unresonable dmag/dflux
        for col, lim in zip(['dmag', 'dflux'],
                            [float(utils.get_constants()['mag_err_lim']),
                             float(utils.get_constants()['flux_err_lim'])]):
            if lim > 0:
                tbl = tbl[tbl[col].astype(float) < lim]
                tbl = tbl[tbl[col].astype(float) != 0]

        return tbl
    def get_details(self, search_radius: float = 2, check_tnskey: bool = True):
        # Check TNS key first
        if check_tnskey:
            TNS_key_tb = utils.check_tnskey(self.coords[0], self.coords[1])
            if TNS_key_tb is not None and len(TNS_key_tb) == 1:
                print(f"[+++] Info found in TNS key! Pulling...")
                self.objname = TNS_key_tb['objname'][0]
                self.z = TNS_key_tb['z'][0]
                self.coords = [TNS_key_tb['ra'][0], TNS_key_tb['dec'][0]]
                self.discdate = TNS_key_tb['discoverydate'][0]
                return

        # Query TNS parameters
        objname = (self.originalname.
                   replace('CSP', '').
                   replace('ATLAS', '').
                   replace('ZTF', ''))
        self.objname = objname
        self.discdate = -99999

        # Query using object name
        print(f"[+++] Querying TNS @ {objname}, with radius={search_radius}...")
        for t in range(10):
            try:
                details = utils.query_tns(objname=self.objname, search_radius=search_radius)
                break
            except Exception as e:
                print(f"[~~~] Warning: TNS found no targets @ ({objname}), with radius={search_radius}... "
                      f"Raising to search_radius={search_radius + 1}")
                search_radius += 1

        # If objname fails, use ra & dec
        if self.discdate < 0:
            print(f"[+++] Querying TNS @ RA, DEC = ({self.coords[0]}, {self.coords[1]}), with radius={search_radius}...")
            for t in range(10):
                try:
                    details = utils.query_tns(coords=[self.coords[0], self.coords[1]], search_radius=search_radius)
                    break
                except Exception as e:
                    print(f"[~~~] Warning: TNS found no targets @ ({self.coords[0]}, {self.coords[1]}), with radius={search_radius}... "
                          f"Raising to search_radius={search_radius + 1}")
                    search_radius += 1

        # Save TNS details
        ## Special case for CSP since CSP-I didn't report all SNe redshifts to TNS
        if details['redshift'] is None or self.origin.split('-')[0] == 'csp':
            print(f'[~~~] Speical CSP case: Keeping CSP heliocentric redshift {self.z}')
            self.z = self.z
        else:
            self.z = details['redshift']
        self.discdate = float(astrotime(details['discoverydate'], format='iso').jd) - 2400000.5  # JD to MJD
        return

        # success = False
        # while not success:
        #     try:
        #         details = utils.query_tns(objname=self.objname, search_radius=search_radius)
        #
        #         # Special case for CSP since CSP-I didn't report all SNe redshifts to TNS
        #         if details['redshift'] is None and self.origin.split('-')[0] == 'csp':
        #             print(f'[~~~] Speical CSP case: Keeping CSP heliocentric redshift {self.z}')
        #             self.z = self.z
        #         else:
        #             self.z = details['redshift']
        #
        #         self.discdate = float(astrotime(details['discoverydate'], format='iso').jd) - 2400000.5  # JD to MJD
        #         success = True
        #     except:
        #         print(f"[~~~] Warning: TNS found no targets @ ({objname}), with radius={r}... Raising to r={r+1}")
        #         r += 1
        #     if r == 10:
        #         # raise RuntimeError('[!!!!!] TNS could not find transient!')
        #         print("[!!!!!] TNS could not find transient!")
        #         self.z_cmb = np.nan
        #         self.objname = self.originalname
        #         self.z = -99999
        #         self.discdate = -99999
        #         return None
        # print(f"[+++] Querying TNS @ ({ra}, {dec}), with radius={r}...")
        # self.coords[0], self.coords[1]
        # # r = utils.query_tns(coords=["174.61858", "20.52622"])
        # # r = utils.query_tns(objname="2006bd")
        # objname = (self.originalname.
        #            replace('CSP', '').
        #            replace('ATLAS', '').
        #            replace('ZTF', ''))
        #
        #
        # success = False
        # while not success:
        #     try:
        #         details = get_TNSDetails(ra, dec, radius=r)
        #
        #         # Special case for CSP since CSP-I didn't report all SNe redshifts to TNS
        #         if details['redshift'] is None and self.origin.split('-')[0] == 'csp':
        #             print(f'[~~~] Speical CSP case: Keeping CSP heliocentric redshift {self.z}')
        #             self.z = self.z
        #         else:
        #             self.z = details['redshift']
        #
        #         self.objname = details['objname']
        #         self.z_cmb = np.nan
        #         self.coords = [details['radeg'], details['decdeg']]
        #         self.discdate = float(astrotime(details['discoverydate'], format='iso').jd) - 2400000.5  # JD to MJD
        #         success = True
        #     except:
        #         print(f"[~~~] Warning: TNS found no targets @ ({ra}, {dec}), with radius={r}... Raising to r={r+1}")
        #         r += 1
        #     if r == 10:
        #         # raise RuntimeError('[!!!!!] TNS could not find transient!')
        #         print("[!!!!!] TNS could not find transient!")
        #         self.z_cmb = np.nan
        #         self.objname = self.originalname
        #         self.z = -99999
        #         self.discdate = -99999
        #         return None
    def check_tarlist(self, source: str):
        # Select target file based on source
        if '91bg' in source:
            tarlist_path = 'txts/target_files/sn1991bglike_tarlist.csv'
        elif 'norm' in source:
            tarlist_path = 'txts/target_files/normal_tarlist.csv'

        # Load data
        data = np.genfromtxt(tarlist_path, delimiter=',', dtype=str, skip_header=1)
        tb = Table(names=data[0, :], data=data[1:, :])
        selc_tb = tb[tb['Name'] == f"SN {self.objname}"].copy()

        if len(selc_tb) >= 1:
            # Set redshift, convert and set discovery date
            self.z = float(selc_tb['Redshift'][0])

            if selc_tb['Discovery Date (UT)'][0] == "Unknown":
                self.discdate = None
                return
            else:
                discdata_utc = selc_tb['Discovery Date (UT)'][0].split(' ')[0].split('-')
                discdata_utc_datetime = datetime.datetime(int(discdata_utc[0]), int(discdata_utc[1]), int(discdata_utc[2]))
                discdata_mjd = astrotime(discdata_utc_datetime, format='datetime', scale='utc').mjd
                self.discdate = float(discdata_mjd)
        else:
            self.z = -99999
            self.discdate = -99999

        return
    def adjust_time_scale(self, tbl: Table):
        # # Find difference between ATLAS and ZTF peaks
        # atlas_peaks, ztf_peaks = [], []
        # for peaks, s_filters in zip([atlas_peaks, ztf_peaks],
        #                             [['o', 'c'], ['ZTF_g', 'ZTF_r', 'ZTF_i']]):
        #     for f in s_filters:
        #         if f in tbl['filter'] and len(tbl[tbl['filter'] == f]) > 3:
        #             f_tbl = tbl[tbl['filter'] == f]
        #             f_tbl_fmax = np.max(f_tbl['flux'])
        #             f_tbl_tmax = f_tbl['mjd'][f_tbl['flux'] == f_tbl_fmax][0]
        #             peaks.append(f_tbl_tmax)
        # if len(atlas_peaks) > 0 or len(ztf_peaks) > 0:
        #     # Calculate difference
        #     atlas_peak = np.average(atlas_peaks)
        #     ztf_peak = np.average(ztf_peaks)
        #     peak_diff = ztf_peak - atlas_peak
        #
        #     # Adjust ATLAS for difference
        #     for f in ['o', 'c']:
        #         if f in tbl['filter'] and len(tbl[tbl['filter'] == f]) > 3:
        #             tbl['mjd'][tbl['filter'] == f] = tbl['mjd'][tbl['filter'] == f] + peak_diff

        # # Adjust bounds of ZTF data
        # tbl = tbl[(tbl['mjd'] > self.discdate-15) & (tbl['mjd'] < self.discdate+60)]

        return tbl
    def save_class(self):
        print(f"[+++] Saving {self.objname} class to {self.path}...")
        with open(self.path, 'w') as f:
            f.write(f"{self.origin},{self.objname},{self.originalname},{self.coords[0]},{self.coords[1]},{self.z},{self.z_cmb},{self.discdate}\n")
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for p in self.params:
                f.write(f"{p},{self.params[p]['value']},{self.params[p]['err']}\n")
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for line in self.covariance:
                p_line = str(line[0])
                for i in range(1, len(line)):
                    p_line += ',' + str(line[i])
                f.write(p_line + '\n')
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for col, name in zip([self.zp, self.filter, self.mjd, self.flux, self.dflux, self.mag, self.dmag],
                                 ['zp', 'filter', 'mjd', 'flux', 'dflux', 'mag', 'dmag']):
                f.write(f"{name}")
                for n_col in col:
                    f.write(f",{n_col}")
                f.write('\n')
    def load_class(self, path: str):
        print(f"[+++] Opening class from {path}...")
        with open(path, 'r') as f:
            # Read header
            details = f.readline().split(',')
            self.path = path
            self.origin = details[0]
            self.objname = details[1]
            self.originalname = details[2]
            self.coords = [float(details[3]), float(details[4])]
            self.z = None if details[5] == "None" else float(details[5])
            self.z_cmb = float(details[6])
            self.discdate = float(details[7][:-1])
            f.readline()  # Skip break line

            # Read params
            self.params = {}
            line = f.readline()
            while '+++' not in line:
                line = line.split(',')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            # Read covariances
            self.covariance = np.array([])
            if 'salt' in path:
                line = f.readline()
                while '+++' not in line:
                    line = line.split(',')
                    line[-1] = line[-1][:-1]
                    self.covariance = np.append(self.covariance, np.array(line).astype(float))
                    line = f.readline()
                if len(self.covariance) > 0:
                    self.covariance = self.covariance.reshape(4, 4)
            else:
                f.readline()  # Skip break line

            # Read arrays
            try:
                self.zp = np.array(f.readline().split(',')[1:])
                self.zp[-1] = self.zp[-1][:-1]
                self.filter = np.array(f.readline().split(',')[1:])
                self.filter[-1] = self.filter[-1][:-1]
                self.mjd = np.array(f.readline().split(',')[1:])
                self.mjd[-1] = self.mjd[-1][:-1]
                self.flux = np.array(f.readline().split(',')[1:])
                self.flux[-1] = self.flux[-1][:-1]
                self.dflux = np.array(f.readline().split(',')[1:])
                self.dflux[-1] = self.dflux[-1][:-1]
                self.mag = np.array(f.readline().split(',')[1:])
                self.mag[-1] = self.mag[-1][:-1]
                self.dmag = np.array(f.readline().split(',')[1:])
                self.dmag[-1] = self.dmag[-1][:-1]
            except IndexError:
                print(f"[~~~] While loading {self.objname}, empty arrays were found!")
                self.zp = np.array([])
                self.filter = np.array([])
                self.mjd = np.array([])
                self.flux = np.array([])
                self.dflux = np.array([])
                self.mag = np.array([])
                self.dmag = np.array([])
            except Exception as e:
                print(e)
        return

    # Display Functions ---------------------------------------------------------------------------------------------- #
    def simple_plot(self):
        for f in np.unique(self.filter):
            if f in ['J', 'H', 'Y', 'Ydw', 'u']: continue
            plt.scatter(self.mjd[self.filter == f].astype(float), self.flux[self.filter == f].astype(float), label=f)
        plt.axvline(self.discdate)
        plt.legend()
        plt.show()
        return
    def plot(self, y_type: str = 'mag', fmt: str = 'o', save_loc: str = ''):
        print(f"[+++] Plotting LC of '{self.objname}' ({y_type})...")
        filter_dict = {'u': 'teal', 'g': 'green', 'r': 'red', 'i': 'indigo', 'B': 'blue',
                       'V0': 'violet', 'V1': 'purple', 'V': 'red', 'Y': 'goldenrod', 'Hdw': 'tomato', 'H': 'salmon',
                       'J': 'aquamarine', 'Jrc2': 'cadetblue', 'Jdw': 'turquoise', 'Ydw': 'olive',
                       'c': 'cyan', 'o': 'orange', 'ZTF_g': 'green', 'ZTF_r': 'red', 'ZTF_i': 'indigo',
                       'zg': 'green', 'zr': 'red', 'zi': 'indigo'}

        # # Select y-axis element
        # y_axis = self.mag if y_type == 'mag' else self.flux
        # y_axis_err = self.dmag if y_type == 'mag' else self.dflux
        #
        # # Plot
        # fig, axs = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
        # for f in np.unique(self.filter):
        #     axs.errorbar(self.mjd[self.filter == f].astype(float), y_axis[self.filter == f].astype(float),
        #                  yerr=y_axis_err[self.filter == f].astype(float),
        #                  color = filter_dict[f], label=f, fmt=fmt, ms=3)

        # Plot
        fig, axs = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
        for f in np.unique(self.filter):
            if y_type == 'mag':
                axs.errorbar(self.mjd[self.filter == f].astype(float),
                             self.mag[self.filter == f].astype(float) - abs(self.mag[self.filter == f].astype(float) - self.zp[self.filter == f].astype(float)),
                             yerr=self.dmag[self.filter == f].astype(float),
                             color=filter_dict[f], label=f, fmt=fmt, ms=3)
            elif y_type == 'flux':
                axs.errorbar(self.mjd[self.filter == f].astype(float),
                             self.flux[self.filter == f].astype(float),
                             yerr=self.dflux[self.filter == f].astype(float),
                             color=filter_dict[f], label=f, fmt=fmt, ms=3)

        # Lines
        axs.axvline(self.discdate, color='black', linestyle='--', label='Discovery Date')
        try: axs.axvline(self.params['Tmax']['value'], color='maroon', linestyle='--', label='Peak Brightness')
        except: pass

        # Formatting
        if y_type == 'mag': axs.invert_yaxis()
        axs.legend()
        axs.set_xlabel('MJD', size=16)
        axs.set_ylabel('Magnitude', size=16) if y_type == 'mag' else axs.set_ylabel('Flux (uJy)', size=16)
        plt.suptitle(f"Lightcurve of '{self.objname}' ({y_type})", size=16)
        if len(save_loc) != 0:
            print('[+++] Saving to '+save_loc)
            plt.savefig(save_loc)
        plt.show()
        return

    # Fitting Functions ---------------------------------------------------------------------------------------------- #
    def write_snpy_ascii(self, save_loc: str):
        filter_dict = {'o': 'ATri', 'c': 'ATgr', 't': 'ATri2',
                       'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i',
                       'zg': 'g', 'zr': 'r', 'zi': 'i',
                       'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'V0': 'V0', 'Y': 'Y', 'Ydw': 'Ydw',
                       'g': 'g', 'i': 'i', 'r': 'r', 'u': 'u'}
        with open(save_loc, 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(f"{self.objname} {self.z} {self.coords[0]} {self.coords[1]}\n")
            for f_w in np.unique(self.filter):
                f_indexes = np.where(self.filter == f_w)[0]
                f.write(f"filter {filter_dict[f_w]}\n")
                for i in f_indexes:
                    # filter photometry block -- Date (JD/MJD), mag, err (i.e. 674.8593 12.94 0.11)
                    f.write(f"{self.mjd[i]} {self.mag[i]} {self.dmag[i]}\n")
        print(f'[+++] Saved file to... {save_loc}')
        return
    def snpy_fit(self):
        ascii_path = 'fitting/snpy-asciis/' + self.objname + '_ascii.txt'
        model_path = 'fitting/snpy-models/' + self.objname + '_model.txt'
        plot_path = 'fitting/snpy-plots/' + self.objname + '_lc.png'
        param_names = ['mu', 'st', 'Tmax', 'EBVhost']
        snpy_param_names = ['DM', 'st', 'Tmax', 'EBVhost']
        show_plots = True

        # Make ascii file for SNooPy to read
        self.write_snpy_ascii(ascii_path)

        # Load Data
        try:
            n_s = snoopyfit.get_sn(ascii_path)
        except Exception as error:
            self.params.update({'mu': {'value': 0.00, 'err': 0.00}})
            print('[!!!] Failed to load ASCII file -- ', error)
            return
        n_s.choose_model('EBV_model2', stype='st')
        n_s.set_restbands()  # Auto pick appropriate rest-bands

        # Remove empty filters -- fix for 'ValueError: attempt to get argmin of an empty sequence'
        for class_filter in list(n_s.data.keys()):
            if len(n_s.data[class_filter].magnitude) == 0:
                del n_s.data[class_filter]
            elif self.origin.split('-')[0] == 'csp' and class_filter in ['u', 'Y', 'J', 'H', 'Jrc2', 'Ydw']:
                print('[~~~] Speical CSP case: Removing ' + class_filter + '...')
                del n_s.data[class_filter]

        # Fit with SNooPy -- gives 20 tries before failing
        for i in range(20):
            print(f"[{i+1}/20] Attempting to fit '{self.objname}' with {list(n_s.data.keys())} =======================")
            try:
                ##### I care recall why we did this in the first place ####
                # if self.origin.split('-')[0] == 'csp':
                #     initial_filters = []
                #     for fil in ['B', 'V', 'g']:
                #         if fil in list(n_s.data.keys()):
                #             initial_filters.append(fil)
                #     print(f'[~~~] Speical CSP case: Fitting as {initial_filters} -> remaining...')
                #
                #     n_s.fit(initial_filters, dokcorr=True, k_stretch=False, reset_kcorrs=True,
                #             **{'mangle': 1, 'calibration': 0})
                #     n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True,
                #             **{'mangle': 1, 'calibration': 0})
                # else:
                n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
                n_s.save(model_path)

                # Save parameters
                for j in range(len(param_names)):
                    self.params.update({param_names[j]: {'value': n_s.parameters[snpy_param_names[j]],
                                                         'err': n_s.errors[snpy_param_names[j]]}})
                self.params.update({'chisquare': {'value': n_s.model.chisquare,
                                                  'err': n_s.model.rchisquare}})
                n_s.plot(outfile=plot_path)
                if show_plots:
                    plt.show()
                    systime.sleep(2)
                plt.close()
                print(f'[+++] Successfully fit {self.objname}!')
                break
            except Exception as error:
                if 'All weights for filter' and 'are zero.' in str(error):
                    print('[!!!] Weights for filter', str(error).split(' ')[4], 'are zero. Removing...')
                    del n_s.data[str(error).split(' ')[4]]
                elif str(error) == 'Error:  to solve for EBVhost, you need to fit more than one filter':
                    print('[!!!] To few filters to fit!')
                    print(f'[---] Could not fit {self.objname}!')
                    self.params.update({'mu': {'value': -124.0, 'err': -124.0}})
                    break
                else:
                    self.params.update({'mu': {'value': -404.0, 'err': -404.0}})
                    print(error)
                    print(f'[---] Could not fit {self.objname}!')
                    break
        self.save_class()
        return
    def salt_fit(self):
        show_plot = True
        plot_path = f"fitting/salt-plots/{self.objname}_lc.png"

        CONSTANTS = utils.get_constants()
        alpha, beta = float(CONSTANTS['salt_alpha']), float(CONSTANTS['salt_beta'])
        mB_const, M0 = float(CONSTANTS['salt_mB_const']), float(CONSTANTS['salt_absolute_mag'])
        try:
            # Make sure more than one filter
            if len(np.unique(self.filter)) < 1:
                self.params.update({'mu': {'value': -124.0, 'err': -124.0}})
                self.save_class()
                print(f"[!!!] Too few filters! Can not fit {self.objname}! {np.unique(self.filter)}")
                return

            # Fix filters
            filter_dict = {'u': 'cspu', 'g': 'cspg', 'r': 'cspr', 'i': 'cspi', 'B': 'cspB',
                           'V0': 'cspv3014', 'V1': 'cspv3009', 'V': 'cspv9844', 'Y': 'cspys',
                           'J': 'cspjs', 'Jrc2': 'cspjd', 'Jdw': 'cspjd', 'Ydw': 'cspyd', 'Hdw': 'csphd', 'H': 'csphs',
                           'c': 'atlasc', 'o': 'atlaso', 'ZTF_g': 'ztfg', 'ZTF_r': 'ztfr', 'ZTF_i': 'ztfi',
                           'zg': 'ztfg', 'zr': 'ztfr', 'zi': 'ztfi'}
            salt_time, salt_filters, salt_flux = np.array([]), np.array([]), np.array([])
            salt_dflux, salt_zp = np.array([]), np.array([])

            for i in range(len(self.filter)):
                if self.origin.split('-')[0] == 'csp' and self.filter[i] in ['u', 'Y', 'J', 'H', 'Jrc2', 'Ydw']:
                    continue
                salt_time = np.append(salt_time, self.mjd[i])
                salt_filters = np.append(salt_filters, filter_dict[self.filter[i]])
                salt_flux = np.append(salt_flux, self.flux[i])
                salt_dflux = np.append(salt_dflux, self.dflux[i])
                salt_zp = np.append(salt_zp, self.zp[i])
            print('[~~~]', np.unique(self.filter), '->', np.unique(salt_filters))

            data = Table([salt_time, salt_filters, salt_flux, salt_dflux, salt_zp, np.full(len(salt_time), 'ab')],
                         names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'))

            if len(np.unique(data['band'])) < 2:
                print(f"[!!!] Not enough filters to fit! ({len(np.unique(data['band']))})")
                self.params.update({'mu': {'value': -124.0, 'err': -124.0}})
                return

            # Create and fit data to model
            model = salt3fit.Model(source='salt3')
            # model.set(z=self.z)  # set the model's redshift.
            model.set(z=self.z, t0=self.discdate)  # set the model's redshift.
            result, fitted_model = salt3fit.fit_lc(data, model,
                                                   ['t0', 'x0', 'x1', 'c'],
                                                   bounds={'x1': (-5, 5)},
                                                   minsnr=1)

            # Save Parameters
            param_names = ['t0', 'x0', 'x1', 'c']
            for i in range(len(param_names)):
                self.params.update({param_names[i]: {'value': result.parameters[i+1],
                                                     'err': result.errors[param_names[i]]}})

            # Save Covariance
            if result['covariance'] is not None:
                self.covariance = result['covariance']

            # Calculate
            pho_mB = -2.5 * np.log10(self.params['x0']['value']) + mB_const
            pho_mB_err = np.abs(-2.5 * (self.params['x0']['err'] / (self.params['x0']['value'] * np.log(10))))

            mu = pho_mB + (alpha * self.params['x1']['value']) - (beta * self.params['c']['value']) - M0
            mu_err = np.sqrt(pho_mB_err ** 2 + (np.abs(alpha) * self.params['x1']['err']) ** 2 + (np.abs(beta) * self.params['c']['err']) ** 2)

            self.params.update({'mu': {'value': mu, 'err': mu_err}})

            # Plot data with fit
            salt3fit.plot_lc(data, model=fitted_model, errors=result.errors)
            plt.savefig(plot_path)
            if show_plot:
                plt.show()
                systime.sleep(2)
            plt.close()

            self.save_class()
            print(f'[+++] Successfully fit {self.objname}!')
        except Exception as error:
            print(f"[!!!] Error: {str(error)}")
            if 'result is NaN for' in str(error):
                print(f"[!!!] SALT3 couldn't fit with current parameter selection! Returning nan for {self.objname}...")
                self.params.update({'mu': {'value': -107.0, 'err': -107.0}})
            elif 'No data points with S/N > 5.0. Initial guessing failed.' in str(error):
                self.params.update({'mu': {'value': -434.0, 'err': -434.0}})
            else:
                print(error)
                print(f'[---] Could not fit {self.objname}!')
                self.params.update({'mu': {'value': -404.0, 'err': -404.0}})
            self.save_class()
            return
        return
    def get_zcmb(self):
        # Set up variables
        CONSTANTS = utils.get_constants()
        local_coords = SkyCoord(self.coords[0], self.coords[1], unit="deg")
        self.z = float(self.z)

        # Get CMB redshift
        galac_coords = local_coords.transform_to(Galactic())
        helio_corr = (float(CONSTANTS['cmb_v_helio']) / float(CONSTANTS['cmb_c']) *
                      ((np.sin(galac_coords.b.deg) * np.sin(float(CONSTANTS['cmb_b_h'])) + np.cos(
                          galac_coords.b.deg) *
                        np.cos(float(CONSTANTS['cmb_b_h'])) * np.cos(
                                  galac_coords.l.deg - float(CONSTANTS['cmb_l_h'])))))
        corr_term = 1 - helio_corr
        self.z_cmb = (1 + self.z) / corr_term - 1

        # Peculiar Velocity Correction -- using 'getVpec.py' from David
        VP = getVpec.VelocityCorrection(f"twomass++_velocity_LH11.npy")
        self.z_cmb = VP.correct_redshift(self.z, 0, local_coords.galactic.l.deg, local_coords.galactic.b.deg)
        vpec, vpec_sys = getVpec.main(self.coords[0], self.coords[1], self.z_cmb)
        self.z_cmb += vpec / 3e5
        return
    def get_hostMass(self):
        # Check Mass Key
        if utils.get_constants()['use_mass_key'] == "True":
            hostMass, hostMassErr = utils.check_mass_key(objname=self.objname, mode='read')
            if ~np.isnan(hostMass) and ~np.isnan(hostMass):
                print("[+++] Host mass found in key! Pulling...")
                if 'hostMass' in self.params.keys():
                    self.params['hostMass']['value'], self.params['hostMass']['err'] = hostMass, hostMassErr
                else:
                    self.params.update({'hostMass': {'value': hostMass, 'err': hostMassErr}})
                return

        # Get GHOST database if not installed
        if os.path.exists('fitting/ghost-database/database/GHOST.csv') == False:
            getGHOST(real=False, verbose=False, installpath='fitting/ghost-database')

        # Run GHOST Transient Algorithm
        hosts = getTransientHosts([f'SN {self.objname}'],
                                  [SkyCoord(self.coords[0] * u.deg, self.coords[1] * u.deg, frame='icrs')],
                                  verbose=False, savepath='fitting/ghost-transients/',
                                  GHOSTpath='fitting/ghost-database')
        if len(hosts) == 0:
            print(f"[+++] GHOST was unable to find {self.objname}! Returning null...")
            hostMass, hostMassErr = -333, -333
        else:
            # Get magnitudes from GHOST
            print('[+++] Identified Host Galaxy:', hosts.loc[0, 'NED_name'])
            if ~np.isnan(hosts.loc[0, 'gKronMag']) and ~np.isnan(hosts.loc[0, 'iKronMag']):
                print('[+++] Calculating mass using gKronMag & iKronMag from GHOST!')
                # Pull values
                gMag, iMag, iAbsMag = (hosts['gKronMag'].loc[0], hosts['iKronMag'].loc[0],
                                       hosts['iKronMag'].loc[0] - utils.current_cosmo().distmod(self.z).value)
                gMagErr, iMagErr, iAbsMagErr = (hosts['gKronMagErr'].loc[0],
                                                hosts['iKronMagErr'].loc[0], hosts['iKronMagErr'].loc[0])
                # Mass Calculation -- Taylor et al. 2011 -- eq. 8
                hostMass = (1.15 + (0.7 * (gMag - iMag)) - (0.4 * iAbsMag))
                # Error Propagation
                giMagErr = np.sqrt((gMagErr ** 2) + (iMagErr ** 2))
                hostMassErr = np.sqrt(((0.7 ** 2) * (giMagErr ** 2)) + ((0.4 ** 2) * (iAbsMagErr ** 2)))
                print(f"[+++] The host galaxy of {self.objname} has a stellar mass of: {hostMass} +/- {hostMassErr} "
                      f"log10(M*/M_sun)])")
            elif np.isnan(hosts.loc[0, 'gKronMag']):
                print("[!!!] GHOST 'gKronMag' returned null! Unable to calculate mass...")
                hostMass, hostMassErr = -333, -333
            elif np.isnan(hosts.loc[0, 'iKronMag']):
                print("[!!!] GHOST 'iKronMag' returned null! Unable to calculate mass...")
                hostMass, hostMassErr = -333, -333

        # If GHOST fails, look manually through SDSS
        if hostMass == -333 or hostMass == np.nan:
            print(f"[~~~] Could not find mass using GHOST, searching SDSS...")
            result = SDSS.query_crossid(SkyCoord(self.coords[0] * u.deg, self.coords[1] * u.deg, frame='icrs'),
                                        photoobj_fields=['modelMag_g', 'modelMagErr_g', 'modelMag_i', 'modelMagErr_i'])
            if result is None:
                print('[!!!] GHOST & SDSS failed to find mass, returning null mass...')
                hostMass, hostMassErr = -333, -333
            else:
                if ~np.isnan(result['modelMag_g'].value[0]) and ~np.isnan(result['modelMag_i'].value[0]):
                    # Pull values
                    gMag, iMag, iAbsMag = (result['modelMag_g'].value[0], result['modelMag_i'].value[0],
                                           result['modelMag_i'].value[0] - utils.current_cosmo().distmod(self.z).value)
                    gMagErr, iMagErr, iAbsMagErr = (result['modelMagErr_g'].value[0],
                                                    result['modelMagErr_i'].value[0], result['modelMagErr_i'].value[0])
                    # Mass Calculation -- Taylor et al. 2011 -- eq. 8
                    hostMass = (1.15 + (0.7 * (gMag - iMag)) - (0.4 * iAbsMag))
                    # Error Propagation
                    giMagErr = np.sqrt((gMagErr ** 2) + (iMagErr ** 2))
                    hostMassErr = np.sqrt(((0.7 ** 2) * (giMagErr ** 2)) + ((0.4 ** 2) * (iAbsMagErr ** 2)))
                    print(
                        f"[+++] The host galaxy of {self.objname} has a stellar mass of: {hostMass} +/- {hostMassErr} "
                        f"log10(M*/M_sun)])")
                elif np.isnan(result['modelMag_g'].value[0]):
                    print("[!!!] SDSS 'modelMag_g' returned null! Unable to calculate mass...")
                    hostMass, hostMassErr = -333, -333
                elif np.isnan(result['modelMag_i'].value[0]):
                    print("[!!!] SDSS 'modelMag_i' returned null! Unable to calculate mass...")
                    hostMass, hostMassErr = -333, -333

        # Save Mass
        if 'hostMass' in self.params.keys():
            self.params['hostMass']['value'], self.params['hostMass']['err'] = hostMass, hostMassErr
        else:
            self.params.update({'hostMass': {'value': hostMass, 'err': hostMassErr}})

        # Saving mass to key
        if utils.get_constants()['use_mass_key'] == "True":
            utils.check_mass_key(objname=self.objname, mode='write', hostMass=hostMass, hostMassErr=hostMassErr)

        return
def get_TNSDetails(ra: str, dec: str, radius: str = '2'):
    """
    :param ra: Right Ascension in degrees
    :param dec: Declination in degrees
    :param radius: Radius to search out for in arcseconds
    :return: dict of TNS details
    """
    APIKEY = utils.get_apikeys()
    tns_bot_id, tns_bot_name, tns_bot_api_key = APIKEY['tns_bot_id'], APIKEY['tns_bot_name'], APIKEY['tns_bot_api_key']
    tns_marker = (
        f'tns_marker{{"tns_id": "{int(tns_bot_id)}",'
        f'"type": "bot", "name": "{tns_bot_name}"}}'
    )
    headers = {"User-Agent": tns_marker}
    search_obj = [
        ("ra", str(ra)),
        ("dec", str(dec)),
        ("radius", str(radius)),
        ("units", "arcsec"),
        ("objname", ""),
        ("objname_exact_match", 0),
        ("internal_name", ""),
        ("internal_name_exact_match ", 0),
        ("objid", ""),
        ("public_timestamp", ""),
    ]
    search_data = {"api_key": tns_bot_api_key, "data": json.dumps(OrderedDict(search_obj))}
    response = requests.post("https://www.wis-tns.org/api/get/search", headers=headers, data=search_data)
    response = json.loads(response.text)
    transients =  response["data"]
    if transients == 'Too Many Requests':
        print('\nToo Many Requests... pausing for 30 seconds: ', end='')
        for i in range(10):
            print(' - ', end='')
            time.sleep(3)
        print('\n')
        transients = response["data"]
    get_obj = [
        ("objname", transients[0]["objname"]),
        ("objid", transients[0]["objid"]),
        ("photometry", "0"),
        ("spectra", "0"),
    ]
    get_data = {"api_key": tns_bot_api_key, "data": json.dumps(OrderedDict(get_obj))}
    response = requests.post("https://www.wis-tns.org/api/get/object", headers=headers, data=get_data)
    response = json.loads(response.text)
    details = response["data"]
    return details
def error_reporting(sne: list[sneObj], source: str, print_out: bool = True, make_report: bool = True):
    """
    :param sne: list of sneObjs
    :param print_out: whether or not to print out error messages
    :return: returns list of sneObjs with no errors
    """
    # Sort error codes and behaved SNe
    valid_sne = []
    fail_reasons = {
        "TNS Faliure! Needs manual TNS": [],
        "Not enough filters to fit!": [],
        "SALT3 couldn't fit with current parameter selection!": [],
        "SALT3: No data points with S/N > 5.0!": [],
        "GHOST failed to intiate!": [],
        "GHOST failed to find mass!": [],
        "Unknown fitting error!": []
    }
    for sn in sne:
        if sn.discdate == -99999:  # TNS Faliure
            fail_reasons["TNS Faliure! Needs manual TNS"].append(sn.objname)
        elif sn.params['mu']['value'] == -124.0:
            fail_reasons["Not enough filters to fit!"].append(sn.objname)
        elif sn.params['mu']['value'] == -107.0:
            fail_reasons["SALT3 couldn't fit with current parameter selection!"].append(sn.objname)
        elif sn.params['mu']['value'] == -434.0:
            fail_reasons["SALT3: No data points with S/N > 5.0!"].append(sn.objname)
        elif 'hostMass' not in sn.params:
            fail_reasons["GHOST failed to intiate!"].append(sn.objname)
        elif np.isnan(sn.params['hostMass']['value']) or sn.params['hostMass']['value'] == -333:
            fail_reasons["GHOST failed to find mass!"].append(sn.objname)
        elif sn.params['mu']['value'] < 0:
            fail_reasons["Unknown fitting error!"].append(sn.objname)
        else:
            valid_sne.append(sn)

    # # Print out SNe and erros
    # if print_out and len(valid_sne) != len(sne):
    #     readme_mode = True
    #     print("[~~~] The following SNe failed fitting/hostMass for the following reasons...")
    #     for reason in fail_reasons.keys():
    #         print(f"{reason}: ", end='')
    #         for i, n in enumerate(fail_reasons[reason]):
    #             if readme_mode and (i + 1) % 5 == 0:
    #                 print(f"SN{n}", end=',<br/> ')
    #             else:
    #                 print(f"SN{n}", end=', ')
    #         print('\n')

    # Update error report txt
    if make_report:
        error_codes, report_tbl = error_report_builder(fail_reasons, source)

    return valid_sne
def error_report_builder(fail_reasons: dict, source: str, report_loc: str = 'txts/error_report.txt'):
    error_codes = {
        "TNS Faliure! Needs manual TNS": "99999",
        "Not enough filters to fit!": "124",
        "SALT3 couldn't fit with current parameter selection!": "107",
        "SALT3: No data points with S/N > 5.0!": "434",
        "GHOST failed to intiate!": "222",
        "GHOST failed to find mass!": "333",
        "Unknown fitting error!": "404"
    }
    # Read current report
    data = np.genfromtxt(report_loc, skip_header=7, dtype=str, delimiter=';')
    report_tbl = Table(names=data[0, :], data=data[1:, :], dtype=[str, list, list, list, list, list, list, list])

    # Update report
    for k in fail_reasons.keys():
        report_tbl[error_codes[k]][np.where(report_tbl['source'] == source)[0][0]] = fail_reasons[k]

    # Print out report
    print("========================================================================================")
    for k, col in zip(error_codes.keys(), report_tbl[report_tbl['source'] == source].colnames[1:]):
        print(f"{error_codes[k]}, {k}: {list(report_tbl[report_tbl['source'] == source][col])[0]}")
    print("========================================================================================")


    # for k in error_codes.keys():
    #     print(f"{error_codes[k]} : {k}")
    # for col in report_tbl[report_tbl['source'] == source].colnames:
    #     print(f"{col}: {list(report_tbl[report_tbl['source'] == source][col])}")
    # # print()  # Prints report

    # Save new report
    with open(report_loc, 'w') as f:
        for k in error_codes.keys():
            print(f"# {error_codes[k]}: {k}", file=f)
        hdr = (str(list(error_codes.values()))[1:-1].replace("'", "").
               replace(" ", "").replace(",", ";"))
        print(f"source;{hdr}", file=f)
        for row in report_tbl:
            print(row[0], end='', file=f)
            for n in range(1, len(row)):
                print(f";{row[n]}", end='', file=f)
            print('\n', file=f)
    return error_codes, report_tbl
def fit_subprocess(dataset: str, path: str, algo: str, rewrite: bool = False):
    """
    :param dataset: Data set to pull light curves from
    :param path: Path to pull lightcurve
    :param algo: Algorithm to use for fitting
    :param rewrite: Even if data is already procced, it will act as if its the not
    :return: sneObj object
    """
    # Check if preveiously fit
    class_save_loc = f"classes/{dataset}/{dataset}_{path.split('/')[-1].split('.txt')[0]}_{algo}_class.txt"
    print(f"[+++] Checking for saved class at {class_save_loc}...")
    if not rewrite and os.path.exists(class_save_loc):
        print(f"[+++] Class found! Loading from {class_save_loc}...")
        sn = sneObj('class', algo, class_save_loc)
    else:
        print(f"[---] No class found! Setting up class light curve at {path}...")
        sn = sneObj(dataset, algo, path)

    # Get CMB Redshift, only if valid z, RA, and DEC
    if np.isnan(sn.z_cmb) and sn.coords[0] > -99999 and sn.coords[1] > -99999:
        print(f"[+++] Calculating CMB Redshift and adjusting for peculiar velocities...")
        try:
            sn.get_zcmb()
        except:
            sn.z_cmb = float(sn.z)

    # Fit SN class - get mu+color+stretch
    if len(sn.params) == 0 or rewrite:
        print(f"[+++] Fitting '{sn.objname}' with '{dataset}' data & the '{algo}' algorithm...")
        if algo.lower() == 'snpy':
            sn.snpy_fit()
        elif algo.lower() == 'salt':
            sn.salt_fit()
    else:
        print(f"[+++] '{sn.objname}' already fit! Loading...")

    # Update TNS Key
    update_TNS = True
    if update_TNS:
        if utils.check_tnskey(sn.coords[0], sn.coords[1]) is None:
            print(f"[+++] Updating TNS with '{sn.objname}' data...")
            utils.append_tnskey(sn.coords[0], sn.coords[1], sn.objname, sn.z, sn.discdate)

    # Get host mass
    if sn.params['mu']['value'] > 0:
        print(f'[+++] Finding host galaxy mass for {sn.objname} using GHOST...')
        sn.get_hostMass()
        if sn.params['hostMass']['value'] < 0:
            print(f"[!!!] Host mass is invalid ({sn.params['hostMass']['value']}+/-{sn.params['hostMass']['err']})")
    else:
        sn.params.update({'hostMass': {'value': np.nan, 'err': np.nan}})

    sn.save_class()
    return sn
def fit(data_loc: str, algo: str, rewrite: bool = False) -> list[sneObj]:
    """
    :param data_loc: Location of data; if single path -> indivisual mode, if directory -> batch mode
    :param algo: Algorithm to fit; either SNooPy or SALT3
    :param rewrite: Even if data is already procced, it will act as if its the not
    :return: sneObj object or list of sneObj objects
    """
    # Check arguments
    paths = glob.glob(data_loc)
    if len(paths) == 0:
        print(f"[!!!] No data found in {data_loc}")
        return []
    dataset = paths[0].split('/')[-2].lower()

    # Adjust for 'combined_91bg' dataset
    if dataset == 'combined_91bg':
        if 'CSP' in paths[0]: dataset = 'csp-91bg'
        elif 'ATLAS' in paths[0]: dataset = 'atlas-91bg'
        elif 'ZTF' in paths[0]: dataset = 'ztf-91bg'
    elif dataset == 'combined_norm':
        if 'CSP' in paths[0]: dataset = 'csp-norm'
        elif 'ATLAS' in paths[0]: dataset = 'atlas-norm'
        elif 'ZTF' in paths[0]: dataset = 'ztf-norm'

    # Verify proper dataset & proper algorithm
    valid_datasets = ['csp-91bg', 'atlas-91bg', 'ztf-91bg', 'combined_91bg',
                      'csp-norm', 'atlas-norm', 'ztf-norm', 'combined_norm']
    valid_algorithms = ['snpy', 'salt']
    if dataset not in valid_datasets: raise ValueError(f"[!!!] Dataset, '{dataset}', not recognized! {valid_datasets}")
    elif algo not in valid_algorithms: raise ValueError(f"[!!!] Algorithm, '{algo}', not recognized! {valid_algorithms}")

    # Select fitting mode
    ## Indivisual mode
    if len(paths) == 1:
        print(f"[+++] Fitting file '{data_loc}'...")  # Indivisual fit
        sn = fit_subprocess(dataset, paths[0], algo, rewrite)
        return [sn]
    ## Batch mode
    elif len(paths) > 1:
        print(f"[+++] Fitting data in '{data_loc}'...")  # Batch fit
        success_counter, sne, fail_reasons  = 0, [], []
        for i, path in enumerate(paths):
            print(f'[{i + 1} / {len(paths)}] ================================================================')
            sn = fit_subprocess(dataset, path, algo, rewrite)
            if (('mu' in sn.params) and ('hostMass' in sn.params) and
                (sn.params['mu']['value'] > 0) and (sn.params['hostMass']['value'] > 0) and (sn.discdate > 0)):
                print(f"\t********************************\n"
                      f"\t{sn.objname}: {sn.params['mu']['value']} +/- {sn.params['mu']['err']}\n"
                      f"\t{' '*len(sn.objname)}: {sn.params['hostMass']['value']} +/- {sn.params['hostMass']['err']}\n"
                      f"\t********************************")
            else:
                print(f"\t---------------[!]--------------")
                if ('hostMass' in sn.params) and (sn.params['hostMass']['value'] < 0):
                      print(f"\t{sn.objname}: Host Mass Faliure!")
                elif ('mu' in sn.params) and (sn.params['mu']['value'] < 0):
                      print(f"\t{sn.objname}: Fitting Faliure!")
                else:
                    print(f"\t{sn.objname}: *Unknown* Faliure!")
                print(f"\t---------------[!]--------------")
            sne.append(sn)

        # Check errors
        if '_' in dataset: s, t = dataset.split('_')
        else: s, t = dataset.split('-')
        good_sne = error_reporting(sne, f"{algo.upper()}_{t.upper()}_{s.upper()}", print_out=True)

        print(f"[+++++] Successfuly fit {len(good_sne)} / {len(sne)} SNe!")
        print("========================================================================================")

        return good_sne
    else:
        print('[!!!] Invalid file/data path!')
    return


if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
