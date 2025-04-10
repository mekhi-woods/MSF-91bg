import os
import json
import glob
import time
import utils  # Import of utils.py
import getVpec # Import of getVpec.py
import datetime
import numpy as np
import time as systime
import snpy as snoopyfit
import sncosmo as salt3fit
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.table import Table, conf
from astroquery.sdss import SDSS
from astropy.time import Time as astrotime
from astropy.coordinates import SkyCoord, Galactic
from astro_ghost.ghostHelperFunctions import getTransientHosts, getGHOST

class sneObj:
    # Construction Functions ----------------------------------------------------------------------------------------- #
    def __init__(self, source: str, algo: str, path: str = ''):
        # Empty Class
        if source == 'empty' or len(path) == 0:
            print('[+++] Creating empty SNe object...')
            pass

        # Loaded Class
        elif source == 'class' and os.path.exists(path):
            print('[+++] Creating SNe object from class file...')
            self.load_class(path)

        # Create new class
        else:
            # Open file
            with open(path, 'r') as f: hdr = f.readline().rstrip('\n').split(',')
            data = np.genfromtxt(path, delimiter=',', skip_header=1, dtype=str)
            tb = Table()
            for h in hdr:
                try: tb[h] = data[:, hdr.index(h)].astype(float)
                except ValueError: tb[h] = data[:, hdr.index(h)]
            # with open(path, 'r') as f: hdr = f.readline().rstrip('\n').split(',')
            # tb = Table(names=hdr,
            #            data=np.genfromtxt(path, delimiter=',', skip_header=1, dtype=str))

            # Set basic elements
            subtype, self.origin, name = path.split('_')
            self.algo = algo
            self.objname = name[:-4]
            self.originalname = self.origin+self.objname
            self.path = f"classes/{subtype.split('/')[-1]}_{self.origin}_{self.objname}_{algo}_class.txt"
            self.z = tb['z'][0]
            self.z_err = tb['z_err'][0]
            self.discdate = tb['discdate'][0]
            self.coords = [np.average(tb['ra']), np.average(tb['dec'])]

            # Intialize for later
            self.z_cmb = np.nan
            self.covariance = []
            self.params = {'mu':       {'value': np.nan, 'err': np.nan},
                           'hostMass': {'value': np.nan, 'err': np.nan}}

            # Set arrays
            self.mjd = np.array(tb['mjd'])
            self.mag = np.array(tb['mag'])
            self.dmag = np.array(tb['dmag'])
            self.flux = np.array(tb['flux'])
            self.dflux = np.array(tb['dflux'])
            self.filter = np.array(tb['filter'])
            self.zp = np.array(tb['zp'])

            # Save class file
            self.save_class()

        #     # Adjust ztf time scale
        #     if source == 'ztf-91bg' or source == 'combined_91bg' or source == 'combined_norm':
        #         var_tbl = var_tbl[(var_tbl['mjd'] > self.discdate - 30) & (var_tbl['mjd'] < self.discdate + 120)]

        return
    def save_class(self):
        print(f"[+++] Saving {self.objname} class to {self.path}...")
        with open(self.path, 'w') as f:
            f.write(f"{self.origin},{self.objname},{self.originalname},{self.coords[0]},{self.coords[1]},{self.z},{self.z_err},{self.z_cmb},{self.discdate},{self.algo}\n")
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for p in self.params:
                f.write(f"{p},{self.params[p]['value']},{self.params[p]['err']}\n")
            # f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            # for line in self.covariance:
            #     p_line = str(line[0])
            #     for i in range(1, len(line)):
            #         p_line += ',' + str(line[i])
            #     f.write(p_line + '\n')
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
            self.z_err = None if details[6] == "None" else float(details[6])
            self.z_cmb = float(details[7])
            self.discdate = details[8]
            self.algo = details[9][:-1]
            f.readline()  # Skip break line

            # Read params
            self.params = {}
            line = f.readline()
            while '+++' not in line:
                line = line.split(',')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            # # Read covariances
            # self.covariance = np.array([])
            # if 'salt' in path:
            #     line = f.readline()
            #     while '+++' not in line:
            #         line = line.split(',')
            #         line[-1] = line[-1][:-1]
            #         self.covariance = np.append(self.covariance, np.array(line).astype(float))
            #         line = f.readline()
            #     if len(self.covariance) > 0:
            #         self.covariance = self.covariance.reshape(4, 4)
            # else:
            #     f.readline()  # Skip break line

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

        # # Lines
        # axs.axvline(self.discdate, color='black', linestyle='--', label='Discovery Date')
        # try: axs.axvline(self.params['Tmax']['value'], color='maroon', linestyle='--', label='Peak Brightness')
        # except: pass

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

                # # Check number of filters
                # if len(n_s.data.keys()) < 3:
                #     print('[!!!] Too few filters to fit (<3)')
                #     self.params.update({'mu': {'value': -124.0, 'err': -124.0}})
                #     return

                # Fit
                n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
                n_s.save(model_path)

                # Save parameters
                for j in range(len(param_names)):
                    self.params.update({param_names[j]: {'value': n_s.parameters[snpy_param_names[j]],
                                                         'err': n_s.errors[snpy_param_names[j]]}})
                self.params.update({'chisquare': {'value': (n_s.model.chisquare/n_s.model.dof),
                                                  'err': 0.00}})
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

            # if len(np.unique(data['band'])) < 3:
            #     print(f"[!!!] Not enough filters to fit! ({len(np.unique(data['band']))})")
            #     self.params.update({'mu': {'value': -124.0, 'err': -124.0}})
            #     return

            # Get JD discdate
            discdate = astrotime(self.discdate, format='iso', scale='utc').mjd

            # Create and fit data to model
            model = salt3fit.Model(source='salt3')
            model.set(z=self.z)  # set the model's redshift.
            # model.set(z=self.z, t0=discdate)  # set the model's redshift.
            result, fitted_model = salt3fit.fit_lc(data, model,
                                                   ['t0', 'x0', 'x1', 'c'],
                                                   bounds={'x1': (-5, 5)},
                                                   minsnr=5)

            # Save Parameters
            param_names = ['t0', 'x0', 'x1', 'c']
            for i in range(len(param_names)):
                self.params.update({param_names[i]: {'value': result.parameters[i+1],
                                                     'err': result.errors[param_names[i]]}})
            self.params.update({'chisquare': {'value': result.chisq/result.ndof,
                                              'err': 0.00}})

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
def fit(data_loc: str, algo: str, rewrite: bool = False):
    # Secure data location to list
    if type(data_loc) == str: data_loc = [data_loc]

    good_sne = []
    hostMass_fail = []
    fitting_fail = []
    for i, path in enumerate(data_loc):
        print(f'[{i + 1} / {len(data_loc)}] ================================================================')
        print(f"[---] Fitting data from '{path}'...")
        # Get info
        subtype, origin, name = path.split('_')
        name = name[:-4]
        subtype = subtype.split('/')[-1]

        # Check if preveiously fit
        if not rewrite and os.path.exists(f"classes/{subtype}_{origin}_{name}_{algo}_class.txt"):
            print(f"[---] Class for {name} already exsists, pulling...")
            sn = sneObj('class', algo, f"classes/{subtype}_{origin}_{name}_{algo}_class.txt")
            if ((sn.params['mu']['value'] < 0) or (sn.params['hostMass']['value'] < 0) or
                (np.isnan(sn.params['mu']['value'])) or (np.isnan(sn.params['hostMass']['value']))):
                print('[!!!] Retrieve class shows previously failed fit...\n')
                continue
        # Start fitting
        else:
            # Make class
            sn = sneObj(subtype, algo, path)

            # Get CMB Redshift, only if valid z, RA, and DEC
            if np.isnan(sn.z_cmb) and sn.coords[0] > -99999 and sn.coords[1] > -99999:
                print(f"[+++] Calculating CMB Redshift and adjusting for peculiar velocities...")
                try:
                    sn.get_zcmb()
                except:
                    sn.z_cmb = float(sn.z)

            # Fit SN class - get mu+color+stretch
            if np.isnan(sn.params['mu']['value']) or rewrite:
                print(f"[+++] Fitting '{sn.objname}' with '{sn.origin}' data & the '{algo}' algorithm...")
                try:
                    if algo.lower() == 'snpy':
                        sn.snpy_fit()
                    elif algo.lower() == 'salt':
                        sn.salt_fit()
                except Exception as e:
                    sn.params['mu']['value'], sn.params['mu']['err'] = -124, -124
                if sn.params['mu']['value'] < 0:
                    print(f"[!!!] Fitting failiure! Skipping {sn.objname}...\n")
                    fitting_fail.append(sn)
                    continue
            else:
                print(f"[+++] '{sn.objname}' already fit! Loading...")

            # Update TNS Key
            update_TNS = True
            if update_TNS:
                if utils.check_tnskey(sn.coords[0], sn.coords[1]) is None:
                    print(f"[+++] Updating TNS with '{sn.objname}' data...")
                    utils.append_tnskey(sn.coords[0], sn.coords[1], sn.objname, sn.z, sn.discdate)

            # Get host mass
            try:
                print(f'[+++] Finding host galaxy mass for {sn.objname} using GHOST...')
                sn.get_hostMass()
                if sn.params['hostMass']['value'] < 0:
                    print(f"[!!!] Host mass is invalid ({sn.params['hostMass']['value']}+/-{sn.params['hostMass']['err']})")
            except Exception as e:
                sn.params['hostMass']['value'], sn.params['hostMass']['err'] = -303, -303
            if sn.params['hostMass']['value'] < 0:
                print(f"[!!!] Host mass failiure! Skipping {sn.objname}...\n")
                hostMass_fail.append(sn)
                continue

            # Save class
            sn.save_class()

            #


        # Add to SNe list
        good_sne.append(sn)
        print(f"\t********************************\n"
              f"\t{sn.objname}: {sn.params['mu']['value']} +/- {sn.params['mu']['err']}\n"
              f"\t{' '*len(sn.objname)}: {sn.params['hostMass']['value']} +/- {sn.params['hostMass']['err']}\n"
              f"\t{' '*len(sn.objname)}: {sn.params['chisquare']['value']}\n"
              f"\t********************************\n")
    return good_sne

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
