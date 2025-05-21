import os
import glob
import utils  # Import of utils.py
import shutil
import datetime
import numpy as np
import time as systime
from astropy.table import Table, vstack
from astropy.stats import sigma_clip, sigma_clipped_stats
from scipy.stats import norm  # For inverse normal CDF (percentile function)

CURRENTDATE = datetime.datetime.now()


def make_objTbl(source: str, path: str) -> Table:
    # Load data and make intial table
    if source == 'atlas':
        # Load data
        with open(path, 'r') as f:
            hdr = f.readline()[1:-1].split(',')
        data = np.genfromtxt(path, dtype='str', delimiter=',', skip_header=1)

        # Check if empty
        if len(data) == 0:
            print("[!!!] Empty file...")
            return None

        # Make table
        var_table = Table()
        for h in hdr:
            try:
                var_table[h] = data[:, hdr.index(h)].astype(float)
            except ValueError:
                var_table[h] = data[:, hdr.index(h)]

        # Add parity with CSP & ZTF
        var_table.remove_columns(['err', 'x', 'y', 'maj', 'min', 'phi', 'apfit', 'Sky', 'Obs'])
        for h_old, h_new in zip(['JD', 'm', 'dm', 'uJy', 'duJy', 'F', 'RA', 'Dec'],
                                ['mjd', 'mag', 'dmag', 'flux', 'dflux', 'filter', 'ra', 'dec']):
            var_table[h_old].name = h_new

        # Recreate Zero-Points using ZP = m + 2.5log_10(flux)
        var_table['zp'] = np.full(var_table['mag'].shape, 23.9)

        # Filter on 5-sigma magnitudes (mag5sig)
        var_table = var_table[var_table['mag'] > var_table['mag5sig']]

        # Filter out chi/N > 2
        var_table = var_table[var_table['chi/N'] < 2.0]

        # Cut errors > 1mag
        var_table = var_table[var_table['dmag'] < 1.0]
    elif source == 'csp':
        var_table = Table(names=['filter', 'zp', 'z', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux'],
                          dtype=[str, float, float, float, float, float, float, float, float, float])
        with open(path, 'r') as f:
            csp_objname, csp_z, csp_ra, csp_dec = f.readline()[2:-1].split(' ')
            for l in f.readlines():
                l = l.split(' ')
                # Filter line
                if len(l) == 2:
                    csp_filter = str(l[1][:-1])
                    csp_zp = float(utils.get_constants()['csp_zpts_' + csp_filter])
                else:
                    csp_mjd, csp_mag, csp_dmag = float(l[-3]) + 53000, float(l[-2]), float(l[-1])
                    csp_flux = 10 ** ((csp_mag - csp_zp) / -2.5)
                    csp_dflux = np.abs(csp_flux) * np.log(10) * ((1 / 2.5) * csp_dmag)

                    var_table.add_row([csp_filter, csp_zp, csp_z, csp_ra, csp_dec, csp_mjd,
                                       csp_mag, csp_dmag, csp_flux, csp_dflux])
        # self.z = var_table['z'][0]  # Special case for CSP since CSP-I didn't report all SNe redshifts to TNS
    elif source == 'ztf':
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

        # Check if empty
        if len(data) == 0:
            print("[!!!] Empty file...")
            return None

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
        var_table = var_table[var_table['forcediffimflux'].astype(float) > 0] # Remove negative values
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
                                  'programid', 'rfid', 'forcediffimsnr',
                                  'forcediffimfluxap', 'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr',
                                  'dnearestrefsrc', 'nearestrefmag', 'nearestrefmagunc', 'nearestrefchi',
                                  'nearestrefsharp', 'refjdstart', 'refjdend', 'procstatu'])
        for h_old, h_new in zip(['zpdiff', 'jd', 'forcediffimflux', 'forcediffimfluxunc'],
                                ['zp', 'mjd', 'flux', 'dflux']):
            var_table[h_old].name = h_new

        # Filter on by magnitude limit (diffmaglim)
        var_table = var_table[var_table['mag'] < var_table['diffmaglim']]

        # Filter out by reduced chi2 > 3.0 -- Might remove too many filters
        # var_table = var_table[var_table['forcediffimchisq'] < 3.0]

        # Cut errors > 1mag
        var_table = var_table[var_table['dmag'] < 1.0]

    # Check for filter >= 2
    if len(np.unique(var_table['filter'])) < 2:
        print("[!!!] <2 Filters...")
        return None

    # Check if too few points in filters
    valid_filters = 0
    for f in np.unique(var_table['filter']):
        if len(var_table[var_table['filter'] == f]) > 10:
            valid_filters += 1
    if valid_filters < 2:
        print("[!!!] Too few filters with +10 data points...")
        return None

    # Check if empty, again
    if len(var_table) == 0:
        print("[!!!] Empty file...")
        return None

    return var_table
def check_combined_stats(subtype: str = '91bg'):
    csp, atlas, ztf, csp_atlas, csp_ztf, atlas_ztf = 0, 0, 0, 0, 0, 0
    for file in glob.glob(f'data/combined_{subtype}/*.txt'):
        if 'CSP' in file.split('/')[-1]: csp += 1
        elif 'ATLAS' in file.split('/')[-1]: atlas += 1
        elif 'ZTF' in file.split('/')[-1]: ztf += 1
        else:
            with open(file, 'r') as f:
                hdr = f.readline().rstrip('\n').replace("'", "").replace(",", "").split(' ')
                combo = f"{hdr[-1]}+{hdr[-2]}"
                if combo == 'CSP+ATLAS' or combo == 'ATLAS+CSP': csp_atlas += 1
                elif combo == 'CSP+ZTF' or combo == 'ZTF+CSP': csp_ztf += 1
                elif combo == 'ATLAS+ZTF' or combo == 'ZTF+ATLAS': atlas_ztf += 1
    print(f"CSP: {csp} | ATLAS: {atlas} | ZTF: {ztf} | CSP+ATLAS: {csp_atlas} | CSP+ZTF: {csp_ztf} | ATLAS+ZTF: {atlas_ztf}")
    return
def combine_like_data(target_file: str, combined_loc: str, subtype: str, clear_old: bool = False):
    # Remove old combined data
    if clear_old:
        if os.path.exists(combined_loc):
            shutil.rmtree(combined_loc)
        os.mkdir(combined_loc)

    # Combine all files into one list
    all_files = (glob.glob(f'data/CSP-{subtype}/*.txt') +
                 glob.glob(f'data/ATLAS-{subtype}/*.txt') +
                 glob.glob(f'data/ZTF-{subtype}/*.txt'))

    # With that list, get the names and origins of data
    names_orgins = Table(names=['name', 'origin'], dtype=[str, str])
    for file in all_files:
        origin = file.split('/')[-2].split('-')[0]
        name = file.split('/')[-1].split('.')[0].replace(origin, '')
        names_orgins.add_row([name, origin])

    # Open target file
    target_file = np.genfromtxt(target_file, delimiter=',', dtype=str, skip_header=1)
    tb_targets = Table(names=target_file[0], data=target_file[1:])

    # Go through list of targets, check if they exsist in data, and take the name of the surveys they're in
    tb_combined = Table(names=list(target_file[0])+['Sources'], dtype=[str]*(len(list(target_file[0]))+1))
    for i, row in enumerate(tb_targets):
        all_sources = ''
        tb_overlaping_sources = names_orgins[names_orgins['name'] == row['Name'][3:]]
        for n_row in tb_overlaping_sources:
            all_sources = all_sources+n_row['origin']+'+'
        if len(tb_overlaping_sources) == 0:
            tb_combined.add_row(list(row) + ['None'])
        else:
            tb_combined.add_row(list(row) + [all_sources[:-1]])

    # Go through overlapping targets and combine data files
    for i, row in enumerate(tb_combined):
        name = row['Name'][3:]

        # Combine data from >1 source
        if row['Sources'] != 'None' and '+' in row['Sources']:
            n_tb = Table(names=['filter', 'zp', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux', 'source'],
                         dtype=[str, float, float, float, float, float, float, float, float, str])
            for source in row['Sources'].split('+'):
                var_table = make_objTbl(source.lower(), f"data/{source}-{subtype}/{source}{name}.txt")
                if var_table is not None:
                    for n_row in var_table:
                        n_tb.add_row([n_row['filter'], n_row['zp'], n_row['ra'], n_row['dec'], n_row['mjd'],
                                      n_row['mag'], n_row['dmag'], n_row['flux'], n_row['dflux'], source])
            with open(f"{combined_loc}{name}.txt", 'w') as f:
                file_hdr = str(n_tb.colnames).strip('[]').replace("'", "")
                print(f"# Combination of '{name}' data from {str(row['Sources'].split('+')).strip('[]')}", file=f)
                print(file_hdr, file=f)
                for r in n_tb:
                    print(str(list(r)).strip('[]').replace("'", ""), file=f)

        # Take copy for later
        elif row['Sources'] != 'None':
          shutil.copyfile(f"data/{row['Sources']}-{subtype}/{row['Sources']}{name}.txt",
                          f"{combined_loc}{row['Sources']}{name}.txt")

    return
def new_combine_like_data(target_file: str, combined_loc: str, subtype: str, clear_old: bool = False):
    # Remove old combined data
    if clear_old:
        if os.path.exists(combined_loc):
            shutil.rmtree(combined_loc)
        os.mkdir(combined_loc)

    # Combine all files into one list
    csp_files = glob.glob(f'data/CSP-{subtype}/*')
    atlas_files = glob.glob(f'data/ATLAS-{subtype}/*')
    ztf_files = glob.glob(f'data/ZTF-{subtype}/*')

    # Get lists of names of SNe in file system
    csp_f_names, atlas_f_names, ztf_f_names = [], [], []
    for names, files in zip([csp_f_names, atlas_f_names, ztf_f_names], [csp_files, atlas_files, ztf_files]):
        for f in files:
            n = (str(f.split('/')[-1][:-4]).
                 replace("CSP", "").replace("ATLAS", "").replace("ZTF", ""))
            names.append(n)
    all_names = csp_f_names+atlas_f_names+ztf_f_names

    # Open target list
    tarlist_data = np.genfromtxt(target_file, delimiter=',', dtype=str)
    tarlist = Table(names=tarlist_data[0, :], data=tarlist_data[1:, :])

    # Sort through target list and find matches
    for row in tarlist:
        n = row['Name'][3:]

        if (n in csp_f_names):
            # Only CSP data
            shutil.copy(f'data/CSP-{subtype}/CSP{n}.txt', f"{combined_loc}CSP{n}.txt")
        if (n in atlas_f_names) and (n in ztf_f_names):
            # Both ATLAS & ZTF data

            # Load ATLAS data from file
            atlas_data = np.genfromtxt(f"data/ATLAS-{subtype}/ATLAS{n}.txt", delimiter=',', dtype=str)
            atlas_hdr, atlas_data = list(atlas_data[0, :]), atlas_data[1:, :]
            atlas_ra = atlas_data[:, atlas_hdr.index('RA')]
            atlas_dec = atlas_data[:, atlas_hdr.index('Dec')]
            atlas_filter = atlas_data[:, atlas_hdr.index('F')]
            atlas_mjd = atlas_data[:, atlas_hdr.index('MJD')]
            atlas_mag = atlas_data[:, atlas_hdr.index('m')]
            atlas_dmag = atlas_data[:, atlas_hdr.index('dm')]
            atlas_flux = atlas_data[:, atlas_hdr.index('uJy')]
            atlas_dflux = atlas_data[:, atlas_hdr.index('duJy')]

            # Load ZTF data from file
            ztf_data = np.genfromtxt(f"data/ZTF-{subtype}/ZTF{n}.csv", delimiter=',', dtype=str)
            ztf_hdr, ztf_data = list(ztf_data[0, :]), ztf_data[1:, :]
            ztf_ra = ztf_data[:, ztf_hdr.index('ra')]
            ztf_dec = ztf_data[:, ztf_hdr.index('dec')]
            ztf_filter = ztf_data[:, ztf_hdr.index('filtercode')]
            ztf_mjd = ztf_data[:, ztf_hdr.index('mjd')]
            ztf_mag = ztf_data[:, ztf_hdr.index('mag')]
            ztf_dmag = ztf_data[:, ztf_hdr.index('magerr')]
            ztf_zp = ztf_data[:, ztf_hdr.index('magzp')]

            # Calculate ATLAS zeropoints using: ZP = m + 2.5log_10(flux)
            # atlas_zp = -2.5 * np.log10(atlas_flux.astype(float)) + 23.9
            atlas_zp = atlas_mag.astype(float) + (2.5 * np.log10(atlas_flux.astype(float)))

            # Calculate ZTF flux values
            ztf_flux = 10 ** ((ztf_mag.astype(float) - ztf_zp.astype(float)) / (-2.5))
            ztf_dflux = np.abs(ztf_flux) * (1 / 2.5) * np.log(10) * ztf_dmag.astype(float)

            # Make table of combined data
            tb = Table(names=['ra', 'dec', 'filter', 'mjd', 'mag', 'dmag', 'flux', 'dflux', 'zp', 'source'],
                       data=[np.hstack((atlas_ra, ztf_ra)), np.hstack((atlas_dec, ztf_dec)),
                             np.hstack((atlas_filter, ztf_filter)), np.hstack((atlas_mjd, ztf_mjd)),
                             np.hstack((atlas_mag, ztf_mag)), np.hstack((atlas_dmag, ztf_dmag)),
                             np.hstack((atlas_flux, ztf_flux)), np.hstack((atlas_dflux, ztf_dflux)),
                             np.hstack((atlas_zp, ztf_zp)),
                             np.full(np.hstack((atlas_ra, ztf_ra)).shape, 'ATLAS-ZTF')])

            # Write data to new location
            with open(f"{combined_loc}{n}.txt", 'w') as f:
                print(f"# Combined data from 'ATLAS-ZTF' for 'SN {n}'", file=f)
                print(str(tb.colnames)[1:-1].replace("'", ""), file=f)  # Column Names
                for n_row in tb:
                    line = ''
                    for n_col in n_row:
                        line += f"{n_col}, "
                    print(line[:-2], file=f)
        if (n in atlas_f_names) and (n not in ztf_f_names):
            # Only ATLAS data
            shutil.copy(f'data/ATLAS-{subtype}/ATLAS{n}.txt', f"{combined_loc}ATLAS{n}.txt")
        if (n in ztf_f_names) and (n not in atlas_f_names):
            # Only ZTF data
            shutil.copy(f'data/ZTF-{subtype}/ZTF{n}.csv', f"{combined_loc}ZTF{n}.txt")
    return
def old_old_combine_snpy_salt(snpy_path: str, salt_path: str, save_loc: str = ''):
    tb_snpy = utils.default_open(snpy_path, 'True')
    tb_salt = utils.default_open(salt_path, 'True')

    # Check if one table is empty
    if len(tb_snpy) == 0:
        tb_combined = tb_salt.copy()
    elif len(tb_salt) == 0:
        tb_combined = tb_snpy.copy()
    else:
        tb_combined = tb_salt.copy()
        for i, n in enumerate(tb_snpy['objname']):
            if n in list(tb_salt['objname']):  # Identify duplicate
                # Check mu_err; if SNPY better -> del SALT row, add SNPY row, else pass
                if tb_snpy['mu_err'][tb_snpy['objname'] == n] < tb_salt['mu_err'][tb_salt['objname'] == n]:
                    tb_combined = tb_combined[tb_combined['objname'] != n]
                    tb_combined.add_row(tb_snpy[i])
            else:  # If no duplicate, add SNooPy row to combined table
                tb_combined.add_row(tb_snpy[i])

    # Save combined table
    if len(save_loc) > 0:
        print(f"[~~~] Saving combined file to '{save_loc}'...")
        with open(save_loc, 'w') as f:
            print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_combined)}", file=f)
            print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
            print(f"objname, origin, subtype, algo, ra, dec, z, z_cmb, MJDs, MJDe, mu, mu_err, hostMass, hostMass_err, "
                  f"Tmax, Tmax_err, stretch, stretch_err, color, color_err, amplitude, amplitude_err, "
                  f"peak_mag, peak_mag_err", file=f)
            for row in tb_combined:
                print(f"{row['objname']}, "
                      f"{row['origin']}, "
                      f"{row['subtype']}, "
                      f"{row['algo']}, "
                      f"{row['ra']}, "
                      f"{row['dec']}, "
                      f"{row['z']}, {row['z_cmb']}, "
                      f"{row['MJDs']}, {row['MJDe']}, "
                      f"{row['mu']}, {row['mu_err']}, "
                      f"{row['hostMass']}, {row['hostMass_err']}, "
                      f"{row['Tmax']}, {row['Tmax_err']}, "
                      f"{row['stretch']}, {row['stretch_err']}, "
                      f"{row['color']}, {row['color_err']}, "
                      f"{row['amplitude']}, {row['amplitude_err']}, "
                      f"{row['peak_mag']}, {row['peak_mag_err']}", file=f)
    return
def old_combine_snpy_salt(snpy_path: str, salt_path: str, save_loc: str = ''):
    tb_snpy = utils.default_open(snpy_path, 'True')
    tb_salt = utils.default_open(salt_path, 'True')

    # Check if one table is empty
    if len(tb_snpy) == 0:
        tb_combined = tb_salt.copy()
    elif len(tb_salt) == 0:
        tb_combined = tb_snpy.copy()
    else:
        tb_combined = tb_salt.copy()
        for i, n in enumerate(tb_snpy['objname']):
            if n in list(tb_salt['objname']):  # Identify duplicate
                # Check mu_err; if SNPY better -> del SALT row, add SNPY row, else pass
                # Combine mu and mu_err from SNPY & SALT
                if tb_snpy['mu_err'][tb_snpy['objname'] == n] < tb_salt['mu_err'][tb_salt['objname'] == n]:
                    mu_snpy, mu_err_snpy = float(tb_snpy['mu'][tb_snpy['objname'] == n]), float(tb_snpy['mu_err'][tb_snpy['objname'] == n])
                    mu_salt, mu_err_salt = float(tb_salt['mu'][tb_salt['objname'] == n]), float(tb_salt['mu_err'][tb_salt['objname'] == n])

                    mu_new = np.average([mu_snpy, mu_salt], weights=[1 / mu_err_snpy ** 2, 1 / mu_err_salt ** 2])
                    mu_err_new = np.average([mu_err_snpy, mu_err_salt])

                    match_index = np.where(tb_combined['objname'] == n)[0][0]
                    tb_combined[match_index]['mu'] = mu_new
                    tb_combined[match_index]['mu_err'] = mu_err_new
                    tb_combined[match_index]['mu_err'] = mu_err_new
            else:  # If no duplicate, add SNooPy row to combined table
                tb_combined.add_row(tb_snpy[i])

    # Save combined table
    if len(save_loc) > 0:
        print(f"[~~~] Saving combined file to '{save_loc}'...")
        with open(save_loc, 'w') as f:
            print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_combined)}", file=f)
            print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
            print(f"objname, origin, subtype, algo, ra, dec, z, z_cmb, MJDs, MJDe, mu, mu_err, hostMass, hostMass_err, "
                  f"Tmax, Tmax_err, stretch, stretch_err, color, color_err, amplitude, amplitude_err, "
                  f"peak_mag, peak_mag_err", file=f)
            for row in tb_combined:
                print(f"{row['objname']}, "
                      f"{row['origin']}, "
                      f"{row['subtype']}, "
                      f"{row['algo']}, "
                      f"{row['ra']}, "
                      f"{row['dec']}, "
                      f"{row['z']}, {row['z_cmb']}, "
                      f"{row['MJDs']}, {row['MJDe']}, "
                      f"{row['mu']}, {row['mu_err']}, "
                      f"{row['hostMass']}, {row['hostMass_err']}, "
                      f"{row['Tmax']}, {row['Tmax_err']}, "
                      f"{row['stretch']}, {row['stretch_err']}, "
                      f"{row['color']}, {row['color_err']}, "
                      f"{row['amplitude']}, {row['amplitude_err']}, "
                      f"{row['peak_mag']}, {row['peak_mag_err']}", file=f)
    return
def combine_snpy_salt(snpy_path: str, salt_path: str, save_loc: str = ''):
    tb_snpy = utils.default_open(snpy_path, True)
    tb_salt = utils.default_open(salt_path, True)
    for name in np.unique(list(tb_snpy["objname"])+list(tb_salt["objname"])):
        if (name in list(tb_snpy["objname"])) & (name in list(tb_salt["objname"])):
            # Prefer SNooPy MU's
            tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])

    if len(tb_snpy) > len(tb_salt):
        tb_combined = tb_snpy.copy()
        for row in tb_salt:
            tb_combined.add_row(row)
    else:
        tb_combined = tb_salt.copy()
        for row in tb_snpy:
            tb_combined.add_row(row)

    with open(save_loc, 'w') as f:
        for comment in [f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_combined)}",
                        f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical."]:
            print(comment, file=f)
        hdr = list(tb_combined.columns)
        print(str(str(hdr).replace("'", "")[1:-1]), file=f)
        for row in tb_combined:
            line = ""
            for c in hdr:
                line += f"{row[c]}, "
            print(line[:-2], file=f)
    print(f"[+++] Saved data to {save_loc}...")
    return
def make_param_file(sne: list, save_loc: str):
    """
    :param sne: list of sneObjs to add to parameter file.
    :param save_loc: location to save parameter file.
    :return: None
    """
    print(f"[+++] Saving list of sneObjs to {save_loc}...")
    with open(save_loc, 'w') as f:
        print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(sne)}", file=f)
        print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
        print("objname, "
              "origin, "
              "subtype, "
              "algo, "
              "ra, "
              "dec, "
              "z, "
              "z_cmb, "
              "MJDs, "
              "MJDe, "
              "mu, mu_err, "
              "hostMass, hostMass_err, "
              "Tmax, Tmax_err, "
              "stretch, stretch_err, "
              "color, color_err, "
              "amplitude, amplitude_err, "
              "peak_mag, peak_mag_err", file=f)
        for sn in sne:
            if 'EBVhost' in list(sn.params.keys()) and 'st' in list(sn.params.keys()):
                param_dict = {'t_max': 'Tmax', 'stretch': 'st', 'color': 'EBVhost'}
            elif 'c' in list(sn.params.keys()) and 'x1' in list(sn.params.keys()):
                param_dict = {'t_max': 't0', 'stretch': 'x1', 'color': 'c'}

            if 'x0' in list(sn.params.keys()):
                amp, amp_err = sn.params['x0']['value'], sn.params['x0']['err']
            else:
                amp, amp_err = np.nan, np.nan
            # if sn.path.split('_')[-2].upper() == "SNPY":
            #     param_dict = {'t_max': 'Tmax', 'stretch': 'st', 'color': 'EBVhost'}
            # elif sn.path.split('_')[-2].upper() == "SALT":
            #     param_dict = {'t_max': 't0', 'stretch': 'x1', 'color': 'c'}
            # if sn.path.split('_')[-2].upper() == "SALT":
            #     amp, amp_err = sn.params['x0']['value'], sn.params['x0']['err']
            # else:
            # amp, amp_err = np.nan, np.nan
            print(f"{sn.objname}, "
                  f"{sn.path.split('/')[-2].replace('_', '-').split('-')[0].upper()}, "
                  f"{sn.path.split('/')[-2].replace('_', '-').split('-')[-1].upper()}, "
                  f"{sn.path.split('_')[-2].upper()}, "
                  f"{sn.coords[0]}, "
                  f"{sn.coords[1]}, "
                  f"{sn.z}, "
                  f"{sn.z_cmb}, "
                  f"{min(sn.mjd)}, "
                  f"{max(sn.mjd)}, "
                  f"{sn.params['mu']['value']}, {sn.params['mu']['err']}, "
                  f"{sn.params['hostMass']['value']}, {sn.params['hostMass']['err']}, "
                  f"{sn.params[param_dict['t_max']]['value']}, {sn.params[param_dict['t_max']]['err']}, "
                  f"{sn.params[param_dict['stretch']]['value']}, {sn.params[param_dict['stretch']]['err']}, "
                  f"{sn.params[param_dict['color']]['value']}, {sn.params[param_dict['color']]['err']}, "
                  f"{amp}, {amp_err}, "
                  f"{np.average(sn.mag[(sn.mjd.astype(float) > sn.params[param_dict['t_max']]['value'] - 5) & (sn.mjd.astype(float) < sn.params[param_dict['t_max']]['value'] + 5)].astype(float))}, "
                  f"{np.average(sn.dmag[(sn.mjd.astype(float) > sn.params[param_dict['t_max']]['value'] - 5) & (sn.mjd.astype(float) < sn.params[param_dict['t_max']]['value'] + 5)].astype(float))}",
                  file=f)
    return
def old_selection_criteria(path: str, subtype: str = '91bg', save_loc: str = ''):
    print(f"[+++] Selection criteria for '{path}'...")
    # 999 = no upper limit, -999 = no lower limit
    if subtype == '91bg':
        # Original cuts
        criteria = {
            'z_cmb':       [0.015, 999],
            'EBVhost':     [-0.3, 0.3],
            'EBVhost_err': [-999, 0.1],
            'st':          [-999, 1.0],
            'st_err':      [-999, 0.1],
            'c':           [-0.6, 0.6],
            'c_err':       [-999, 0.1],
            'x1':          [-3.2, 3.2],
            'x1_err':      [-999, 1.0],
            'Tmax_err':    [-999, 1.0],
            't0_err':      [-999, 1.0],
            'mu_err':      [-999, 0.2]
        }
    elif subtype == 'norm':
        criteria = {
            'z_cmb':       [0.015, 999],
            'EBVhost':     [-0.2, 0.3],
            'EBVhost_err': [-999, 0.1],
            'st':          [-999, 1.0],
            'st_err':      [-999, 0.1],
            'c':           [-0.6, 0.6],
            'c_err':       [-999, 0.1],
            'x1':          [-3.2, 3.2],
            'x1_err':      [-999, 1.0],
            'Tmax_err':    [-999, 1.0],
            't0_err':      [-999, 1.0],
            'mu_err':      [-999, 0.2]
        }

    # Load data
    tb = utils.default_open(path, True)

    # Basic cuts
    prev_tb = tb.copy()
    for p in ['z_cmb', 'Tmax_err', 'mu_err']:
        print(f"{p}: {len(tb)} -> ", end='')
        tb = tb[(tb[p] > criteria[p][0]) & (tb[p] <= criteria[p][1])]
        print(f"{len(tb)} ", end='')
        print(f"(σ={round(np.std(tb['mu'].astype(float) - utils.current_cosmo().distmod(tb['z_cmb'].astype(float)).value), 4)}) ", end='')
        print(list(set(list(prev_tb['objname'])).symmetric_difference(set(list(tb['objname'])))))
        # dif_tb = list(set(list(prev_tb['objname'])).symmetric_difference(set(list(tb['objname']))))
        # print(f"(x{len(dif_tb)})", dif_tb)
        prev_tb = tb.copy()

    # Seperate SNPY & SALT
    tb_snpy = tb[tb['algo'] == 'SNPY']
    tb_salt = tb[tb['algo'] == 'SALT']

    # Cut on SNPY params
    prev_tb = vstack([tb_snpy, tb_salt])
    for f_name, d_name in zip(['EBVhost', 'st'], ['color', 'stretch']):
        for n_end in ['', '_err']:
            print(f"{f_name+n_end}: {len(tb_snpy)+len(tb_salt)} -> ", end='')
            tb_snpy = tb_snpy[(tb_snpy[d_name+n_end] >= criteria[f_name+n_end][0]) &
                              (tb_snpy[d_name+n_end] <= criteria[f_name+n_end][1])]
            print(f"{len(tb_snpy)+len(tb_salt)} ", end='')

            tb_temp = vstack([tb_snpy, tb_salt])
            print(f"(σ={round(np.std(tb_temp['mu'].astype(float) - utils.current_cosmo().distmod(tb_temp['z_cmb'].astype(float)).value))}) ", end='')
            print(list(set(list(prev_tb['objname'])).symmetric_difference(set(list(tb_temp['objname'])))))
            # dif_tb = list(set(list(prev_tb['objname'])).symmetric_difference(set(list(tb_temp['objname']))))
            # print(f"(x{len(dif_tb)})", dif_tb)
            prev_tb = tb_temp.copy()

    # Cut on SALT params
    prev_tb = vstack([tb_snpy, tb_salt])
    for f_name, d_name in zip(['c', 'x1'], ['color', 'stretch']):
        for n_end in ['', '_err']:
            print(f"{f_name+n_end}: {len(tb_snpy)+len(tb_salt)} -> ", end='')
            tb_salt = tb_salt[(tb_salt[d_name+n_end] >= criteria[f_name+n_end][0]) &
                              (tb_salt[d_name+n_end] <= criteria[f_name+n_end][1])]
            print(f"{len(tb_snpy)+len(tb_salt)} ", end='')

            tb_temp = vstack([tb_snpy, tb_salt])
            print(f"(σ={round(np.std(tb_temp['mu'].astype(float) - utils.current_cosmo().distmod(tb_temp['z_cmb'].astype(float)).value), 4)}) ", end='')
            print(list(set(list(prev_tb['objname'])).symmetric_difference(set(list(tb_temp['objname'])))))
            # dif_tb = list(set(list(prev_tb['objname'])).symmetric_difference(set(list(tb_temp['objname']))))
            # print(f"(x{len(dif_tb)})", dif_tb)
            prev_tb = tb_temp.copy()

    # Recombine SNPY & SALT
    tb_combined = vstack([tb_snpy, tb_salt])

    # Chauvenet’s Criterion
    prev_tb = tb_combined.copy()
    print(f"Chauvenet’s Criterion: {len(tb_combined)} -> ", end='')
    resid = tb_combined['mu'].astype(float) - utils.current_cosmo().distmod(tb_combined['z_cmb'].astype(float)).value
    limit = norm.ppf(1 - (1.0 / (4 * len(resid))))
    tb_combined = tb_combined[np.abs((resid - np.average(resid, weights=(1/(tb_combined['mu_err'].astype(float)**2)))) / np.std(resid)) < limit]
    limit_num, limit_num_prev = len(tb_combined), len(prev_tb)
    print(f"{len(tb_combined)} "
          f"(σ={round(np.std(tb_combined['mu'].astype(float) - utils.current_cosmo().distmod(tb_combined['z_cmb'].astype(float)).value), 4)}) ", end='')
    print(list(set(list(prev_tb['objname'])).symmetric_difference(set(list(tb_combined['objname'])))))
    # dif_tb = list(set(list(prev_tb['objname'])).symmetric_difference(set(list(tb_combined['objname']))))
    # print(f"(x{len(dif_tb)})", dif_tb)

    # Adding intrinsic dispersion (0.1 mag) in quadrature for mass (taylor+11) & mu
    tb_combined['mu_err'] = np.sqrt(tb_combined['mu_err'] ** 2.0 + 0.1 ** 2.0)
    tb_combined['hostMass_err'] = np.sqrt(tb_combined['hostMass_err'] ** 2.0 + 0.1 ** 2.0)

    # Remaining SNe
    print(f"[~~~] Selected SNe: {list(tb_combined['objname'])}")

    # Save new parameter file
    if len(save_loc) != 0:
        print(f"[~~~] Saving new parameter file to '{save_loc}'...")
        with open(save_loc, 'w') as f:
            print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_combined)}", file=f)
            print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
            print(f"objname, origin, subtype, algo, ra, dec, z, z_cmb, MJDs, MJDe, mu, mu_err, hostMass, hostMass_err, "
                  f"Tmax, Tmax_err, stretch, stretch_err, color, color_err, amplitude, amplitude_err, "
                  f"peak_mag, peak_mag_err", file=f)
            for row in tb_combined:
                print(f"{row['objname']}, "
                      f"{row['origin']}, "
                      f"{row['subtype']}, "
                      f"{row['algo']}, "
                      f"{row['ra']}, "
                      f"{row['dec']}, "
                      f"{row['z']}, {row['z_cmb']}, "
                      f"{row['MJDs']}, {row['MJDe']}, "
                      f"{row['mu']}, {row['mu_err']}, "
                      f"{row['hostMass']}, {row['hostMass_err']}, "
                      f"{row['Tmax']}, {row['Tmax_err']}, "
                      f"{row['stretch']}, {row['stretch_err']}, "
                      f"{row['color']}, {row['color_err']}, "
                      f"{row['amplitude']}, {row['amplitude_err']}, "
                      f"{row['peak_mag']}, {row['peak_mag_err']}", file=f)

    return limit, limit_num_prev-limit_num
def selection_criteria(snpy_path: str = '', salt_path: str = '', save_loc: str = '', criteria: dict or None = None):
    save_toggle = [True, True, True]
    if len(snpy_path) > 0 and len(salt_path) > 0:
        tb_snpy = utils.default_open(snpy_path, True)
        tb_salt = utils.default_open(salt_path, True)
    elif len(snpy_path) > 0:
        save_toggle = [False, True, False]
        tb_snpy = utils.default_open(snpy_path, True)
        tb_salt = Table(names=tb_snpy.colnames)
    elif len(salt_path) > 0:
        save_toggle = [False, False, True]
        tb_salt = utils.default_open(salt_path, True)
        tb_snpy = Table(names=tb_salt.colnames)

    # tb_snpy = utils.default_open(snpy_path, True)
    # tb_salt = utils.default_open(salt_path, True)

    if criteria == None:
        criteria = {
            'z_cmb': [0.015, 999],
            'Tmax_err': [-999, 1.0],
            'mu_err': [-999, 0.2],
            'EBVhost': [-0.3, 0.3],
            'EBVhost_err': [-999, 0.1],
            'st': [-999, 1.0],
            'st_err': [-999, 0.1],
            'c': [-0.6, 0.6],
            'c_err': [-999, 0.1],
            'x1': [-3.2, 3.2],
            'x1_err': [-999, 1.0]
        }

    # Base Cuts
    for c in ['z_cmb', 'mu_err', 'Tmax_err']:
        print(f"[{criteria[c][0]} < {c} < {criteria[c][1]}]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
        tb_snpy = tb_snpy[(tb_snpy[c] > criteria[c][0]) & (tb_snpy[c] < criteria[c][1])]
        tb_salt = tb_salt[(tb_salt[c] > criteria[c][0]) & (tb_salt[c] < criteria[c][1])]
        print(len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname']))))

    # Visual Inspection
    # bad_snpy_sne = "2008bd, 2019ahh, 2019ecx, 2022ubt, 2022vse, 2022vse, 2022xhh, 2023jah, 2023mkp, 2024yhg".split(', ')
    # bad_salt_sne = "2006bd, 2006gt, 2007N, 2008bd, 2008bt, 2021mab, 2022vse, 2022vxf, 2023dk, 2025nn".split(', ')
    # print(f"[Visual Cuts]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
    # for name in bad_snpy_sne:
    #     if name in tb_snpy['objname']: tb_snpy.remove_row(np.where(tb_snpy["objname"] == name)[0][0])
    # for name in bad_salt_sne:
    #     if name in tb_salt['objname']: tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])
    # print(len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname']))))

    # SNooPy Cuts
    for c1, c2 in zip(['color', 'color_err', 'stretch', 'stretch_err'], ['EBVhost', 'EBVhost_err', 'st', 'st_err']):
        print(f"[{criteria[c2][0]} < {c2} < {criteria[c2][1]}]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
        tb_snpy = tb_snpy[(tb_snpy[c1] > criteria[c2][0]) & (tb_snpy[c1] < criteria[c2][1])]
        print(len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname']))))

    # SALT Cuts
    for c1, c2 in zip(['color', 'color_err', 'stretch', 'stretch_err'], ['c', 'c_err', 'x1', 'x1_err']):
        print(f"[{criteria[c2][0]} < {c2} < {criteria[c2][1]}]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
        tb_salt = tb_salt[(tb_salt[c1] > criteria[c2][0]) & (tb_salt[c1] < criteria[c2][1])]
        print(len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname']))))

    # Combine on mu
    if all(save_toggle) == True:
        for name in np.unique(list(tb_snpy["objname"])+list(tb_salt["objname"])):
            if (name in list(tb_snpy["objname"])) & (name in list(tb_salt["objname"])):
                # # Prefer SNooPy MU's
                # tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])

                # Take SNooPy data, remove SALT
                if tb_snpy[tb_snpy["objname"] == name]["mu_err"] < tb_salt[tb_salt["objname"] == name]["mu_err"]:
                    tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])
                # Take SALT data, remove SNooPy
                else:
                    tb_snpy.remove_row(np.where(tb_snpy["objname"] == name)[0][0])

    # Combine data
    if len(tb_snpy) > len(tb_salt):
        tb_combined = tb_snpy.copy()
        for row in tb_salt:
            tb_combined.add_row(row)
    else:
        tb_combined = tb_salt.copy()
        for row in tb_snpy:
            tb_combined.add_row(row)

    # Chauvenet’s Criterion
    resid = tb_combined['mu'].astype(float) - utils.current_cosmo().distmod(tb_combined['z_cmb'].astype(float)).value
    limit = norm.ppf(1 - (1.0 / (4 * len(resid))))
    tb_combined = tb_combined[np.abs((resid - np.average(resid, weights=(1/(tb_combined['mu_err'].astype(float)**2)))) / np.std(resid)) < limit]
    og_len = len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))
    print(f"[Chauvenet: {limit:.2f}]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
    tb_snpy = tb_combined[tb_combined['algo'] == 'SNPY']
    tb_salt = tb_combined[tb_combined['algo'] == 'SALT']
    print(len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname']))))
    chuv_lim, chuv_lim_num = limit, og_len-len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname'])))

    # Adding intrinsic dispersion (0.1 mag) in quadrature for mass (taylor+11) & mu
    tb_combined['mu_err'] = np.sqrt(tb_combined['mu_err'] ** 2.0 + 0.1 ** 2.0)
    tb_snpy['mu_err'] = np.sqrt(tb_snpy['mu_err'] ** 2.0 + 0.1 ** 2.0)
    tb_salt['mu_err'] = np.sqrt(tb_salt['mu_err'] ** 2.0 + 0.1 ** 2.0)
    tb_combined['hostMass_err'] = np.sqrt(tb_combined['hostMass_err'] ** 2.0 + 0.1 ** 2.0)
    tb_snpy['hostMass_err'] = np.sqrt(tb_snpy['hostMass_err'] ** 2.0 + 0.1 ** 2.0)
    tb_salt['hostMass_err'] = np.sqrt(tb_salt['hostMass_err'] ** 2.0 + 0.1 ** 2.0)

    # Get Mass Step
    tb_low = tb_combined[tb_combined['hostMass'] < 10]
    tb_high = tb_combined[tb_combined['hostMass'] > 10]
    resid = tb_combined['mu'].astype(float) - utils.current_cosmo().distmod(tb_combined['z_cmb'].astype(float)).value
    resid_high = tb_high['mu'].astype(float) - utils.current_cosmo().distmod(tb_high['z_cmb'].astype(float)).value
    resid_low = tb_low['mu'].astype(float) - utils.current_cosmo().distmod(tb_low['z_cmb'].astype(float)).value
    mass_step = (np.average(resid_high, weights=(1/(tb_high['mu_err']**2))) -
                 np.average(resid_low,  weights=(1/(tb_low['mu_err']**2))))
    mass_step_err = np.sqrt(((np.std(resid_low) / np.sqrt(len(resid_low))) ** 2) +
                            ((np.std(resid_high) / len(resid_high))) ** 2) # (sigma/N)**2

    # Scatter
    print(f"[Scatter]: {np.std(resid):.3f}")

    # Nums
    print(f"[# SNooPy]: {len(tb_snpy)}")
    print(f"[# SALT3]: {len(tb_salt)}")

    # Write data to files
    for tb, save_loc, comments, st in zip([tb_combined, tb_snpy, tb_salt],
                                          [save_loc, f"{snpy_path[:-4]}_cut.txt", f"{salt_path[:-4]}_cut.txt"],
                                          [[f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_combined)}",
                                            f"# Scatter: {np.std(resid):.4f}, MASS_STEP: {mass_step:.4f}+/-{mass_step_err:.4f}",
                                            f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical."],
                                          [f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_snpy)}"],
                                          [f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_salt)}"]],
                                          save_toggle):
        if st:
            with open(save_loc, 'w') as f:
                for comment in comments:
                    print(comment, file=f)
                hdr = list(tb.columns)
                print(str(str(hdr).replace("'", "")[1:-1]), file=f)
                for row in tb:
                    line = ""
                    for c in hdr:
                        line += f"{row[c]}, "
                    print(line[:-2], file=f)
            print(f"[+++] Saved data to {save_loc}...")
    return chuv_lim, chuv_lim_num
def visual_inspection(path: str, save_loc: str = ''):
    print(f"[+++] Visual inspection for '{path}'...")

    # Load data
    hdr = []
    with open(path, 'r') as f:
        for i in range(3): hdr.append(f.readline())
        lines = f.readlines()


    bad_csp = "2006bd, 2007N, 2008bd, 2008bt".split(', ')
    no_atlas_ztf = "2023jah, 2018hkw, 2021gel, 2022ubt, 2022vse, 2022vxf, 2022yv, 2023dk, 2024jhk, 2025nn".split(', ')
    selected_cuts = bad_csp.copy() + no_atlas_ztf.copy()

    with open(save_loc, 'w') as f:
        for line in hdr:
            print(line[:-1], file=f)
        for line in lines:
            n = line.split(', ')[0]
            if n not in selected_cuts:
                print(line[:-1], file=f)
    print(f"      SNe removed ({len(selected_cuts)}): {selected_cuts}")
    return
def old_seperate_cut_data(mixed_path: str, snpy_path: str, salt_path: str, snpy_path_new: str, salt_path_new: str):
    print(f"[~~~] Spliting data from '{mixed_path}' into SNooPy & SALT3 cut files...")
    tb_mixed = utils.default_open(mixed_path, True)
    tb_snpy = utils.default_open(snpy_path, True)
    tb_salt = utils.default_open(salt_path, True)

    tb_snpy_new = Table(names=tb_mixed.colnames, dtype=tb_mixed.dtype)
    tb_salt_new = Table(names=tb_mixed.colnames, dtype=tb_mixed.dtype)
    for name in list(tb_mixed['objname']):
        if name in list(tb_snpy['objname']):
            row = tb_snpy[tb_snpy['objname'] == name][0]
            mixed_row = tb_mixed[tb_mixed['objname'] == name][0]
            row['mu'] = mixed_row['mu']
            row['mu_err'] = mixed_row['mu_err']
            tb_snpy_new.add_row(row)
        if name in list(tb_salt['objname']):
            row = tb_salt[tb_salt['objname'] == name][0]
            mixed_row = tb_mixed[tb_mixed['objname'] == name][0]
            row['mu'] = mixed_row['mu']
            row['mu_err'] = mixed_row['mu_err']
            tb_salt_new.add_row(row)

    # Save SNPY file
    print(f"      Updating '{snpy_path_new}' ({len(tb_snpy_new)}/{len(tb_mixed)})...")
    with open(snpy_path_new, 'w') as f:
        print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_snpy_new)}", file=f)
        print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
        print(str(tb_snpy_new.colnames)[1:-1].replace("'", ""), file=f)
        for row in tb_snpy_new:
            print(str(list(row))[1:-1].replace("'", ""), file=f)

    # Save SALT file
    print(f"      Updating '{salt_path_new}' ({len(tb_salt_new)}/{len(tb_mixed)})...")
    with open(salt_path_new, 'w') as f:
        print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_salt_new)}", file=f)
        print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
        print(str(tb_salt_new.colnames)[1:-1].replace("'", ""), file=f)
        for row in tb_salt_new:
            print(str(list(row))[1:-1].replace("'", ""), file=f)
    return
def seperate_cut_data(mixed_path: str, snpy_path: str, salt_path: str):
    print(f"[~~~] Spliting data from '{mixed_path}' into SNooPy & SALT3 cut files...")
    tb_mixed = utils.default_open(mixed_path, True)

    tb_snpy = tb_mixed[tb_mixed['algo'] == 'SNPY']
    tb_salt = tb_mixed[tb_mixed['algo'] == 'SALT']

    # Save SNPY file
    print(f"      Updating '{snpy_path}' ({len(tb_snpy)}/{len(tb_mixed)})...")
    with open(snpy_path, 'w') as f:
        print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_snpy)}", file=f)
        print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
        print(str(tb_snpy.colnames)[1:-1].replace("'", ""), file=f)
        for row in tb_snpy:
            print(str(list(row))[1:-1].replace("'", ""), file=f)

    # Save SALT file
    print(f"      Updating '{salt_path}' ({len(tb_salt)}/{len(tb_mixed)})...")
    with open(salt_path, 'w') as f:
        print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_salt)}", file=f)
        print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
        print(str(tb_salt.colnames)[1:-1].replace("'", ""), file=f)
        for row in tb_salt:
            print(str(list(row))[1:-1].replace("'", ""), file=f)
    return
def save_tb(tb: Table, save_loc: str = '', comments: list[str] = []):
    with open(save_loc, 'w') as f:
        for comment in comments:
            print(comment, file=f)
        hdr = list(tb.columns)
        print(str(str(hdr).replace("'", "")[1:-1]), file=f)
        for row in tb:
            line = ""
            for c in hdr:
                line += f"{row[c]}, "
            print(line[:-2], file=f)
    print(f"[+++] Saved data to {save_loc}...")
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
