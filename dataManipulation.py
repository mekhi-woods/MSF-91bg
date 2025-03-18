import os
import glob
import utils  # Import of utils.py
import shutil
import datetime
import numpy as np
import time as systime
from astropy.table import Table, vstack
from astropy.stats import sigma_clip, sigma_clipped_stats

CURRENTDATE = datetime.datetime.now()


def make_objTbl(source: str, path: str) -> Table:
    # Load data and make intial table
    if source == 'atlas':
        # Load data
        with open(path, 'r') as f:
            hdr = f.readline()[1:-1].split(',')
        data = np.genfromtxt(path, dtype='str', delimiter=',', skip_header=1)

        # Make table
        var_table = Table()
        for h in hdr:
            try:
                var_table[h] = data[:, hdr.index(h)].astype(float)
            except ValueError:
                var_table[h] = data[:, hdr.index(h)]

        # Add parity with CSP & ZTF
        var_table.remove_columns(['err', 'chi/N', 'x', 'y', 'maj', 'min', 'phi', 'apfit', 'mag5sig', 'Sky', 'Obs'])
        var_table['zp'] = np.full(len(var_table), np.nan)
        var_table['z'] = np.full(len(var_table), np.nan)
        for h_old, h_new in zip(['JD', 'm', 'dm', 'uJy', 'duJy', 'F', 'RA', 'Dec'],
                                ['mjd', 'mag', 'dmag', 'flux', 'dflux', 'filter', 'ra', 'dec']):
            var_table[h_old].name = h_new
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
    elif source == 'ztf':
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

        # Add temp redshift
        var_table['z'] = np.full(len(var_table), np.nan)

        # Fix time, JD to MJD
        var_table['jd'] = np.array(var_table['jd']).astype(float) - 2400000.5
        var_table['forcediffimflux'][var_table['forcediffimflux'] == 'null'] = 'nan'
        var_table['forcediffimfluxunc'][var_table['forcediffimfluxunc'] == 'null'] = 'nan'
        var_table['mag'] = (-2.5 * np.log10(np.array(var_table['forcediffimflux']).astype(float))) + var_table['zpdiff']
        var_table['dmag'] = np.abs(-1.08573620476 * (np.array(var_table['forcediffimfluxunc']).astype(float)
                                                     / np.array(var_table['forcediffimflux']).astype(float)))

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
    else:
        raise ValueError(f'[!!!] Unknown source, {source}! Must be [atlas/csp/ztf]...')
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
def old_ztf_combine_like_data(target_file: str, combined_loc: str, subtype: str, clear_old: bool = False):
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
            n_tb = Table(names=['filter', 'zp', 'z', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux', 'source'],
                         dtype=[str, float, float, float, float, float, float, float, float, float, str])
            for source in row['Sources'].split('+'):
                var_table = make_objTbl(source.lower(), f"data/{source}-{subtype}/{source}{name}.txt")
                for n_row in var_table:
                    n_tb.add_row([n_row['filter'], n_row['zp'], n_row['z'], n_row['ra'], n_row['dec'], n_row['mjd'],
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
def combine_like_data(target_file: str, combined_loc: str, subtype: str, clear_old: bool = False):
    # Remove old combined data
    if clear_old:
        if os.path.exists(combined_loc):
            shutil.rmtree(combined_loc)
        os.mkdir(combined_loc)

    # Combine all files into one list
    all_files = (glob.glob(f'data/CSP-{subtype}/*') +
                 glob.glob(f'data/ATLAS-{subtype}/*') +
                 glob.glob(f'data/ZTF-{subtype}/*'))

    # Open target list
    tarlist_data = np.genfromtxt(target_file, delimiter=',', dtype=str)
    tarlist = Table(names=tarlist_data[0, :], data=tarlist_data[1:, :])

    for row in tarlist:
        print(row['Name'])

        break


    # for file in all_files:
    #     f_origin = file.split('/')[-1][:-4].split('2')[0]
    #     f_name = "SN " + (file.split('/')[-1][:-4].
    #               replace('CSP', '').replace('ATLAS', '').replace('ZTF', ''))
    #     f_tar_index = np.where(tarlist['Name'] == f_name)[0][0]
    #     print(tarlist[f_tar_index])
    #     break

    #
    # # With that list, get the names and origins of data
    # names_orgins = Table(names=['name', 'origin'], dtype=[str, str])
    # for file in all_files:
    #     origin = file.split('/')[-2].split('-')[0]
    #     name = file.split('/')[-1].split('.')[0].replace(origin, '')
    #     names_orgins.add_row([name, origin])
    #
    # # Open target file
    # target_file = np.genfromtxt(target_file, delimiter=',', dtype=str, skip_header=1)
    # tb_targets = Table(names=target_file[0], data=target_file[1:])
    #
    # # Go through list of targets, check if they exsist in data, and take the name of the surveys they're in
    # tb_combined = Table(names=list(target_file[0])+['Sources'], dtype=[str]*(len(list(target_file[0]))+1))
    # for i, row in enumerate(tb_targets):
    #     all_sources = ''
    #     tb_overlaping_sources = names_orgins[names_orgins['name'] == row['Name'][3:]]
    #     for n_row in tb_overlaping_sources:
    #         all_sources = all_sources+n_row['origin']+'+'
    #     if len(tb_overlaping_sources) == 0:
    #         tb_combined.add_row(list(row) + ['None'])
    #     else:
    #         tb_combined.add_row(list(row) + [all_sources[:-1]])
    #
    # # Go through overlapping targets and combine data files
    # for i, row in enumerate(tb_combined):
    #     name = row['Name'][3:]
    #
    #     # Combine data from >1 source
    #     if row['Sources'] != 'None' and '+' in row['Sources']:
    #         n_tb = Table(names=['filter', 'zp', 'z', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux', 'source'],
    #                      dtype=[str, float, float, float, float, float, float, float, float, float, str])
    #         for source in row['Sources'].split('+'):
    #             var_table = make_objTbl(source.lower(), f"data/{source}-{subtype}/{source}{name}.txt")
    #             for n_row in var_table:
    #                 n_tb.add_row([n_row['filter'], n_row['zp'], n_row['z'], n_row['ra'], n_row['dec'], n_row['mjd'],
    #                               n_row['mag'], n_row['dmag'], n_row['flux'], n_row['dflux'], source])
    #         with open(f"{combined_loc}{name}.txt", 'w') as f:
    #             file_hdr = str(n_tb.colnames).strip('[]').replace("'", "")
    #             print(f"# Combination of '{name}' data from {str(row['Sources'].split('+')).strip('[]')}", file=f)
    #             print(file_hdr, file=f)
    #             for r in n_tb:
    #                 print(str(list(r)).strip('[]').replace("'", ""), file=f)
    #
    #     # Take copy for later
    #     elif row['Sources'] != 'None':
    #       shutil.copyfile(f"data/{row['Sources']}-{subtype}/{row['Sources']}{name}.txt",
    #                       f"{combined_loc}{row['Sources']}{name}.txt")

    return
def combine_snpy_salt(snpy_path: str, salt_path: str, save_loc: str = ''):
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
            if sn.path.split('_')[-2].upper() == "SNPY":
                param_dict = {'t_max': 'Tmax', 'stretch': 'st', 'color': 'EBVhost'}
            elif sn.path.split('_')[-2].upper() == "SALT":
                param_dict = {'t_max': 't0', 'stretch': 'x1', 'color': 'c'}
            if sn.path.split('_')[-2].upper() == "SALT":
                amp, amp_err = sn.params['x0']['value'], sn.params['x0']['err']
            else:
                amp, amp_err = np.nan, np.nan
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
def selection_criteria(path: str, subtype: str = '91bg', save_loc: str = ''):
    print(f"[+++] Selection criteria for '{path}'...")
    # 999 = no upper limit, -999 = no lower limit
    if subtype == '91bg':
        # Original cuts
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
            'mu_err':      [-999, 0.2],
        }
        # criteria = {
        #     'z_cmb':        [0.015, 999],
        #     'EBVhost':      [-0.2, 0.3],
        #     'EBVhost_err':  [-999, 0.1],
        #     'st':           [-999, 1.0],
        #     'st_err':       [-999, 0.1],
        #     'c':            [-0.6, 0.6],
        #     'c_err':        [-999, 0.1],
        #     'x1':           [-3.2, 3.2],
        #     'x1_err':       [-999, 1.0],
        #     'Tmax_err':     [-999, 1.0],
        #     't0_err':       [-999, 1.0],
        #     'mu_err':       [-999, 0.2],
        # }
    elif subtype == 'norm':
       criteria = {
            'z_cmb':       [0.015, 999],
            'EBVhost':     [-999, 999],
            'EBVhost_err': [-999, 1.0],
            'st':          [-999, 999],
            'st_err':      [-999, 1.0],
            'c':           [-999, 999],
            'c_err':       [-999, 1.0],
            'x1':          [-999, 999],
            'x1_err':      [-999, 1.0],
            'Tmax_err':    [-999, 1.0],
            't0_err':      [-999, 1.0],
            'mu_err':      [-999, 0.2],
        }

    # Load data
    tb = utils.default_open(path, True)

    # Sigma Clipping
    resid = (tb['mu'].astype(float) - utils.current_cosmo().distmod(tb['z_cmb'].astype(float)).value)
    clipped_data = sigma_clip(resid, sigma=3, maxiters=5)
    tb_clipped = Table(names=tb.colnames, dtype=[object] * len(tb.colnames))
    for i in range(len(clipped_data)):
        if ~clipped_data.mask[i]:
            tb_clipped.add_row(tb[i])
    resid_clipped = (tb_clipped['mu'].astype(float) - utils.current_cosmo().distmod(tb_clipped['z_cmb'].astype(float)).value)
    print(f"[~~~] Sigma Clipping on Hubble Residual: σ = 3... SNe = {len(resid)} ---> {len(resid_clipped)} "
          f"(Hubble Residual Scatter = {round(np.std(resid), 3)} --> {round(np.std(resid_clipped), 3)})")
    tb = tb_clipped.copy()

    # Basic cuts
    for p in ['z_cmb', 'Tmax_err', 'mu_err']:
        print(f"{p}: {len(tb)} -> ", end='')
        tb = tb[(tb[p] > criteria[p][0]) & (tb[p] <= criteria[p][1])]
        print(f"{len(tb)} ", end='')
        print(f"({np.std(tb['mu'].astype(float) - utils.current_cosmo().distmod(tb['z_cmb'].astype(float)).value)})")

    # Seperate SNPY & SALT
    tb_snpy = tb[tb['algo'] == 'SNPY']
    tb_salt = tb[tb['algo'] == 'SALT']

    # Cut on SNPY params
    for f_name, d_name in zip(['EBVhost', 'st'], ['color', 'stretch']):
        for n_end in ['', '_err']:
            print(f"{f_name+n_end}: {len(tb_snpy)+len(tb_salt)} -> ", end='')
            tb_snpy = tb_snpy[(tb_snpy[d_name+n_end] >= criteria[f_name+n_end][0]) &
                              (tb_snpy[d_name+n_end] <= criteria[f_name+n_end][1])]
            print(f"{len(tb_snpy)+len(tb_salt)} ", end='')

            tb_temp = vstack([tb_snpy, tb_salt])
            print(f"({np.std(tb_temp['mu'].astype(float) - utils.current_cosmo().distmod(tb_temp['z_cmb'].astype(float)).value)})")

    # Cut on SALT params
    for f_name, d_name in zip(['c', 'x1'], ['color', 'stretch']):
        for n_end in ['', '_err']:
            print(f"{f_name+n_end}: {len(tb_snpy)+len(tb_salt)} -> ", end='')
            tb_salt = tb_salt[(tb_salt[d_name+n_end] >= criteria[f_name+n_end][0]) &
                              (tb_salt[d_name+n_end] <= criteria[f_name+n_end][1])]
            print(f"{len(tb_snpy)+len(tb_salt)} ", end='')

            tb_temp = vstack([tb_snpy, tb_salt])
            print(f"({np.std(tb_temp['mu'].astype(float) - utils.current_cosmo().distmod(tb_temp['z_cmb'].astype(float)).value)})")

    # Recombine SNPY & SALT
    tb_combined = vstack([tb_snpy, tb_salt])

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

    # Old cutting
    # # Load data & Seperate tables
    # tb = utils.default_open(path, True)
    #
    # # Preform general cuts
    # for cut in ['z_cmb', 'mu_err']:
    #     print(f"      Selecting on {cut}: {criteria[cut][0]} < {cut} < {criteria[cut][1]}... ", end='')
    #     print(f"SNe = {len(tb)} ---> ", end='')
    #     tb = tb[(tb[cut] >= criteria[cut][0]) &
    #             (tb[cut] <= criteria[cut][1])]
    #     print(len(tb), end='\n')
    #     # print(f"({round(min(tb[cut]), 4)} < "
    #     #       f"{round(np.median(tb[cut]), 4)} < "
    #     #       f"{round(max(tb[cut]), 4)})",
    #     #       end='\n')
    #
    # # Preform algorithm specific cuts
    # tb_snpy = tb[tb['algo'] == 'SNPY']
    # tb_salt = tb[tb['algo'] == 'SALT']
    # for snpy_param, salt_param, hdr_name in zip(['c', 'c_err', 'x1', 'x1_err', 't0_err'],
    #                                             ['EBVhost', 'EBVhost_err', 'st', 'st_err', 'Tmax_err'],
    #                                             ['color', 'color_err', 'stretch', 'stretch_err', 'Tmax_err']):
    #     # Cut on SNooPy Parameter
    #     print(f"      Selecting on {hdr_name}: {criteria[snpy_param][0]} < {snpy_param} < {criteria[snpy_param][1]}... ", end='')
    #     print(f"SNe = {len(tb_snpy) + len(tb_salt)} ---> ", end='')
    #     tb_snpy = tb_snpy[(tb_snpy[hdr_name] >= criteria[snpy_param][0]) &
    #                       (tb_snpy[hdr_name] <= criteria[snpy_param][1])]
    #     print(len(tb_snpy) + len(tb_salt), end='\n')
    #     # print(f"({round(min(tb_snpy[hdr_name]), 4)} < "
    #     #       f"{round(np.median(tb_snpy[hdr_name]), 4)} < "
    #     #       f"{round(max(tb_snpy[hdr_name]), 4)})",
    #     #       end='\n')
    #
    #     # Cut on SALT Parameter
    #     print(f"      Selecting on {hdr_name}: {criteria[salt_param][0]} < {salt_param} < {criteria[salt_param][1]}... ", end='')
    #     print(f"SNe = {len(tb_snpy) + len(tb_salt)} ---> ", end='')
    #     tb_salt = tb_salt[(tb_salt[hdr_name] >= criteria[salt_param][0]) &
    #                       (tb_salt[hdr_name] <= criteria[salt_param][1])]
    #     print(len(tb_snpy) + len(tb_salt), end='\n')
    #     # print(f"({round(min(tb_salt[hdr_name]), 4)} < "
    #     #       f"{round(np.median(tb_salt[hdr_name]), 4)} < "
    #     #       f"{round(max(tb_salt[hdr_name]), 4)})",
    #     #       end='\n')
    #
    # # Combine tables
    # tb_combined = tb_salt.copy()
    # for row in tb_snpy:
    #     tb_combined.add_row(row)
    #
    # # Sigma Clipping
    # resid = (tb_combined['mu'].astype(float) - utils.current_cosmo().distmod(tb_combined['z_cmb'].astype(float)).value)
    # clipped_data = sigma_clip(resid, sigma=3, maxiters=5)
    # tb_clipped = Table(names=tb_combined.colnames, dtype=[object]*len(tb_combined.colnames))
    # for i in range(len(clipped_data)):
    #     if ~clipped_data.mask[i]:
    #         tb_clipped.add_row(tb_combined[i])
    # resid_clipped = (tb_clipped['mu'].astype(float) - utils.current_cosmo().distmod(tb_clipped['z_cmb'].astype(float)).value)
    # print(f"      Sigma Clipping on Hubble Residual: σ = 3... SNe = {len(resid)} ---> {len(resid_clipped)} "
    #       f"(Hubble Residual Scatter = {round(np.std(resid), 3)} --> {round(np.std(resid_clipped), 3)})")
    #
    # # old sigma clipping
    # # r_mn, r_md, r_std = sigma_clipped_stats(resid)
    # # print(f'      Pre-Sigma-Clipping Hubble Residual Scatter: {round(r_std, 3)}')
    # # tb_combined = tb_combined[abs(resid - r_mn) < 3 * r_std]
    # # new_resid = (tb_combined['mu'].astype(float) - utils.current_cosmo().distmod(tb_combined['z_cmb'].astype(float)).value)
    # # mn, md, std = sigma_clipped_stats(new_resid)
    # # print(f'      Post-Sigma-Clipping (3-sigma) Hubble Residual Scatter: {round(std, 3)}')
    #
    # # Save new parameter file
    # if len(save_loc) != 0:
    #     print(f"[~~~] Saving new parameter file to '{save_loc}'...")
    #     with open(save_loc, 'w') as f:
    #         print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_snpy) + len(tb_salt)}", file=f)
    #         print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
    #         print(f"objname, origin, subtype, algo, ra, dec, z, z_cmb, MJDs, MJDe, mu, mu_err, hostMass, hostMass_err, "
    #               f"Tmax, Tmax_err, stretch, stretch_err, color, color_err, amplitude, amplitude_err, "
    #               f"peak_mag, peak_mag_err", file=f)
    #         for row in tb_clipped:
    #             print(f"{row['objname']}, "
    #                   f"{row['origin']}, "
    #                   f"{row['subtype']}, "
    #                   f"{row['algo']}, "
    #                   f"{row['ra']}, "
    #                   f"{row['dec']}, "
    #                   f"{row['z']}, {row['z_cmb']}, "
    #                   f"{row['MJDs']}, {row['MJDe']}, "
    #                   f"{row['mu']}, {row['mu_err']}, "
    #                   f"{row['hostMass']}, {row['hostMass_err']}, "
    #                   f"{row['Tmax']}, {row['Tmax_err']}, "
    #                   f"{row['stretch']}, {row['stretch_err']}, "
    #                   f"{row['color']}, {row['color_err']}, "
    #                   f"{row['amplitude']}, {row['amplitude_err']}, "
    #                   f"{row['peak_mag']}, {row['peak_mag_err']}", file=f)

    # # old Save new parameter file
    # if len(save_loc) != 0:
    #     print(f"[~~~] Saving new parameter file to '{save_loc}'...")
    #     with open(save_loc, 'w') as f:
    #         print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_snpy) + len(tb_salt)}", file=f)
    #         print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
    #         print(f"objname, origin, subtype, algo, ra, dec, z, z_cmb, MJDs, MJDe, mu, mu_err, hostMass, hostMass_err, "
    #               f"Tmax, Tmax_err, stretch, stretch_err, color, color_err, amplitude, amplitude_err, "
    #               f"peak_mag, peak_mag_err", file=f)
    #         for tb in [tb_snpy, tb_salt]:
    #             for row in tb:
    #                 print(f"{row['objname']}, "
    #                       f"{row['origin']}, "
    #                       f"{row['subtype']}, "
    #                       f"{row['algo']}, "
    #                       f"{row['ra']}, "
    #                       f"{row['dec']}, "
    #                       f"{row['z']}, {row['z_cmb']}, "
    #                       f"{row['MJDs']}, {row['MJDe']}, "
    #                       f"{row['mu']}, {row['mu_err']}, "
    #                       f"{row['hostMass']}, {row['hostMass_err']}, "
    #                       f"{row['Tmax']}, {row['Tmax_err']}, "
    #                       f"{row['stretch']}, {row['stretch_err']}, "
    #                       f"{row['color']}, {row['color_err']}, "
    #                       f"{row['amplitude']}, {row['amplitude_err']}, "
    #                       f"{row['peak_mag']}, {row['peak_mag_err']}", file=f)
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
