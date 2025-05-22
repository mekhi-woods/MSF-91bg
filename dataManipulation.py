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

def make_param_file(sne: list, subtype: str, save_loc: str):
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
              "z_err, "
              "z_cmb, "
              "MJDs, "
              "MJDe, "
              # "chisquare, chisquare_err, "
              "mu, mu_err, "
              "hostMass, hostMass_err, "
              "Tmax, Tmax_err, "
              "stretch, stretch_err, "
              "color, color_err, "
              "amplitude, amplitude_err, "
              "peak_mag, peak_mag_err", file=f)
        for sn in sne:
            if sn.algo == 'snpy':
                param_dict = {'t_max': 'Tmax', 'stretch': 'st', 'color': 'EBVhost'}
                amp, amp_err = np.nan, np.nan
            elif sn.algo == 'salt':
                param_dict = {'t_max': 't0', 'stretch': 'x1', 'color': 'c'}
                amp, amp_err = sn.params['x0']['value'], sn.params['x0']['err']
            if sn.z_err == None or sn.z_err == 'None':
                sn.z_err = 0.000
            print(f"{sn.objname}, "
                  f"{sn.origin}, "
                  f"{subtype}, "
                  f"{sn.algo.upper()}, "
                  f"{sn.coords[0]}, "
                  f"{sn.coords[1]}, "
                  f"{sn.z}, "
                  f"{sn.z_err}, "
                  f"{sn.z_cmb}, "
                  f"{min(sn.mjd)}, "
                  f"{max(sn.mjd)}, "
                  # f"{sn.params['chisquare']['value']}, {sn.params['chisquare']['err']}, "
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
def combine_snpy_salt(snpy_path: str, salt_path: str, save_loc: str = ''):
    tb_snpy = utils.default_open(snpy_path, True)
    tb_salt = utils.default_open(salt_path, True)
    for name in np.unique(list(tb_snpy["objname"])+list(tb_salt["objname"])):
        if (name in list(tb_snpy["objname"])) & (name in list(tb_salt["objname"])):
            # # Prefer SNooPy MU's
            # tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])

            # Reduced Chi2 Comparision
            # Take SNooPy data, remove SALT
            # if tb_snpy[tb_snpy["objname"] == name]["chisquare"] < tb_salt[tb_salt["objname"] == name]["chisquare"]:
            #     tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])
            # # Take SALT data, remove SNooPy
            # else:
            #     tb_snpy.remove_row(np.where(tb_snpy["objname"] == name)[0][0])

            # MU_err Comparision
            # Take SNooPy data, remove SALT
            if tb_snpy[tb_snpy["objname"] == name]["mu_err"] < tb_salt[tb_salt["objname"] == name]["mu_err"]:
                tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])
            # Take SALT data, remove SNooPy
            else:
                tb_snpy.remove_row(np.where(tb_snpy["objname"] == name)[0][0])

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

    if criteria == None:
        criteria = {
            'z_cmb': [0.015, 999],
            'z_err': [-999, 999],
            'Tmax_err': [-999, 1.0],
            'mu_err': [-999, 0.2],
            'chisquare': [-999, 999],
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
    for c in ['z_cmb', 'mu_err', 'Tmax_err']: # 'z_err', 'chisquare'
        print(f"[{criteria[c][0]} < {c} < {criteria[c][1]}]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
        tb_snpy = tb_snpy[(tb_snpy[c] > criteria[c][0]) & (tb_snpy[c] < criteria[c][1])]
        tb_salt = tb_salt[(tb_salt[c] > criteria[c][0]) & (tb_salt[c] < criteria[c][1])]
        # print(len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname']))))
        print(f"{len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname'])))} [SNooPy: {len(list(tb_snpy['objname']))}| SALT: {len(list(tb_salt['objname']))}]")

    # Visual Inspection
    # bad_snpy_sne = "2008bd, 2019ahh, 2019ecx, 2022ubt, 2022vse, 2022vse, 2022xhh, 2023jah, 2023mkp, 2024yhg".split(', ')
    # bad_salt_sne = "2006bd, 2006gt, 2007N, 2008bd, 2008bt, 2021mab, 2022vse, 2022vxf, 2023dk, 2025nn".split(', ')
    bad_snpy_sne = "2008bd,2007ba,2022vse,2023jah,2023mkp,2024yhg,2019ecx,2006bd,2007N,2025nn,2021mab,2023dk,2022vxf".split(',') # 2022zsp,2019ahh,2024zaj
    bad_salt_sne = "2008bd,2007ba,2022vse,2023jah,2023mkp,2024yhg,2019ecx,2006bd,2007N,2025nn,2021mab,2023dk,2022vxf".split(',')
    print(f"[Visual Cuts]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
    for name in bad_snpy_sne:
        if name in tb_snpy['objname']: tb_snpy.remove_row(np.where(tb_snpy["objname"] == name)[0][0])
    for name in bad_salt_sne:
        if name in tb_salt['objname']: tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])
    print(len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname']))))

    # SNooPy Cuts
    for c1, c2 in zip(['color', 'color_err', 'stretch', 'stretch_err'], ['EBVhost', 'EBVhost_err', 'st', 'st_err']):
        print(f"[{criteria[c2][0]} < {c2} < {criteria[c2][1]}]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
        tb_snpy = tb_snpy[(tb_snpy[c1] > criteria[c2][0]) & (tb_snpy[c1] < criteria[c2][1])]
        # print(len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname']))))
        print(f"{len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname'])))} [SNooPy: {len(list(tb_snpy['objname']))}]")

    # SALT Cuts
    for c1, c2 in zip(['color', 'color_err', 'stretch', 'stretch_err'], ['c', 'c_err', 'x1', 'x1_err']):
        print(f"[{criteria[c2][0]} < {c2} < {criteria[c2][1]}]: {len(np.unique(list(tb_snpy['objname'])+list(tb_salt['objname'])))}", end=' ---> ')
        tb_salt = tb_salt[(tb_salt[c1] > criteria[c2][0]) & (tb_salt[c1] < criteria[c2][1])]
        # print(len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname']))))
        print(f"{len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname'])))} [SALT: {len(list(tb_salt['objname']))}]")

    # Combine on mu
    if all(save_toggle) == True:
        for name in np.unique(list(tb_snpy["objname"])+list(tb_salt["objname"])):
            if (name in list(tb_snpy["objname"])) & (name in list(tb_salt["objname"])):
                # # Prefer SNooPy MU's
                # tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])

                # MU_err Comparision
                # Take SNooPy data, remove SALT
                if tb_snpy[tb_snpy["objname"] == name]["mu_err"] < tb_salt[tb_salt["objname"] == name]["mu_err"]:
                    tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])
                # Take SALT data, remove SNooPy
                else:
                    tb_snpy.remove_row(np.where(tb_snpy["objname"] == name)[0][0])

                # # Reduced Chi2 Comparision
                # # Take SNooPy data, remove SALT
                # if tb_snpy[tb_snpy["objname"] == name]["chisquare"] < tb_salt[tb_salt["objname"] == name]["chisquare"]:
                #     tb_salt.remove_row(np.where(tb_salt["objname"] == name)[0][0])
                # # Take SALT data, remove SNooPy
                # else:
                #     tb_snpy.remove_row(np.where(tb_snpy["objname"] == name)[0][0])
    print(f"Combining on MU [SNooPy: {len(list(tb_snpy['objname']))}| SALT: {len(list(tb_salt['objname']))}]")

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
    # print(len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname']))))
    print(f"{len(np.unique(list(tb_snpy['objname']) + list(tb_salt['objname'])))} [SNooPy: {len(list(tb_snpy['objname']))}| SALT: {len(list(tb_salt['objname']))}]")

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
def check_peak(sne: object, tol: float = 10.0, limit: int = 10):
    valid_sne, invalid_sne = [], []
    for sn in sne:
        if sn.algo == 'snpy': tmax_key = 'Tmax'
        else: tmax_key = 't0'
        flux_in_range = sn.mag[(sn.mjd.astype(float) < sn.params[tmax_key]['value'] + tol) &
                               (sn.mjd.astype(float) > sn.params[tmax_key]['value'] - tol)]
        if len(flux_in_range) >= limit:
            valid_sne.append(sn)
        else:
            invalid_sne.append(sn.objname)
    print(f"===========================================\n"
          f"{len(valid_sne)}/{len(sne)} SNe with '{limit}' data points around Tmax \n"
          f"Missing: {invalid_sne}\n"
          f"===========================================")
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
