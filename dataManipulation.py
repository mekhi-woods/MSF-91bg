import os
import glob
import utils  # Import of utils.py
import shutil
import datetime
import numpy as np
import time as systime
from astropy.table import Table

CURRENTDATE = datetime.datetime.now()

def combine_like_data(atlas_loc: str, ztf_loc: str, atlas_ztf_loc: str, clear_old: bool = False):
    """
    :param atlas_loc: Location of atlas data set
    :param ztf_loc: Location of ztf data set
    :param atlas_ztf_loc: Location of where to store combined data set
    :return:
    """
    data_atlas_unique_loc = atlas_loc + "-unique/"
    data_ztf_unique_loc = ztf_loc + "-unique/"

    # Clear old combination data
    if clear_old:
        print(f"Removing old ATLAS+ZTF combinations..."
              f"\n\t{atlas_ztf_loc}\n\t{data_atlas_unique_loc}\n\t{data_ztf_unique_loc}")
        systime.sleep(3)
        shutil.rmtree(atlas_ztf_loc, ignore_errors=True)
        os.makedirs(atlas_ztf_loc)
        shutil.rmtree(data_atlas_unique_loc, ignore_errors=True)
        os.makedirs(data_atlas_unique_loc)
        shutil.rmtree(data_ztf_unique_loc, ignore_errors=True)
        os.makedirs(data_ztf_unique_loc)

    # Get RA & DEC from ATLAS
    atlas_ra_dec = Table(names=["RA", "DEC", "path"], dtype=[float, float, str])
    for file in glob.glob(atlas_loc + "/*"):
        with open(file, 'r') as f:
            hdr = f.readline().rstrip('\n').split(',')
            tb = Table(names=hdr, dtype=[str]*len(hdr))
            for line in f.readlines():
                tb.add_row(line.rstrip('\n').split(','))
            ra, dec = np.average(np.array(tb['RA']).astype(float)), np.average(np.array(tb['Dec']).astype(float))
            atlas_ra_dec.add_row([ra, dec, file])

    # Get RA & DEC from ZTF
    ztf_ra_dec = Table(names=["RA", "DEC", "path"], dtype=[float, float, str])
    for file in glob.glob(ztf_loc + "/*"):
        with open(file, 'r') as f:
            for i in range(3): f.readline()
            ra, dec = float(f.readline().rstrip('\n').split(' ')[-2]), float(f.readline().rstrip('\n').split(' ')[-2])
            ztf_ra_dec.add_row([ra, dec, file])

    # File combination: match ra & dec then combine data
    sens = 0.01
    overlap_atlas_names, overlap_ztf_names = [], []
    for atlas_ra, atlas_dec, atlas_path in atlas_ra_dec.iterrows():
        for ztf_ra, ztf_dec, ztf_path in ztf_ra_dec.iterrows():
            # Unique names for files
            atlas_name = atlas_path.split('/')[-1].split('.')[0]
            ztf_name = ztf_path.split('/')[-1].split('.')[0]

            # Combine files and save to ATLAS-ZTF data folder
            if (abs(atlas_ra-ztf_ra) < sens) and (abs(atlas_dec-ztf_dec) < sens):
                print(f"[+++] Combining data for {atlas_name} ({atlas_path}) and {ztf_name} ({ztf_path})...")
                # ZTF data ============================================================================================
                data = np.genfromtxt(ztf_path, delimiter=' ', dtype=str, skip_header=54)
                hdr, data = list(data[0]), data[1:]
                for i in range(len(hdr)): hdr[i] = hdr[i][:-1]

                # Make table
                ztf_table = Table()
                for h in hdr:
                    try:
                        ztf_table[h] = data[:, hdr.index(h)].astype(float)
                    except ValueError:
                        ztf_table[h] = data[:, hdr.index(h)]

                # Add RA & DEC
                ztf_table['ra'] = np.full(len(ztf_table), ztf_ra)
                ztf_table['dec'] = np.full(len(ztf_table), ztf_dec)

                # Fix time, JD to MJD
                ztf_table['jd'] = np.array(ztf_table['jd']).astype(float) - 2400000.5
                ztf_table['forcediffimflux'][ztf_table['forcediffimflux'] == 'null'] = 'nan'
                ztf_table['forcediffimfluxunc'][ztf_table['forcediffimfluxunc'] == 'null'] = 'nan'
                ztf_table['mag'] = (-2.5 * np.log10(np.array(ztf_table['forcediffimflux']).astype(float))) + ztf_table[
                    'zpdiff']
                ztf_table['dmag'] = np.abs(-1.08573620476 * (np.array(ztf_table['forcediffimfluxunc']).astype(float)
                                                             / np.array(ztf_table['forcediffimflux']).astype(float)))

                # Add parity with CSP & ATLAS
                ztf_table.remove_columns(['index', 'field', 'ccdid', 'qid', 'pid', 'infobitssci', 'sciinpseeing',
                                          'scibckgnd', 'scisigpix', 'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms',
                                          'clrcoeff', 'clrcoeffunc', 'ncalmatches', 'exptime', 'adpctdif1', 'adpctdif2',
                                          'diffmaglim', 'programid', 'rfid', 'forcediffimsnr', 'forcediffimchisq',
                                          'forcediffimfluxap', 'forcediffimfluxuncap', 'forcediffimsnrap',
                                          'aperturecorr',
                                          'dnearestrefsrc', 'nearestrefmag', 'nearestrefmagunc', 'nearestrefchi',
                                          'nearestrefsharp', 'refjdstart', 'refjdend', 'procstatu'])
                for h_old, h_new in zip(['zpdiff', 'jd', 'forcediffimflux', 'forcediffimfluxunc'],
                                        ['zp', 'mjd', 'flux', 'dflux']):
                    ztf_table[h_old].name = h_new
                # ATLAS data ==========================================================================================
                with open(atlas_path, 'r') as f:
                    hdr = f.readline()[1:-1].split(',')
                data = np.genfromtxt(atlas_path, dtype='str', delimiter=',', skip_header=1)

                # Make table
                atlas_table = Table()
                for h in hdr:
                    try:
                        atlas_table[h] = data[:, hdr.index(h)].astype(float)
                    except ValueError:
                        atlas_table[h] = data[:, hdr.index(h)]

                # Add parity with CSP & ZTF
                atlas_table.remove_columns(
                    ['err', 'chi/N', 'x', 'y', 'maj', 'min', 'phi', 'apfit', 'mag5sig', 'Sky', 'Obs'])
                atlas_table['zp'] = np.full(len(atlas_table), np.nan)
                for h_old, h_new in zip(['JD', 'm', 'dm', 'uJy', 'duJy', 'F', 'RA', 'Dec'],
                                        ['mjd', 'mag', 'dmag', 'flux', 'dflux', 'filter', 'ra', 'dec']):
                    atlas_table[h_old].name = h_new
                # Save combined file ==================================================================================
                combined_name = (atlas_name) + '|' + (ztf_name)
                if os.path.exists(f"{atlas_ztf_loc}/{combined_name}.txt") == False:
                    with open(f"{atlas_ztf_loc}/{combined_name}.txt", 'w') as f:
                        print(f"# Comined data from ATLAS '{atlas_name}' ({atlas_path}) and ZTF '{ztf_name}' ({ztf_path})",
                              file=f)
                        print_hdr = str(list(ztf_table.columns)).replace("'", '')[1:-1]
                        print(print_hdr, file=f)
                        for tb in [atlas_table, ztf_table]:
                            for i in range(len(tb)):
                                print(f"{tb[i]['filter']}, {tb[i]['zp']}, {tb[i]['mjd']}, "
                                      f"{tb[i]['flux']}, {tb[i]['dflux']}, {tb[i]['ra']}, "
                                      f"{tb[i]['dec']}, {tb[i]['mag']}, {tb[i]['mag']}", file=f)
                    print(f"      Data save to '{atlas_ztf_loc}/{combined_name}.txt'...")
                overlap_atlas_names.append(atlas_name)
                overlap_ztf_names.append(ztf_name)

    # Save copy of files that don't overlapp
    for atlas_ra, atlas_dec, atlas_path in atlas_ra_dec.iterrows():
        atlas_name = atlas_path.split('/')[-1].split('.')[0]
        if atlas_name not in overlap_atlas_names:
            print(f"[~~~] Saving 'unique' copy of {atlas_name} @ '{data_atlas_unique_loc}{atlas_name}.txt'...")
            shutil.copy(atlas_path, f"{data_atlas_unique_loc}{atlas_name}.txt")
    for ztf_ra, ztf_dec, ztf_path in ztf_ra_dec.iterrows():
        ztf_name = ztf_path.split('/')[-1].split('.')[0]
        if ztf_name not in overlap_ztf_names:
            print(f"[~~~] Saving 'unique' copy of {ztf_name} @ '{data_ztf_unique_loc}{ztf_name}.txt'...")
            shutil.copy(ztf_path, f"{data_ztf_unique_loc}{ztf_name}.txt")

    # Data Stats
    print("===========================================================================================================")
    print(f"[~~~] ATLAS files '{atlas_loc}': {len(glob.glob(atlas_loc+'/*.txt'))}")
    print(f"      'Unique' ATLAS files '{data_atlas_unique_loc}': {len(glob.glob(data_atlas_unique_loc+'*.txt'))}")
    print(f"[~~~] ZTF files '{ztf_loc}': {len(glob.glob(ztf_loc+'/*.txt'))}")
    print(f"      'Unique' ZTF files '{data_ztf_unique_loc}': {len(glob.glob(data_ztf_unique_loc+'*.txt'))}")
    print(f"[~~~] ATLAS-ZTF files '{atlas_ztf_loc}': {len(glob.glob(atlas_ztf_loc+'/*.txt'))}")

    return
def combine_snpy_salt(snpy_path: str, salt_path: str, save_loc: str = ''):
    tb_snpy = utils.default_open(snpy_path, 'True')
    tb_salt = utils.default_open(salt_path, 'True')

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
                  f"Tmax, Tmax_err, stretch, stretch_err, color, color_err, peak_mag, peak_mag_err", file=f)
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
              "peak_mag, peak_mag_err", file=f)
        for sn in sne:
            if sn.path.split('_')[-2].upper() == "SNPY":
                param_dict = {'t_max': 'Tmax', 'stretch': 'st', 'color': 'EBVhost'}
            elif sn.path.split('_')[-2].upper() == "SALT":
                param_dict = {'t_max': 't0', 'stretch': 'x1', 'color': 'c'}
            print(f"{sn.objname}, "
                  f"{sn.path.split('/')[1].split('-')[0].upper()}, "
                  f"{sn.path.split('/')[1].split('-')[1].lower()}, "
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
                  f"{np.average(sn.mag[(sn.mjd.astype(float) > sn.params[param_dict['t_max']]['value'] - 5) & (sn.mjd.astype(float) < sn.params[param_dict['t_max']]['value'] + 5)].astype(float))}, "
                  f"{np.average(sn.dmag[(sn.mjd.astype(float) > sn.params[param_dict['t_max']]['value'] - 5) & (sn.mjd.astype(float) < sn.params[param_dict['t_max']]['value'] + 5)].astype(float))}",
                  file=f)
    return
def selection_criteria(path: str, save_loc: str = ''):
    print(f"[+++] Selection criteria for '{path}'...")
    # 999 = no upper limit, -999 = no lower limit
    criteria = {
        'z':            [0.015, 999],
        'EBVhost':      [-0.2, 0.3],
        'EBVhost_err':  [-999, 0.1],
        'st':           [-999, 1.0],
        'st_err':       [-999, 0.1],
        'c':            [-0.6, 0.6],
        'c_err':        [-999, 0.1],
        'x1':           [-3.2, 3.2],
        'x1_err':       [-999, 1.0],
        'Tmax_err':     [-999, 1.0],
        't0_err':       [-999, 1.0],
        'mu_err':       [-999, 0.2],
    }

    # Load data & Seperate tables
    tb = utils.default_open(path, True)

    # Preform general cuts
    for cut in ['z', 'mu_err']:
        print(f"      Selecting on {cut}: {criteria[cut][0]} < {cut} < {criteria[cut][1]}... ", end='')
        print(f"SNe = {len(tb)} ---> ", end='')
        tb = tb[(tb[cut] >= criteria[cut][0]) &
                (tb[cut] <= criteria[cut][1])]
        print(len(tb), end='\n')

    # Preform algorithm specific cuts
    tb_snpy = tb[tb['algo'] == 'SNPY']
    tb_salt = tb[tb['algo'] == 'SALT']
    for snpy_param, salt_param, hdr_name in zip(['c', 'c_err', 'x1', 'x1_err', 't0_err'],
                                                ['EBVhost', 'EBVhost_err', 'st', 'st_err', 'Tmax_err'],
                                                ['color', 'color_err', 'stretch', 'stretch_err', 'Tmax_err']):
        # Cut on SNooPy Parameter
        print(f"      Selecting on {hdr_name}: {criteria[snpy_param][0]} < {snpy_param} < {criteria[snpy_param][1]}... ", end='')
        print(f"SNe = {len(tb_snpy) + len(tb_salt)} ---> ", end='')
        tb_snpy = tb_snpy[(tb_snpy[hdr_name] >= criteria[snpy_param][0]) &
                          (tb_snpy[hdr_name] <= criteria[snpy_param][1])]
        print(len(tb_snpy) + len(tb_salt), end='\n')

        # Cut on SALT Parameter
        print(f"      Selecting on {hdr_name}: {criteria[salt_param][0]} < {salt_param} < {criteria[salt_param][1]}... ", end='')
        print(f"SNe = {len(tb_snpy) + len(tb_salt)} ---> ", end='')
        tb_salt = tb_salt[(tb_salt[hdr_name] >= criteria[salt_param][0]) &
                          (tb_salt[hdr_name] <= criteria[salt_param][1])]
        print(len(tb_snpy) + len(tb_salt), end='\n')

    # Save new parameter file
    if len(save_loc) != 0:
        print(f"[~~~] Saving new parameter file to '{save_loc}'...")
        with open(save_loc, 'w') as f:
            print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(tb_snpy) + len(tb_salt)}", file=f)
            print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
            print(f"objname, origin, subtype, algo, ra, dec, z, z_cmb, MJDs, MJDe, mu, mu_err, hostMass, hostMass_err, "
                  f"Tmax, Tmax_err, stretch, stretch_err, color, color_err, peak_mag, peak_mag_err", file=f)
            for tb in [tb_snpy, tb_salt]:
                for row in tb:
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
                          f"{row['peak_mag']}, {row['peak_mag_err']}", file=f)
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
