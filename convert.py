import os
import sys
import glob
import datetime
import numpy as np
from astropy.table import Table

CSP_CONSTANTS = {
    'csp_zpts_u':    12.986,
    'csp_zpts_g':    15.111,
    'csp_zpts_r':    14.903,
    'csp_zpts_i':    14.545,
    'csp_zpts_B':    14.328,
    'csp_zpts_V':    14.439,
    'csp_zpts_V1':   14.388,
    'csp_zpts_V0':   14.437,
    'csp_zpts_Y':    13.921,
    'csp_zpts_Jrc1': 13.836,
    'csp_zpts_Jrc2': 13.836,
    'csp_zpts_H':    13.510,
    'csp_zpts_Ydw':  13.770,
    'csp_zpts_J':    13.866,
    'csp_zpts_H':    13.502,
    'csp_zpts_Ks':   11.968
}
CURRENTDATE = datetime.datetime.now()

def covert_file_format(file_loc, origin, save_loc):
    if origin == 'CSP':
        var_table = Table(names=['filter', 'zp', 'z', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux'],
                          dtype=[str, float, float, float, float, float, float, float, float, float])
        with open(file_loc, 'r') as f:
            csp_objname, csp_z, csp_ra, csp_dec = f.readline()[2:-1].split(' ')
            for l in f.readlines():
                l = l.split(' ')
                # Filter line
                if len(l) == 2:
                    csp_filter = str(l[1][:-1])
                    csp_zp = float(CSP_CONSTANTS['csp_zpts_' + csp_filter])
                else:
                    csp_mjd, csp_mag, csp_dmag = float(l[-3]) + 53000, float(l[-2]), float(l[-1])
                    csp_flux = 10 ** ((csp_mag - csp_zp) / -2.5)
                    csp_dflux = np.abs(csp_flux) * np.log(10) * ((1 / 2.5) * csp_dmag)

                    if csp_filter not in ['u', 'Y', 'J', 'H', 'Jrc2', 'Ydw']: # Unwanted filters
                        var_table.add_row([csp_filter, csp_zp, csp_z, csp_ra, csp_dec, csp_mjd,
                                           csp_mag, csp_dmag, csp_flux, csp_dflux])
    elif origin == 'ATLAS':
        ###
        # needs redshift
        ###
        # Load data
        with open(file_loc, 'r') as f: hdr = f.readline()[1:-1].split(',')
        data = np.genfromtxt(file_loc, dtype='str', delimiter=',', skip_header=1)

        # Check if empty
        if len(data) == 0:
            print("[!!!] Empty file...")
            return False

        # Make table
        var_table = Table()
        for h in hdr:
            try:
                var_table[h] = data[:, hdr.index(h)].astype(float)
            except ValueError:
                var_table[h] = data[:, hdr.index(h)]

        # Update coloumns
        var_table.remove_columns(['err', 'x', 'y', 'maj', 'min', 'phi', 'apfit', 'Sky', 'Obs'])
        for h_old, h_new in zip(['JD', 'm', 'dm', 'uJy', 'duJy', 'F', 'RA', 'Dec'],
                                ['mjd', 'mag', 'dmag', 'flux', 'dflux', 'filter', 'ra', 'dec']):
            var_table[h_old].name = h_new
        var_table['zp'] = np.full(var_table['mag'].shape, 23.9) # Recreate Zero-Points using ZP = m + 2.5log_10(flux)

        # # Filter on 5-sigma magnitudes (mag5sig)
        # var_table = var_table[var_table['mag'] < var_table['mag5sig']]
        #
        # # Filter out chi/N > 2
        # var_table = var_table[var_table['chi/N'] < 3.0]

        # Filter out unresonable mags
        var_table = var_table[var_table['mag'] > 0.0]

        # Cut errors > 1mag
        var_table = var_table[var_table['dmag'] < 1.0]
    elif origin == 'ZTF':
        ###
        # needs redshift
        ###
        # Load hdr & data
        with open(file_loc, 'r') as f:
            for i in range(3): f.readline()
            ztf_ra = float(f.readline().split(' ')[-2])
            ztf_dec = float(f.readline().split(' ')[-2])
        data = np.genfromtxt(file_loc, delimiter=' ', dtype=str, skip_header=54)
        hdr, data = list(data[0]), data[1:]
        for i in range(len(hdr)): hdr[i] = hdr[i][:-1]

        # Check if empty
        if len(data) == 0:
            print("[!!!] Empty file...")
            return False

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

        # # Add temp redshift
        # var_table['z'] = np.full(len(var_table), np.nan)

        var_table['jd'] = np.array(var_table['jd']).astype(float) - 2400000.5 # JD to MJD
        var_table['forcediffimflux'][var_table['forcediffimflux'] == 'null'] = 'nan' # Replace 'null' with 'nan'
        var_table['forcediffimfluxunc'][var_table['forcediffimfluxunc'] == 'null'] = 'nan' # Replace 'null' with 'nan'
        var_table = var_table[var_table['forcediffimflux'].astype(float) > 0] # Remove negative values
        var_table['mag'] = (-2.5 * np.log10(np.array(var_table['forcediffimflux']).astype(float))) + var_table['zpdiff']
        var_table['dmag'] = np.abs(-1.08573620476 * (np.array(var_table['forcediffimfluxunc']).astype(float)
                                                     / np.array(var_table['forcediffimflux']).astype(float)))

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

        # # Filter on by magnitude limit (diffmaglim)
        # var_table = var_table[var_table['mag'] < var_table['diffmaglim']]

        # # Filter out by reduced chi2 > 3.0 -- Might remove too many filters
        # var_table = var_table[var_table['forcediffimchisq'] < 3.0]

        # Filter out unresonable mags
        var_table = var_table[var_table['mag'] > 0.0]

        # Cut errors > 1mag
        var_table = var_table[var_table['dmag'] < 1.0]

    elif origin == 'HSF':
        # Open data
        data = np.genfromtxt(file_loc, delimiter=', ', dtype=str, skip_header=1)
        hdr, data = list(data[0]), data[1:]

        # Check if empty
        if len(data) == 0:
            print("[!!!] Empty file...")
            return False

        # Make table
        var_table = Table()
        for h in hdr:
            try:
                var_table[h] = data[:, hdr.index(h)].astype(float)
            except ValueError:
                var_table[h] = data[:, hdr.index(h)]

        # Cut errors > 1mag
        var_table = var_table[var_table['dmag'] < 1.0]

    # # Check for filter >= 2
    # if len(np.unique(var_table['filter'])) < 2:
    #     return False

    # Check if less than 2 points
    if len(var_table) <= 2:
        # print("[!!!] Empty file...")
        return False

    # Save data to new location
    with open(save_loc, 'w') as f:
        print('ra,dec,mjd,mag,dmag,flux,dflux,zp,filter', file=f)
        for row in var_table:
            print(f"{row['ra']},"
                  f"{row['dec']},"
                  f"{row['mjd']},"
                  f"{row['mag']},"
                  f"{row['dmag']},"
                  f"{row['flux']},"
                  f"{row['dflux']},"
                  f"{row['zp']},"
                  f"{row['filter']}",
                  file=f)

    return True
def combine_atlas_ztf(subtype: str = '91bg'):

    # Get lists of ATLAS & ZTF names
    atlas_sne, ztf_sne = [], []
    for file in glob.glob(f"data/{subtype.upper()}_*.txt"):
        if file.split('_')[-2] == 'ATLAS':
            atlas_sne.append(file.split('_')[-1][:-4])
        elif file.split('_')[-2] == 'ZTF':
            ztf_sne.append(file.split('_')[-1][:-4])

    # Check overlapp
    for i, n in enumerate(atlas_sne):
        if n in ztf_sne:
            # print(f"[{i+1} / {len(atlas_sne)}] Combining ATLAS & ZTF for {n}...")
            # Open data
            atlas_tb = Table(names=["ra","dec","mjd","mag","dmag","flux","dflux","zp","filter"],
                             data=np.genfromtxt(f"data/91BG_ATLAS_{n}.txt", delimiter=',', skip_header=1, dtype=str),
                             dtype=[float, float, float, float, float, float, float, float, str])
            ztf_tb = Table(names=["ra","dec","mjd","mag","dmag","flux","dflux","zp","filter"],
                           data=np.genfromtxt(f"data/91BG_ZTF_{n}.txt", delimiter=',', skip_header=1, dtype=str),
                           dtype=[float, float, float, float, float, float, float, float, str])

            # Save data to new location
            with open(f"data/{subtype.upper()}_ATLAS-ZTF_{n}.txt", 'w') as f:
                print('ra,dec,mjd,mag,dmag,flux,dflux,zp,filter', file=f)
                for row in atlas_tb:
                    print(f"{row['ra']},"
                          f"{row['dec']},"
                          f"{row['mjd']},"
                          f"{row['mag']},"
                          f"{row['dmag']},"
                          f"{row['flux']},"
                          f"{row['dflux']},"
                          f"{row['zp']},"
                          f"{row['filter']}",
                          file=f)
                for row in ztf_tb:
                    print(f"{row['ra']},"
                          f"{row['dec']},"
                          f"{row['mjd']},"
                          f"{row['mag']},"
                          f"{row['dmag']},"
                          f"{row['flux']},"
                          f"{row['dflux']},"
                          f"{row['zp']},"
                          f"{row['filter']}",
                          file=f)

            # Remove original files
            os.remove(f"data/{subtype.upper()}_ATLAS_{n}.txt")
            os.remove(f"data/{subtype.upper()}_ZTF_{n}.txt")
    return
def add_z_and_discdate(subtype: str = '91bg'):
    if subtype == '91bg':
        tarlist = np.genfromtxt('txts/target_files/new_sn1991bglike_tarlist.csv', delimiter=',', dtype=str, skip_header=2)
    elif subtype == 'norm':
        tarlist = np.genfromtxt('txts/target_files/new_normal_tarlist.csv', delimiter=',', dtype=str, skip_header=2)

    tarlist_hdr = "Name,RA,DEC,Redshift,Redshift Error,Discovery Date (UT),Data Source(s)".split(',')

    files = glob.glob(f'data/{subtype.upper()}*.txt')
    for i, file in enumerate(files):
        glob_name = f"SN {file.split('_')[-1][:-4]}"

        if glob_name not in list(tarlist[:, 0]):
            os.remove(file)
            continue

        row = tarlist[tarlist[:, 0] == glob_name][0]
        z, z_err, discdate = row[3], row[4], row[5]

        with open(file, 'r') as f: hdr = f.readline().rstrip('\n').split(',')
        tb = Table(names=hdr, data=np.genfromtxt(file, delimiter=',', skip_header=1, dtype=str))

        # Save data to new location
        with open(file, 'w') as f:
            print('ra,dec,z,z_err,discdate,mjd,mag,dmag,flux,dflux,zp,filter', file=f) # , file=f
            for r in tb:
                print(f"{r['ra']},"
                      f"{r['dec']},"
                      f"{z},"
                      f"{z_err},"
                      f"{discdate},"
                      f"{r['mjd']},"
                      f"{r['mag']},"
                      f"{r['dmag']},"
                      f"{r['flux']},"
                      f"{r['dflux']},"
                      f"{r['zp']},"
                      f"{r['filter']}",
                      file=f)
    return
def data_stats(subtype: str = '91bg'):
    csp, atlas, ztf, atlasztf = [], [], [], []
    for file in glob.glob(f"data/{subtype.upper()}_*.txt"):
        if file.split('_')[-2] == "CSP": csp.append(file)
        elif file.split('_')[-2] == "ATLAS": atlas.append(file)
        elif file.split('_')[-2] == "ZTF": ztf.append(file)
        elif file.split('_')[-2] == "ATLAS-ZTF": atlasztf.append(file)
    print("======================================================")
    print(f"{subtype.upper()} Statistics...")
    print("======================================================")
    print(f"CSP: #{len(csp)}")
    print(f"ATLAS: #{len(atlas)}")
    print(f"ZTF: #{len(ztf)}")
    print(f"ATLAS-ZTF: #{len(atlasztf)}")
    print("======================================================")

    return
def convert_dust(data_loc: str = 'txts/global_dust.txt', save_loc: str = 'txts/global_dust_params.txt'):
    # Clean dust files
    data = np.genfromtxt(data_loc, delimiter=' ', skip_header=1, dtype=str)
    with open(save_loc, 'w') as f:
        print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(data[:, 0])}", file=f)
        print(f"# Avg. 'av_16': {data[:, 1].astype(float).mean():.4f}", file=f)
        print(f"# Avg. 'av_50': {data[:, 2].astype(float).mean():.4f}", file=f)
        print(f"# Avg. 'av_84': {data[:, 3].astype(float).mean():.4f}", file=f)
        print("objname, av_16, av_50, av_84", file=f)
        for n in range(len(data[:, 0])):
            print(f"{data[n, 0]}, {data[n, 1]}, {data[n, 2]}, {data[n, 3]}", file=f)
    return
def dust_check():
    # Needed dust files
    local_dust = utils.default_open('txts/local_dust_params.txt', True)
    global_dust = utils.default_open('txts/global_dust_params.txt', True)
    norms_cut = utils.default_open('results/norm_snpy-salt_params_cut.txt', True)
    sne_91bg_cut = utils.default_open('results/91bg_snpy-salt_params_cut.txt', True)

    all_names = list(norms_cut['objname']) + list(sne_91bg_cut['objname'])
    with open('needed_dust.txt', 'w') as f:
        print("objname,missing", file=f)
        for n in all_names:
            if (n not in list(local_dust['objname'])) & (n not in list(global_dust['objname'])):
                print(f"{n},local+global", file=f)
            elif (n not in list(local_dust['objname'])):
                print(f"{n},local", file=f)
            elif (n not in list(global_dust['objname'])):
                print(f"{n},global", file=f)
    return
def run():

    # tar_list = np.genfromtxt('../txts/target_files/sn1991bglike_tarlist.csv', delimiter=',', dtype=str, skip_header=2)
    # with open('../txts/target_files/new_sn1991bglike_tarlist.csv', 'w') as f:
    #     print(f"# Object Type: 1991bg-like SNe Ia; Date: 2005 to 2025\n"
    #           f"Name,RA,DEC,Redshift,Redshift Error,Discovery Date (UT),Data Source(s)", file=f)
    #     for n in range(len(tar_list)):
    #         print(f"{tar_list[n, 0]},"
    #               f"{tar_list[n, 1]},"
    #               f"{tar_list[n, 2]},"
    #               f"{tar_list[n, 3]},"
    #               f"None,"
    #               f"{tar_list[n, 4]},"
    #               f"{tar_list[n, 5]}", file=f)




    # # Convert CSP norms
    # csp_files = glob.glob('../data/CSP-norm/*')
    # n_csp = 0
    # for file in csp_files:
    #     s = covert_file_format(file, 'CSP', f"data/NORM_CSP_{file.split('/')[-1][3:]}")
    #     if s: n_csp += 1
    # print(f"{n_csp}/{len(csp_files)}... CSP")
    #
    # # Convert HSF files
    # hsf_files = glob.glob('../data/combined_norm/2*')
    # n_hsf = 0
    # for file in hsf_files:
    #     s = covert_file_format(file, 'HSF', f"data/NORM_ATLAS-ZTF_{file.split('/')[-1]}")
    #     if s: n_hsf += 1
    # print(f"{n_hsf}/{len(hsf_files)}... HSF")
    #
    # # Combine ATLAS-ZTF
    # combine_atlas_ztf('norm')
    #
    # # Add in redshift and discovery date
    # add_z_and_discdate('norm')
    # #
    # # # Stats
    # data_stats('norm')

    # =================================================================================================================

    # Convert CSP 1991bgs
    csp_files = glob.glob('old_data/CSP-91bg/*')
    n_csp = 0
    for file in csp_files:
        s = covert_file_format(file, 'CSP', f"data/91BG_CSP_{file.split('/')[-1][3:]}")
        if s: n_csp += 1
    print(f"{n_csp}/{len(csp_files)}... CSP")

    # Convert ATLAS 1991bgs
    atlas_files = glob.glob('old_data/ATLAS-91bg/*')
    n_atlas = 0
    for file in atlas_files:
        s = covert_file_format(file, 'ATLAS', f"data/91BG_ATLAS_{file.split('/')[-1][5:]}")
        if s: n_atlas += 1
    print(f"{n_atlas}/{len(atlas_files)}... ATLAS")

    # Convert ZTF 1991bgs
    ztf_files = glob.glob('old_data/ZTF-91bg/*')
    n_ztf = 0
    for file in ztf_files:
        s = covert_file_format(file, 'ZTF', f"data/91BG_ZTF_{file.split('/')[-1][3:]}")
        if s: n_ztf += 1
    print(f"{n_ztf}/{len(ztf_files)}... ZTF")

    # Combine ATLAS-ZTF
    combine_atlas_ztf('91bg')

    # Add in redshift and discovery date
    add_z_and_discdate('91bg')

    # Stats
    data_stats('91bg')

    return
