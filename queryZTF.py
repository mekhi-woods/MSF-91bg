import os
import wget
import time
import glob
import utils  # Import of utils.py
import shutil
import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time as astrotime

import requests
import time

ZTF_API_URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?" # Light Curve URL
ZTF_USERNAME = utils.get_apikeys()["ztf_user"]
ZTF_PASSWORD = utils.get_apikeys()["ztf_pass"]

# File Modification ===================================================================================================
def verify():
    tarlist = np.genfromtxt("txts/target_files/tarlist_CSP-ATLAS-ZTF_91bg.csv",
                            delimiter=',', dtype=str, skip_header=1)
    tb = Table(names=tarlist[0], data=tarlist[1:])
    tb_ztf = tb[tb['Source'] == 'ZTF']

    # Check RA & DEC with target list
    overlap = 0
    missing_sne = []
    for i, row in enumerate(tb_ztf):
        c = SkyCoord(row['RA'], row['DEC'], frame='icrs', unit=(u.hourangle, u.deg))
        tarlist_ra, tarlist_dec = c.ra.deg, c.dec.deg

        target_found = False
        for file in glob.glob('data/ZTF-91bg/*.txt'):
            with open(file, 'r') as f:
                # Get ZTF-name, RA, DEC
                og_name = file.split('/')[-1].split('_')[-2]
                for i in range(3):
                    f.readline()
                file_ra, file_dec = float(f.readline().split(' ')[-2]), float(f.readline().split(' ')[-2])

            if np.abs(file_ra - tarlist_ra) < 0.01 and np.abs(file_dec - tarlist_dec) < 0.01:
                overlap += 1
                target_found = True
                print(f"[{overlap}]", row['Name'], tarlist_ra, tarlist_dec, '--->', og_name, file_ra, file_dec)
                shutil.copy(file, f"data/new-ZTF-91bg/ZTF{row['Name'][3:]}.txt")
        if target_found == False:
            missing_sne.append(row['Name'])
    if len(missing_sne) > 0:
        print(f"[~~~] The following SNe are missing ({len(missing_sne)})...\n{missing_sne}")
    return
def get_TNS_names():
    from fitter import get_TNSDetails
    import os
    for i, file in enumerate(glob.glob('data/ZTF-91bg/*.txt')):
        original_name = file.split('/')[-1].split('_')[-2]
        with open(file, 'r') as f:
            for k in range(3):
                f.readline()
            ra = float(f.readline().split(' ')[-2])
            dec = float(f.readline().split(' ')[-2])
        if os.path.exists(f'data/new-ZTF-91bg/{file.split("/")[-1]}') == False:
            print(f"[{i}] {original_name} [{round(ra, 3)}, {round(dec, 3)}] --- > ", end='')
            details = get_TNSDetails(ra, dec)
            objname = details['objname']
            print(f"{objname} ({file} ---> {f'data/new-ZTF-91bg/ZTF{objname}.txt'})")
            shutil.copy(file, f"data/new-ZTF-91bg/ZTF{objname}.txt")
            shutil.copy(file, f'data/new-ZTF-91bg/{file.split("/")[-1]}')
        # break
    return

# Querying ============================================================================================================
def get_targets(path: str):
    names, ra, dec, disc_date = [], [], [], []
    tarlist = np.genfromtxt(path, delimiter=',', dtype=str, skip_header=1)
    tb = Table(names=tarlist[0], data=tarlist[1:])
    for i, row in enumerate(tb):
        names.append(row['Name'].split('SN ')[-1])
        c = SkyCoord(row['RA'], row['DEC'], frame='icrs', unit=(u.hourangle, u.deg))
        ra.append(c.ra.deg)
        dec.append(c.dec.deg)
        if row['Discovery Date (UT)'] == 'Unknown':
            disc_date.append(None)
        else:
            disc_date.append(row['Discovery Date (UT)'])
    return names, ra, dec, disc_date
def download(tar_list: str, save_loc: str):
    objnames, ra, dec, disc_date = get_targets(tar_list)

    # Get known names
    known_names = []
    # for k in glob.glob(save_loc + "*.csv"):
    #     known_names.append(k.split('/')[-1].split('ZTF')[-1].split('.csv')[0])

    # Download light curves for each set of coordinates
    failed_names = []
    for i, n in enumerate(zip(objnames, ra, dec, disc_date)):
        hdr = f"[{i+1} / {len(objnames)}] Downloading {n[0]} ({round(n[1], 3)}, {round(n[2],3)})... ==================="
        csr = f"{'='*len(hdr)}"
        #######################
        print(hdr)

        # Check if already downloaded
        if n[0] in known_names:
            print(f'[---] Already downloaded! Skipping...\n{csr}')
            continue

        # Construct data for request call
        raidus = 0.0028  # in deg
        t = astrotime(n[3], format='iso', scale='utc').mjd
        # post_url = f"{ZTF_API_URL}POS=CIRCLE+{round(float(n[1]), 4)}+{round(float(n[2]), 4)}+{0.0028}&NOBS_MIN=3&TIME={t - 100}+{t + 100}&FORMAT=CSV"
        post_url = f"{ZTF_API_URL}POS=CIRCLE+{float(n[1])}+{float(n[2])}+{0.01}&TIME={t - 100}+{t + 100}&FORMAT=CSV"
        print(post_url)

        try:
            # Commit wget call to temp file
            if os.path.exists(f"{save_loc}temp.csv"):
                os.remove(f"{save_loc}temp.csv")
            ztf_save_loc = f"{save_loc}ZTF{n[0]}.csv"
            wget.download(url=post_url, out=f"{save_loc}temp.csv")

            # Verify if there is data in file
            temp_data = np.genfromtxt(f"{save_loc}temp.csv", delimiter=',', dtype=str)
            if len(temp_data.shape) > 1 and len(temp_data[:, 0]) > 3:
                # Test if date range makes sense
                stripped_year = int(n[0][:4])
                utc_time = astrotime(float(temp_data[1, 3]) + 2400000.5, format='jd', scale='utc').to_datetime()
                download_year = int(str(utc_time)[:4])
                if download_year != stripped_year:
                    print(f"[!!!] Invalid time space! Downloaded file is in '{download_year}' & actual year is '{stripped_year}'")
                    failed_names.append(n[0])
                else:
                    shutil.copy(f"{save_loc}temp.csv", ztf_save_loc)
                    print(f"[+++] Success! {n[0]} saved to '{ztf_save_loc}'...\n{csr}")
            else:
                print(f"[!!!] Downloaded file is empty! \n{csr}")
                failed_names.append(n[0])

        except Exception as e:
            print(f"[!!!] Error downloading data: \n\t{e}\n{csr}")
            failed_names.append(n[0])

    # Remove temp file
    if os.path.exists(f"{save_loc}temp.csv"):
        os.remove(f"{save_loc}temp.csv")

    # Report failed downloads
    print(f"[~~~] The following resulting in ({len(failed_names)}) errors...\n[", end='')
    for n in failed_names:
        print(f"{n}", end=', ')
    print(']')

    return


if __name__ == "__main__":
    start = time.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(time.time() - start, 4), 'seconds\n|---------------------------|')
