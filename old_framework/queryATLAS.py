import io
import os
import re
import sys
import time
import glob
import utils  # Import of utils.py
import requests
import numpy as np
import pandas as pd
from datetime import datetime
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time as astrotime

# Functions ===========================================================================================================
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
def initate_download(ra: str, dec: str, disc_date: str, headers: dict[str, str], tol: int = 100):
    # Construct data for call
    if disc_date is not None:
        t = astrotime(disc_date, format='iso', scale='utc').mjd
        post_data = {'ra': ra, 'dec': dec, 'mjd_min': str(t-tol), 'mjd_max': str(t+tol)}
    else:
        post_data = {'ra': ra, 'dec': dec}

    task_url = None
    while not task_url:
        with requests.Session() as s:
            resp = s.post(f"https://fallingstar-data.com/forcedphot/queue/",
                          headers=headers, data=post_data)

            if resp.status_code == 201:  # successfully queued
                task_url = resp.json()['url']
                print(f'The task URL is {task_url}')
            elif resp.status_code == 429:  # throttled
                message = resp.json()["detail"]
                print(f'{resp.status_code} {message}')
                t_sec = re.findall(r'available in (\d+) seconds', message)
                t_min = re.findall(r'available in (\d+) minutes', message)
                if t_sec:
                    waittime = int(t_sec[0])
                elif t_min:
                    waittime = int(t_min[0]) * 60
                else:
                    waittime = 10
                print(f'Waiting {waittime} seconds')
                time.sleep(waittime)
            else:
                print(f'ERROR {resp.status_code}')
                print(resp.json())
                sys.exit()
    return task_url
def check_download(task_url: str, headers: dict[str, str]):
    result_url = None
    while not result_url:
        with requests.Session() as s:
            resp = s.get(task_url, headers=headers)

            if resp.status_code == 200:  # HTTP OK
                if resp.json()['finishtimestamp']:
                    result_url = resp.json()['result_url']
                    print(f"Task is complete with results available at {result_url}")
                    break
                elif resp.json()['starttimestamp']:
                    print(f"Task is running (started at {resp.json()['starttimestamp']})")
                else:
                    print("Waiting for job to start. Checking again in 10 seconds...")
                time.sleep(10)
            else:
                print(f'ERROR {resp.status_code}')
                print(resp.json())
                sys.exit()
    return result_url
def save_download(result_url: str, headers: dict[str, str], save_path: str):
    with requests.Session() as s:
        # Download data from url
        textdata = s.get(result_url, headers=headers).text
        dfresult = pd.read_csv(io.StringIO(textdata.replace("###", "")), sep='\s+')

        # Test if date range makes sense
        stripped_year = int(save_path.split('/')[-1][:-4].replace('ATLAS', '')[:4])
        utc_time = astrotime(np.average(dfresult['MJD']) + 2400000.5, format='jd', scale='utc').to_datetime()
        download_year = int(str(utc_time)[:4])
        if download_year != stripped_year:
            print(f"[!!!] Invalid time space! Downloaded file is in '{download_year}' & actual year is '{stripped_year}'")
            return

        # Download data
        dfresult.to_csv(save_path, index=False)
        print('[+++] Saved data to... ', save_path)
    return

# Main Call ===========================================================================================================
def download(tar_list: str, save_loc: str):
    # Get TNS CSV of RA/DEC Information
    objnames, ra, dec, disc_date = get_targets(tar_list)

    # Create TNS Header
    headers = {'Authorization': f'Token {utils.get_apikeys()["atlas_key"]}', 'Accept': 'application/json'}

    # Get known names
    known_names = []
    for k in glob.glob(save_loc + "*.txt"):
        known_names.append(k.split('/')[-1].split('ATLAS')[-1].split('.txt')[0])

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

        # # Known issue target
        # known_not_in_atlas = ("SN 2009F SN 2008bt SN 2008bi SN 2008bd SN 2007N SN 2007ba SN 2007ax SN 2007al SN 2006mr "
        #                       "SN 2006gt SN 2006bd SN 2005ke SN 2005bl SN 2025nn SN 2023sps SN 2023ex SN 2022zvu "
        #                       "SN 2021zqs SN 2021afur SN 2020yo SN 2020vae SN 2020ecn SN 2018lph SN 2018baz SN 2018ast "
        #                       "SN 2017fzw SN 2016iuh SN 2016ije SN 2016brx").replace("SN ", "").split(" ")
        # if n[0] in known_not_in_atlas:
        #     print(f'[---] Known to not be in ATLAS! Skipping...\n{csr}')
        #     continue

        # Attempt download
        task_url = initate_download(n[1], n[2], n[3], headers)
        result_url = check_download(task_url, headers)
        if result_url is not None:
            # Successful Download
            print("[+++] ", end='')
            save_download(result_url, headers, f'{save_loc}ATLAS{n[0]}.txt')
            print(csr)
        else:
            print(f"[!!!] Result unavailable! Skipping...\n{csr}")
            failed_names.append(n[0])
            continue

    # Report failed downloads
    print(f"[~~~] The following resulting in errors...\n[", end='')
    for n in failed_names:
        print(f"{n}", end=', ')
    print(']')
    return


if __name__ == "__main__":
    start = time.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(time.time() - start, 4), 'seconds\n|---------------------------|')
