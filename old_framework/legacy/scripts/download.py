import io
import re
import sys
import time
import glob
import requests
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from old_file_system.scripts.general import get_APIkeys

API_KEYS = get_APIkeys('../api_keys.txt')
BASEURL = "https://fallingstar-data.com/forcedphot"
TOKEN = API_KEYS['atlas_key']

def retrive_objname_ra_dec(path: str = "../txts/norm_Ia_tns_24-25.csv"):
    names, ra, dec = [], [], []
    for i, line in enumerate(open(path)):
        if i == 0:
            continue  # Skip header
        n_line = line.split('",')
        n_ra, n_dec = n_line[2][1:], n_line[3][1:]
        c = SkyCoord(n_ra, n_dec, frame='icrs', unit=(u.hourangle, u.deg))
        ra.append(c.ra.deg)
        dec.append(c.dec.deg)
        names.append(n_line[1][4:])
    return names, ra, dec
def initate_download(ra: str, dec: str, headers: dict[str, str]):
    task_url = None
    while not task_url:
        with requests.Session() as s:
            resp = s.post(f"{BASEURL}/queue/", headers=headers, data={'ra': ra, 'dec': dec})

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
        textdata = s.get(result_url, headers=headers).text
        dfresult = pd.read_csv(io.StringIO(textdata.replace("###", "")), sep='\s+')
        dfresult.to_csv(save_path, index=False)
        print('Saved data to... ', save_path)
    return
def main():
    objnames, ra, dec = retrive_objname_ra_dec("../norm_Ia_tns_24-25.csv")
    save_loc = "../data/ATLASnorms/"  # Output directory for light curves
    headers = {'Authorization': f'Token {TOKEN}', 'Accept': 'application/json'}

    # Get known names
    known_names = []
    for k in glob.glob(save_loc + "*.txt"):
        known_names.append(k.split('/')[-1][5:-4])

    # Download light curves for each set of coordinates
    for i, n in enumerate(zip(objnames, ra, dec)):
        hdr = f"[{i+1} / {len(objnames)}] Downloading {n[0]} ({round(n[1], 3)}, {round(n[2],3)})... ==================="
        csr = f"{'='*len(hdr)}"
        #######################
        print(hdr)

        # Check if already downloaded
        if n[0] in known_names:
            print(f'[---] Already downloaded! Skipping...\n{csr}')
            continue

        # Attempt download
        task_url = initate_download(n[1], n[2], headers)
        result_url = check_download(task_url, headers)
        if result_url is None:
            print(f"[!!!] Result unavailable! Skipping...\n{csr}")
            continue

        # Successful Download
        print("[+++] ", end='')
        save_download(result_url, headers, f'{save_loc}ATLAS{n[0]}.txt')
        print(csr)


if __name__ == "__main__":
    start = time.time()  # Runtime tracker
    main()
    print('|---------------------------|\n Run-time: ', round(time.time() - start, 4), 'seconds\n|---------------------------|')
