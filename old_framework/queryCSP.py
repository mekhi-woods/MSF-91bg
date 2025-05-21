import os
import glob
import time
import utils  # Import of utils.py
import shutil
import tarfile
import requests
import numpy as np
from astropy.table import Table

def download(save_loc: str):
    """
    :param save_loc: Location to save tarball and 'DR3/' folder with all CSP-I Third Release data.
    :return: None
    """
    url = "https://csp.obs.carnegiescience.edu/data/CSP_Photometry_DR3.tgz"
    tar_name = "allCSP-I_tarball.tgz"
    path = save_loc+tar_name

    # Check if already downloded
    if os.path.exists(save_loc+tar_name):
        print("[+++] File previous downloaded! Skipping...")
    else:
        # Use request to get file
        response = requests.get(url)

        # Save file
        if response.status_code == 200:
            with open(path, "wb") as f:
                f.write(response.content)
            print(f"[+++] File '{tar_name}' downloaded successfully to {path}.")
        else:
            raise(f"[!!!] Failed to download the file. Status code: {response.status_code}")

    # Unpack tarball
    extract_path = save_loc
    if os.path.exists(extract_path+'DR3/'):
        print("[+++] File previous extracted! Skipping...")
    else:
        print(f"[+++] Extracting tarball to {extract_path}...")
        with tarfile.open(path, "r:gz") as tar:
            tar.extractall(extract_path)
    return
def separate_data(subtype: str, data_loc: str, save_loc: str):
    """
    :param tar_list: List of target names to seperate
    :param data_loc:
    :param save_loc:
    :return:
    """
    targets = np.genfromtxt('txts/DR3_subtypes.csv', delimiter=',', skip_header=1, dtype=str)
    hdr, targets = targets[0], targets[1:]
    tb_targets = Table(names=hdr, data=targets)
    target_names = tb_targets[tb_targets['Subtype'] == subtype]['Name'].tolist()

    # Sort through data
    print(f"[+++] Relocating '{subtype}' SNe from '{data_loc}' to '{save_loc}'...")
    for path in glob.glob(data_loc+'*.txt'):
        name = 'CSP' + path.split('/')[-1][2:-9]
        if f"SN {name[3:]}" in target_names:
            shutil.copy(path, save_loc+name+'.txt')
    return


if __name__ == "__main__":
    start = time.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(time.time() - start, 4), 'seconds\n|---------------------------|')
