import time
import glob
import utils  # Import of utils.py
import shutil
import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

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
def download(tns_targets_path: str, save_loc: str):

    return


if __name__ == "__main__":
    start = time.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(time.time() - start, 4), 'seconds\n|---------------------------|')
