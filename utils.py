import os
import requests
import numpy as np
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM

def get_apikeys(apikeys_loc: str = 'txts/api_keys.txt') -> dict:
    """
    :param apikeys_loc: Location of API keys file.
    :return: dictionary with API keys.
    """
    # Check file exsists
    if not os.path.isfile(apikeys_loc):
        raise FileNotFoundError('[!!!] API keys file not found!')

    # Read api_keys.txt
    APIKEY = {}
    with open(apikeys_loc, 'r') as f:
        # Check if file empty
        temp = f.readlines()
        if len(temp) > 0:
            for line in temp:
                line = line[:-1].split(', ')
                if len(line) != 2:
                    continue
                APIKEY.update({line[0]: line[1]})

    # If empty write user input to api_keys.txt
    if len(APIKEY) == 0:
        with open(apikeys_loc, 'w') as f:
            print('[!!!] API keys file empty! Please enter the following API keywords...')
            f.write(f"tns_bot_id, {input('tns_bot_id: ')}\n")
            f.write(f"tns_bot_name, {input('tns_bot_name: ')}\n")
            f.write(f"tns_bot_api_key, {input('tns_bot_api_key: ')}\n")
            f.write(f"atlas_key, {input('atlas_key [Optional, Press enter to skip...]: ')}\n")
    return APIKEY
def get_constants(constant_loc: str = 'txts/constants.txt') -> dict:
    """
    :param constant_loc: Location of constants.txt file.
    :return: dictionary with constants.
    """
    CONSTANTS = {}
    with open(constant_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            CONSTANTS.update({line[0]: line[1]})
    return CONSTANTS
def get_tnskey(tnskey_loc: str = 'txts/TNS_key.txt') -> Table:
    """
    :param tnskey_loc: Location of TNS_key.txt file.
    :return: An astropy Table with TNS information.
    """
    with open(tnskey_loc, 'r') as f:
        hdr = f.readline().rstrip('\n').split(', ')
        TNS_KEY = Table(names=hdr, dtype=[float, float, str, str, float])
        for line in f.readlines():
            n_line = line.rstrip('\n').split(', ')
            if len(n_line) < 5:
                continue
            TNS_KEY.add_row(n_line)
    return TNS_KEY
def check_tnskey(target_ra: float, target_dec: float, sens: float = 0.01) -> Table or None:
    """
    :param target_ra: RA to check.
    :param target_dec: DEC to check.
    :param sens: sensativity of RA and DEC detection; default = 0.01.
    :return: astropy Table of TNS key data.
    """
    TNS_KEY = get_tnskey()
    # Check RA first
    temp_table = TNS_KEY[(TNS_KEY['ra'] > target_ra - sens) & (TNS_KEY['ra'] < target_ra + sens)]
    if len(temp_table) == 1:
        pass
    elif len(temp_table) == 0:
        print(f"[~~~] Target: '{target_ra}', '{target_dec}' not found with {sens} sensitivity...")
        return None
    else:
        # If ambiguous check DEC
        temp_table = TNS_KEY[(TNS_KEY['dec'] > target_dec - sens) & (TNS_KEY['dec'] < target_dec + sens)]
        if len(temp_table) == 1:
            pass
        else:
            print(f"[~~~] Target: '{target_ra}', '{target_dec}' not found with {sens} sensitivity...")
            return None
    return temp_table
def append_tnskey(ra: float, dec: float, objname: str, z: float, discdate: float, tnskey_loc: str = 'txts/TNS_key.txt'):
    with open(tnskey_loc, 'a') as f:
        f.write(f"{ra}, {dec}, {objname}, {z}, {discdate}\n")
    return
def check_mass_key(objname: str, mode: str, hostMass: float = np.nan, hostMassErr: float = np.nan,
                   mass_key_loc: str = 'txts/mass_key.txt') -> (float, float):
    """
    :param objname: Name of SN to check for in mass key.
    :param mode: 'read' = check mass key; 'write' = append info to mass key.
    :param mass_key_loc: Location of mass_key.txt file.
    :return:
    """
    if mode == 'read':
        with open(mass_key_loc, 'r') as f:
            f.readline()  # Skip header
            for line in f.readlines():
                n_line = line.rstrip('\n').split(',')
                if n_line[0] == objname:
                    hostMass, hostMassErr = float(n_line[1]), float(n_line[2])
                    break
    elif mode == 'write':
        with open(mass_key_loc, 'a') as f:
            f.write(f"{objname},{hostMass},{hostMassErr}\n")
    return hostMass, hostMassErr
def current_cosmo(H0=70, O_m=0.3):
    return FlatLambdaCDM(H0, O_m)
def default_open(path: str, table_mode: bool = False, delimiter: str = ', '):
    """
    Open file that has the typical format of this project's 'txt' files.
    :param path: str; location of file.
    :param table_mode: bool; whether or not to return the data as an astropy table.
    :param delimiter: str; delimiter character.
    :return: (list, np.array[str]) | astropy.table.Table.
    """
    data = np.genfromtxt(path, dtype='str', delimiter=delimiter)
    hdr, data = data[0, :], data[1:, :]
    for i in range(len(hdr)):
        if hdr[i] in ['objname', 'origin', 'algo', 'subtype', 'hostMass', 'hostMass_err']: continue
        data[:, i] = data[:, i].astype(float)
    if table_mode:
        var_table = Table()
        hdr = hdr.tolist()
        for h in hdr:
            try: var_table[h] = data[:, hdr.index(h)].astype(float)
            except ValueError: var_table[h] = data[:, hdr.index(h)]
        return var_table
    else:
        return hdr.tolist(), data
def get_twomass():
    if not os.path.exists('twomass++_velocity_LH11.npy'):
        raise FileNotFoundError('[!!!] TwoMass velocities file "twomass++_velocity_LH11.npy" not found! '
                                'Please download from the follow...\n'
                                "https://drive.usercontent.google.com/download?id=1DGcWQPgmI2ZoHJm_zCqyscogmWwY7lQu&"
                                "export=download&authuser=0")
    return
def verify_downloads(combined_tarlist_path: str, source: str, subtype: str):
    """
    :param combined_tarlist_path: Target list of several sources (CSP/ATLAS/ZTF)
    :param source: Selected source to check for (CSP/ATLAS)
    :param subtype: Wether to check for 1991bg-like or Normals
    :return: None
    """
    print(f"[+++] Verifying downloaded '{source}' '{subtype}' files via {combined_tarlist_path}...")
    tarlist = np.genfromtxt(combined_tarlist_path, delimiter=',', dtype=str, skip_header=1)
    tb = Table(names=tarlist[0], data=tarlist[1:])
    num_missing = 0
    for i, row in enumerate(tb[tb['Source'] == source]):
        n = row['Name'].split(' ')[-1]

        # Path verification
        if subtype.lower() not in ['91bg', 'norm']:
            print(f"[!!!] '{subtype}' is an invalid 'subtype' selection!")
            return
        if source.lower() == 'csp':
            path = f"data/CSP-{subtype.lower()}/SN{n}_snpy.txt"
        elif source.lower() == 'atlas':
            path = f"data/ATLAS-{subtype.lower()}/ATLAS{n}.txt"
        elif source.lower() == 'ztf' and subtype.lower() == '91bg':
            path = f"data/ZTF-91bg/ZTF{n}.txt"
        else:
            print(f"[!!!] '{source}' is an invalid 'source' selection!")
            return

        # Check path
        if os.path.exists(path) == False:
            num_missing += 1
            print(f"[{num_missing}] {n} has not been downloaded!")
    if num_missing > 0:
        print(f"[~~~] Missing '{num_missing}' files, it is recommended to download them.")
    return
def verifiy_overlap():
    tarlist_91bg = np.genfromtxt('txts/target_files/tarlist_CSP-ATLAS-ZTF_91bg.csv',
                                 delimiter=',', skip_header=1, dtype=str)
    tb_91bg = Table(names=tarlist_91bg[0], data=tarlist_91bg[1:])
    tarlist_norm = np.genfromtxt('txts/target_files/tarlist_CSP-ATLAS_norm.csv',
                                 delimiter=',', skip_header=1, dtype=str)
    tb_norm = Table(names=tarlist_norm[0], data=tarlist_norm[1:])

    for tb, lb in zip([tb_91bg, tb_norm], ['1991bg-like', 'Normal']):
        print(f'[~~~] Verifying overlap for {lb} SNe...')
        tb_CSP = tb[tb['Source'] == 'CSP']
        tb_ATLAS = tb[tb['Source'] == 'ATLAS']
        tb_ZTF = tb[tb['Source'] == 'ZTF']

        for n in tb['Name']:
            if (n in list(tb_CSP['Name'])) and (n in list(tb_ATLAS['Name'])):
                print(n, "CSP+ATLAS")
            elif (n in list(tb_CSP['Name'])) and (n in list(tb_ZTF['Name'])):
                print(n, "CSP+ZTF")
            elif (n in list(tb_ATLAS['Name'])) and (n in list(tb_ZTF['Name'])):
                print(n, "ATLAS+ZTF")
            else:
                pass
    return