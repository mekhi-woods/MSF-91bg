import os
import glob
import json
import time
import requests
import numpy as np
from astropy.table import Table
from collections import OrderedDict
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
def query_tns(coords: list or None = None, objname: str or None = None, search_radius: str = '2') -> dict:
    """
    :param coords: [Right Ascension in degrees, Declination in degrees] or None when not in use
    :param objname: Name of object or None when not in use
    :param search_radius: Radius to search out for in arcseconds
    :return: dict of TNS details
    """
    # Set up keys and default arguments
    APIKEY = get_apikeys()
    tns_bot_id, tns_bot_name, tns_bot_api_key = APIKEY['tns_bot_id'], APIKEY['tns_bot_name'], APIKEY[
        'tns_bot_api_key']
    tns_marker = (
        f'tns_marker{{"tns_id": "{int(tns_bot_id)}",'
        f'"type": "bot", "name": "{tns_bot_name}"}}'
    )
    headers = {"User-Agent": tns_marker}

    # Check mode
    ## if coords == None, then use objname
    if type(objname) == str:
        # search_obj = [
        #     ("ra", ""),
        #     ("dec", ""),
        #     ("radius", str(search_radius)),
        #     ("units", "arcsec"),
        #     ("objname", str(objname)),
        #     ("objname_exact_match", 1),
        #     ("internal_name", ""),
        #     ("internal_name_exact_match ", 0),
        #     ("objid", ""),
        #     ("public_timestamp", ""),
        # ]
        search_obj = [
            ("radius", str(search_radius)),
            ("units", "arcsec"),
            ("objname", str(objname)),
            ("objname_exact_match", 1),
        ]
    elif type(coords) == list:
        # search_obj = [
        #     ("ra", str(coords[0])),
        #     ("dec", str(coords[1])),
        #     ("radius", str(search_radius)),
        #     ("units", "arcsec"),
        #     ("objname", ""),
        #     ("objname_exact_match", 0),
        #     ("internal_name", ""),
        #     ("internal_name_exact_match ", 0),
        #     ("objid", ""),
        #     ("public_timestamp", ""),
        # ]
        search_obj = [
            ("ra", str(coords[0])),
            ("dec", str(coords[1])),
            ("radius", str(search_radius)),
            ("units", "arcsec"),
        ]
    else:
        print("[!!!] Invalid selection! Enter objname as str {'objname'} or coordinates as list {['ra', 'dec']}...")
        return {}

    # Query TNS
    search_data = {"api_key": tns_bot_api_key, "data": json.dumps(OrderedDict(search_obj))}
    response = requests.post("https://www.wis-tns.org/api/get/search", headers=headers, data=search_data)
    response = json.loads(response.text)
    transients = response["data"]
    if len(transients) > 1:
        print("[~~~] WARNING: More than one transient object found! Returning first one...")
    if transients == 'Too Many Requests':
        print('\nToo Many Requests... pausing for 30 seconds: ', end='')
        for i in range(10):
            print(' - ', end='')
            # time.sleep(5)
        print('\n')
        transients = response["data"]
    get_obj = [
        ("objname", transients[0]["objname"]),
        ("objid", transients[0]["objid"]),
        ("photometry", "0"),
        ("spectra", "0"),
    ]
    get_data = {"api_key": tns_bot_api_key, "data": json.dumps(OrderedDict(get_obj))}
    response = requests.post("https://www.wis-tns.org/api/get/object", headers=headers, data=get_data)
    response = json.loads(response.text)
    details = response["data"]
    return details
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

    # Check if file empty
    if len(data.shape) == 1:
        hdr, data = data, np.array([])
        if table_mode:
            return Table(names=hdr, dtype=[str]*len(hdr))
        else:
            return hdr.tolist(), data
    else:
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
def old_verify_downloads(combined_tarlist_path: str, source: str, subtype: str):
    """
    :param combined_tarlist_path: Target list of several sources (CSP/ATLAS/ZTF)
    :param source: Selected source to check for (CSP/ATLAS)
    :param subtype: Wether to check for 1991bg-like or Normals
    :return: None
    """
    print(f"[+++] Verifying downloaded '{source}' '{subtype}' files via {combined_tarlist_path}...")
    tarlist = np.genfromtxt(combined_tarlist_path, delimiter=',', dtype=str, skip_header=1)
    tb = Table(names=tarlist[0], data=tarlist[1:])
    list_sne = []
    for i, row in enumerate(tb):
        n = row['Name'][3:]
        # Path verification
        if subtype.lower() not in ['91bg', 'norm']:
            print(f"[!!!] '{subtype}' is an invalid 'subtype' selection!")
            return
        path = f"data/{source.upper()}-{subtype.lower()}/{source.upper()}{n}.txt"

        # Check path
        if os.path.exists(path):
            list_sne.append(n)
            # print(f"[{len(list_sne)}] {n} exists in {source} download path!")
    if len(list_sne) > 0:
        print(f"      {len(list_sne)} '{source}' files identified out of {len(tb)} '{subtype}' files in target list. "
              f"({len(list_sne)} / {len(tb)})")
        return True
    else:
        print(f"[!!!] Missing data for '{source}' '{subtype}'!")
        return False
def verify_downloads(combined_tarlist_path: str, source: str, subtype: str):
    print(f"[+++] Verifying downloaded '{source}' '{subtype}' files via {combined_tarlist_path}...")

    tarlist = np.genfromtxt(combined_tarlist_path, delimiter=',', dtype=str, skip_header=1)
    tb_tarlist = Table(names=tarlist[0], data=tarlist[1:])
    tb_source = tb_tarlist[tb_tarlist['Data Source(s)'] == source]

    count = 0
    missing = []
    for row in tb_source:
        name = row['Name'][3:]

        if source == 'ATLAS+ZTF':
            atlas_file = f"data/ATLAS-{subtype}/ATLAS{name}.txt"
            ztf_file = f"data/ZTF-{subtype}/ZTF{name}.txt"
            if os.path.exists(atlas_file) and os.path.exists(ztf_file):
                count += 1
            elif os.path.exists(atlas_file):
                missing.append(f"ZTF{name}")
            else:
                missing.append(f"ATLAS{name}")
        else:
            n_file = f"data/{source}-{subtype}/{source}{name}.txt"
            if os.path.exists(n_file):
                count += 1
            else:
                missing.append(name)

    print(f"      ({count}/{len(tb_source)}) {count} files found out of {len(tb_source)} expected from target file.")
    if len(missing) > 0:
        print(f"      ({len(missing)}/{len(tb_source)}) {len(missing)} files missing: {missing}")
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