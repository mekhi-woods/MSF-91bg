import os
import sys
import time as systime
import numpy as np

from astropy.cosmology import FlatLambdaCDM
from astropy.time import Time as astrotime
from astropy.stats import sigma_clipped_stats
from astropy.table import Table


def current_cosmo(H0=70, O_m=0.3):
    return FlatLambdaCDM(H0, O_m)
def get_constants(cosntant_loc='txts/constants.txt'):
    CONSTANTS = {}
    with open(cosntant_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            CONSTANTS.update({line[0]: line[1]})
    return CONSTANTS
def get_APIkeys(apikeys_loc='api_keys.txt'):
    if not os.path.isfile(apikeys_loc):
        raise FileNotFoundError('[!!!] API keys file not found!')
    APIKEY = {}
    with open(apikeys_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            APIKEY.update({line[0]: line[1]})
    return APIKEY
def TNS_details(obj_ra, obj_dec, attempts=10, radius=2, use_keys=True):
    from old_file_system.scripts import tnsAPI
    APIkey = get_APIkeys()
    print('-----------------------------------------------------------------------------------------------------------')
    print('Searching TNS for ['+str(obj_ra)+', '+str(obj_dec)+'] ['+str(radius)+'"]...')
    tns_bot_id, tns_bot_name, tns_bot_api_key = APIkey['tns_bot_id'], APIkey['tns_bot_name'], APIkey['tns_bot_api_key']
    tns_key = dict_handler(choice='unpack', path=get_constants()['tns_key_txt'])
    obj_name, obj_z, obj_discdate = '', 0.00, ''

    # Check key
    print('Checking TNS key...')
    if use_keys and (str(obj_ra) in tns_key):
        obj_name, obj_z, obj_discdate = tns_key[str(obj_ra)]['objname'], tns_key[str(obj_ra)]['z'], tns_key[str(obj_ra)]['discoverydate']
    else:
        # Query TNS key
        print('Not in TNS key, querying TNS...')
        try:
            # Code abridged from David's code
            headers = tnsAPI.build_tns_header(tns_bot_id, tns_bot_name)
            tns_api_url = f"https://www.wis-tns.org/api/get"
            search_tns_url = tnsAPI.build_tns_url(tns_api_url, mode="search")
            get_tns_url = tnsAPI.build_tns_url(tns_api_url, mode="get")
            search_data = tnsAPI.build_tns_search_query_data(tns_bot_api_key, obj_ra, obj_dec, radius=radius)
            transients = tnsAPI.rate_limit_query_tns(search_data, headers, search_tns_url)
            get_data = tnsAPI.build_tns_get_query_data(tns_bot_api_key, transients[0])
            details = tnsAPI.rate_limit_query_tns(get_data, headers, get_tns_url)
            obj_name, obj_z = details['objname'], details['redshift']
            obj_discdate = float(astrotime(details['discoverydate'], format='iso').jd) - 2400000.5 # JD to MJD
            if use_keys:
                with open(get_constants()['tns_key_txt'], 'a') as f:
                    f.write(str(obj_ra) + ', ' + str(obj_dec) + ', ' + obj_name + ', ' + str(obj_z) +
                            ', ' + str(obj_discdate) + '\n') # Save to key
        except Exception as error:
            print('Error:', error)
            if attempts <= 0:
                print('TNS completly timed out.')
                return None, None, None
            print(f'{attempts} attempts remaining')
            systime.sleep(5)
            obj_name, obj_z, obj_discdate = TNS_details(obj_ra, obj_dec, attempts=attempts-1, radius=radius+1, use_keys=use_keys)
    return obj_name, obj_z, obj_discdate
def TNS_get_RA_DEC(objname, attempts=10, radius=2, use_keys=True):
    from old_file_system.scripts import tnsAPI
    APIkey = get_APIkeys()
    print('-----------------------------------------------------------------------------------------------------------')
    print(f'Searching TNS for {objname} [{radius}"]...')
    tns_bot_id, tns_bot_name, tns_bot_api_key = APIkey['tns_bot_id'], APIkey['tns_bot_name'], APIkey['tns_bot_api_key']
    tns_key = dict_handler(choice='unpack', path=get_constants()['tns_key_txt'])
    obj_ra, obj_dec, obj_z, obj_discdate = '', '', 0.00, ''

    known_objname, known_ra = [], []
    for o in tns_key:
        known_objname.append(tns_key[o]['objname'])
        known_ra.append(o)

    # Check key
    print('Checking TNS key...')
    if use_keys and (objname in known_objname):
        obj_ra = known_ra[known_objname.index(objname)]
        obj_dec, obj_z, obj_discdate = tns_key[obj_ra]['dec'], tns_key[obj_ra]['z'], tns_key[obj_ra]['discoverydate']
    else:
        # Query TNS key
        print('Querying TNS...')
        try:
            # Code abridged from David's code
            headers = tnsAPI.build_tns_header(tns_bot_id, tns_bot_name)
            tns_api_url = f"https://www.wis-tns.org/api/get"
            search_tns_url = tnsAPI.build_tns_url(tns_api_url, mode="search")
            get_tns_url = tnsAPI.build_tns_url(tns_api_url, mode="get")
            search_data = tnsAPI.build_tns_search_query_data_objname(tns_bot_api_key, objname, radius=radius)
            transients = tnsAPI.rate_limit_query_tns(search_data, headers, search_tns_url)
            get_data = tnsAPI.build_tns_get_query_data(tns_bot_api_key, transients[0])
            details = tnsAPI.rate_limit_query_tns(get_data, headers, get_tns_url)
            obj_ra, obj_dec, obj_z = details['radeg'], details['decdeg'], details['redshift']
            obj_discdate = float(astrotime(details['discoverydate'], format='iso').jd) - 2400000.5 # JD to MJD
            if use_keys:
                with open(get_constants()['tns_key_txt'], 'a') as f:
                    f.write(str(obj_ra) + ', ' + str(obj_dec) + ', ' + objname + ', ' + str(obj_z) +
                            ', ' + str(obj_discdate) + '\n') # Save to key
        except Exception as error:
            print('Error:', error)
            if attempts <= 0:
                print('TNS completly timed out.')
                return None, None, None, None
            print(f'{attempts} attempts remaining')
            systime.sleep(5)
            obj_name, obj_z, obj_discdate = TNS_get_RA_DEC(objname=objname, attempts=attempts-1, radius=radius+1, use_keys=use_keys)
    return obj_ra, obj_dec, obj_z, obj_discdate
def dict_handler(data_dict={}, choice='', path='default/dict.txt', delimiter=', '):
    if choice == 'unpack':
        print('[+++] Unpacking objects from '+path+'...')
        with open(path, 'r') as f:
            hdr = f.readline()[:-1].split(delimiter)
        data = np.genfromtxt(path, delimiter=delimiter, dtype=str, skip_header=1)
        if len(data) == 0:
            return {}
        temp_objs = {}

        if len(np.shape(data)) == 1:
            iter_period = len(data[0])
        else:
            iter_period = len(data[:, 0])

        for i in range(iter_period):
            obj = data[:, 0][i]
            temp_objs.update({obj: {}})
            for j in range(len(hdr)):
                temp_objs[obj].update({hdr[j]: data[i, j]})
        return temp_objs
    elif choice == 'pack':
        print('[+++] Packing objects into text file...')
        if len(data_dict) == 0:
            print('[!!!] Attempted to pack an empty dictionary!')
            return
        catagories = list(data_dict[list(data_dict.keys())[0]].keys())
        with open(path, 'w') as f:
            f.write('objname')
            for category in catagories:
                f.write(delimiter+category)
            f.write('\n')
            for objname in data_dict:
                f.write(objname)
                for category in catagories:
                    f.write(delimiter+str(data_dict[objname][category]))
                f.write('\n')
        print('[+++] Files packed to', path)
        return
    elif choice == 'arrays':
        print('[+++] Unpacking objects from '+path+' into arraies...')
        with open(path, 'r') as f:
            hdr = f.readline()[:-1].split(delimiter)
        data = np.genfromtxt(path, delimiter=delimiter, dtype=str, skip_header=1)
        if len(data) == 0:
            return {}
        temp_arrs = {}
        for i in range(len(hdr)):
            temp_arrs.update({hdr[i]: data[:, i]})
        return temp_arrs
    else:
        print("[!!!] Invalid packing option! ['arrays'/'pack'/'unpack']")
        return
def default_open(path: str, table_mode: bool = False):
    """
    Open file that has the typical format of this project's 'txt' files.
    :param path: str; location of file
    :param table_mode: bool; whether or not to return the data as an astropy table
    :return: (list, np.array[str]) | astropy.table.Table
    """
    data = np.genfromtxt(path, dtype='str', delimiter=', ')
    hdr, data = data[0, :], data[1:, :]
    for i in range(len(hdr)):
        if hdr[i] in ['objname', 'origin', 'algo']: continue
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
def rm_outliers(mu, z):
    resid = mu - current_cosmo().distmod(z).value
    mn, md, std = sigma_clipped_stats(resid)
    return np.where(abs(resid - mn) < 3 * std)[0]
def get_resid(mu, z):
    return mu.astype(float) - current_cosmo().distmod(z.astype(float)).value
def get_chi2(intercept, x, y, sigma, slope):
    b = intercept[0]  # Extract intercept
    model = slope * x + b
    return np.sum(((y - model) / sigma) ** 2)
def new_get_chi2(m, y, x, b, y_err):
    model_y = (m*x)+b
    return np.sum(((y - model_y) / y_err) ** 2)
def save_figure(save_loc: str):
    if len(save_loc) > 0:
        print(f"Saved figure to...  {save_loc}")
        plt.savefig(save_loc, dpi=300)
    return
def quiet_mode(enabled: bool):
    if enabled: sys.stdout = open(os.devnull, 'w')
    else: sys.stdout = sys.__stdout__
    return
def get_bin_num(data: np.ndarray) -> int:
    # Using Freedman Diaconis Rule: bn = (max-min)/(2*IQR*n^(1/3))
    q1, q3 = np.percentile(data, 25), np.percentile(data, 75)
    iqr = q3 - q1
    return int((np.max(data) - np.min(data)) / (2 * iqr * (len(data)**(-1/3))))

# import numpy as np
# import requests
# from requests.auth import HTTPBasicAuth
# import urllib.request
# import os
#
# API_TOKEN = "7f4e1dee8f52cf0c8cdaf55bf29d04bef4810fb4"
#
# def main():
#
#         #pickle = np.load('tmp.npz', allow_pickle=True)
#         #data = pickle['data']
#
#         data = requests.post(
#             'https://star.pst.qub.ac.uk/sne/atlas4/api/objectlist/',
#             headers={'Authorization': f'Token {API_TOKEN}'},
#             data={'objectlistid':2}
#         ).json()
#
#         np.savez('tmp.npz', data=data)
#
#         count = 0
#         for d in data:
#             if d['observation_status'] is not None and d['observation_status'].startswith('SN Ia'):
#                 print(d['atlas_designation'],d['observation_status'].replace(' ',''),d['ra'],d['dec'])
#                 count += 1
#
#                 ids = d['id']
#                 base_url = 'https://star.pst.qub.ac.uk/sne/atlas4/lightcurveforced/1161048951013729300/'
#                 new_url = base_url.replace('1161048951013729300/',str(ids))
#                 print(new_url)
#
#                 idfile = '/content/ATLAS_ids/' + str(ids)+'.txt'
#                 if os.path.exists(idfile):
#                     continue
#                 urllib.request.urlretrieve(str(new_url), str(idfile))
#                 print(idfile)
#
#             if count > 3000:
#                 break