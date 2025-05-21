#!/usr/bin/env python
# Getting redshifts from TNS for SNe with a given RA/dec

from collections import OrderedDict
import requests
import json
import time
from old_file_system.scripts import general as gen

APIkey = gen.get_APIkeys()

tns_bot_id = APIkey['tns_bot_id']
tns_bot_name = APIkey['tns_bot_name']
tns_bot_api_key = APIkey['tns_bot_api_key']

def build_tns_header(tns_bot_id, tns_bot_name):
    """
    Builds the TNS header dictionary.

    Args:
        tns_bot_id (int): Transient name server bot id.
        tns_bot_name (str): Transient name server bot name.
    Returns:
        (dict): Transient name server header dictionary.
    """
    tns_marker = (
        f'tns_marker{{"tns_id": "{int(tns_bot_id)}",'
        f'"type": "bot", "name": "{tns_bot_name}"}}'
    )
    return {"User-Agent": tns_marker}

def query_tns(data, headers, search_url):
    """
    Query the TNS server
    """

    response = requests.post(search_url, headers=headers, data=data)
    response = json.loads(response.text)

    response_message = response.get("id_message")
    response_id_code = response.get("id_code")

    response_status_good = response_id_code == 200

    data = response['data'] if response_status_good else []
    response_reset_time = 10

    response_return = {
        "response_message": response_message,
        "response_id_code": response_id_code,
        "data": data,
        "response_reset_time": response_reset_time,
    }
    return response_return

def rate_limit_query_tns(data, headers, search_url):
    """
    Query TNS but wait if we have reached too many api requests.
    """
    response = query_tns(data, headers, search_url)
    too_many_requests = response["response_id_code"] == 429
    while too_many_requests:
        time_util_rest = response["response_reset_time"]
        time.sleep(time_util_rest + 1)
        response = query_tns(data, headers, search_url)
        too_many_requests = response["response_id_code"] == 429
    return response["data"]

def build_tns_url(tns_api_url, mode=None):
    """
    Builds the url to the tns api service

    Args:
        tns_api_url (str): URL of the Transient name server API.
        mode (str): Which endpoint to access the API. Options are search and get
    Returns:
        (str) Full transient name server api url
    """
    if mode == "search":
        url_end_point = "/search"
    elif mode == "get":
        url_end_point = "/object"
    else:
        raise ValueError("Mode invalid, provide a valid mode (search or get)")
    return tns_api_url + url_end_point

def build_tns_query_data(tns_bot_api_key, data_obj):
    """
    Builds tns search data dictionary.

    Args:
        tns_bot_api_key (str): Transient name server bot api key.
        data_obj (list): List of data representing the tns query.
    Returns:
        (dict): Transient name server query data.

    """
    data_obj = OrderedDict(data_obj)
    return {"api_key": tns_bot_api_key, "data": json.dumps(data_obj)}

def build_tns_search_query_data(tns_bot_api_key, ra, dec, radius=2):
    """
    Build the the search query to find tns transients with a public timestamp
    after time after.

    Args:
        tns_bot_api_key (str): Transient name server bot api key.
        time_after (datetime.datetime): Time to search the transient name
            server for new transients.
    Returns:
        (dict): Transient name server query data.

    """
    search_obj = [
        ("ra", ra),
        ("dec", dec),
        ("radius", str(radius)),
        ("units", "arcsec"),
        ("objname", ""),
        ("objname_exact_match", 0),
        ("internal_name", ""),
        ("internal_name_exact_match ", 0),
        ("objid", ""),
        ("public_timestamp", ""),
    ]
    return build_tns_query_data(tns_bot_api_key, search_obj)

def build_tns_search_query_data_objname(tns_bot_api_key, objname, radius=2):
    """
    Build the the search query to find tns transients with a public timestamp
    after time after.

    Args:
        tns_bot_api_key (str): Transient name server bot api key.
        time_after (datetime.datetime): Time to search the transient name
            server for new transients.
    Returns:
        (dict): Transient name server query data.

    """
    search_obj = [
        ("radius", str(radius)),
        ("units", "arcsec"),
        ("objname", objname),
        ("objname_exact_match", 0),
        ("internal_name", ""),
        ("internal_name_exact_match ", 0),
        ("objid", ""),
        ("public_timestamp", ""),
    ]
    return build_tns_query_data(tns_bot_api_key, search_obj)

def build_tns_get_query_data(tns_bot_api_key, transient):
    """
    Build the the get query data for a TNS transient.

    Args:
        tns_bot_api_key (str): Transient name server bot api key.
        transient (dict): Transient name server transient information.
    Returns:
        (dict): Transient name server query data.

    """
    get_obj = [
        ("objname", transient["objname"]),
        ("objid", transient["objid"]),
        ("photometry", "0"),
        ("spectra", "0"),
    ]
    return build_tns_query_data(tns_bot_api_key, get_obj)

def get_transient_redshift(ra,dec):

    headers = build_tns_header(tns_bot_id, tns_bot_name)
    tns_api_url = f"https://www.wis-tns.org/api/get"

    # get the API URLs
    search_tns_url = build_tns_url(tns_api_url, mode="search")
    get_tns_url = build_tns_url(tns_api_url, mode="get")
    
    search_data = build_tns_search_query_data(tns_bot_api_key, ra, dec)
    transients = rate_limit_query_tns(search_data, headers, search_tns_url)

    get_data = build_tns_get_query_data(tns_bot_api_key, transients[0])
    transient_detail = rate_limit_query_tns(get_data, headers, get_tns_url)

    return transient_detail['redshift']
