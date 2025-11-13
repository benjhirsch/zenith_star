import sys
import os
import json

addr_cache_path = os.path.abspath('./cache/addr_cache.json')
datetime_cache_path = os.path.abspath('./cache/datetime_cache.json')

def addr_parse(addr, cache):
    try:
        from geopy.geocoders import Nominatim
    except ImportError:
        print('Cannot parse address without geopy. Run "pip install geopy" from the command line.')
        sys.exit()

    addr_cache = {}
    if os.path.exists(addr_cache_path):
        with open(addr_cache_path) as f:
            addr_cache = json.load(f)
        
    if addr in addr_cache:
        return addr_cache[addr]
    else:
        print('\nCalculating location...')
        geolocator = Nominatim(user_agent='zenith_star')
        obs_loc = geolocator.geocode(addr)
        try:
            obs_lat = obs_loc.latitude
            obs_lon = obs_loc.longitude
        except:
            print('Invalid or unrecognized address.')
            sys.exit()

        if cache:
            addr_cache[addr] = (obs_lat, obs_lon)
            os.makedirs(os.path.dirname(addr_cache_path), exist_ok=True)
            with open(os.path.abspath(addr_cache_path), 'w') as f:
                q = json.dump(addr_cache, f, indent=4)

        return obs_lat, obs_lon

def dt_parse(datetime_arg, cache):
    try:
        from dateutil.parser import parse
    except ImportError:
        print('Cannot parse non-standard date and time without dateutil. Run "pip install python-dateutil" from the command line.')
        sys.exit()

    datetime_cache = {}
    if os.path.exists(datetime_cache_path):
        with open(datetime_cache_path) as f:
            datetime_cache = json.load(f)
        
    if datetime_arg in datetime_cache:
        from datetime import datetime
        return datetime.fromisoformat(datetime_cache[datetime_arg])
    else:
        try:
            local_datetime = parse(datetime_arg)
        except:
            print('Invalid or unrecognized date.')
            sys.exit()

        if cache:
            datetime_cache[datetime_arg] = local_datetime.isoformat()
            os.makedirs(os.path.dirname(datetime_cache_path), exist_ok=True)
            with open(os.path.abspath(datetime_cache_path), 'w') as f:
                q = json.dump(datetime_cache, f, indent=4)

        return local_datetime
