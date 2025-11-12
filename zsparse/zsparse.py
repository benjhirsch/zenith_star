import sys

def addr_parse(addr):
    try:
        from geopy.geocoders import Nominatim
    except ImportError:
        print('Cannot parse address without geopy. Run "pip install geopy" from the command line.')
        sys.exit()

    print('\nCalculating location...')
    geolocator = Nominatim(user_agent='zenith_star')
    obs_loc = geolocator.geocode(addr)
    try:
        obs_lat = obs_loc.latitude
        obs_lon = obs_loc.longitude
        return obs_lat, obs_lon
    except:
        print('Invalid or unrecognized address.')
        sys.exit()

def dt_parse(datetime, timezone):
    try:
        from dateutil.parser import parse
    except ImportError:
        print('Cannot parse non-standard date and time without dateutil. Run "pip install python-dateutil" from the command line.')
        sys.exit()
    
    try:
        local_datetime = parse(datetime).replace(tzinfo=timezone)
        return local_datetime
    except:
        print('Invalid or unrecognized date.')
        sys.exit()