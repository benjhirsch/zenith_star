import os
import json
import argparse
import sys
from geopy.geocoders import Nominatim
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle
from timezonefinder import TimezoneFinder
from zoneinfo import ZoneInfo
from datetime import datetime, timezone
from dateutil.parser import parse
from astropy.time import Time
from astroquery.vizier import Vizier
import simpleeval
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning

def is_float(x):
    try:
        y = float(x) 
        return True
    except:
        return False
    
def to_angle(theta, theta_format):
    return Angle(theta if isinstance(theta, float) else tuple(float(a) for a in theta.split()), unit=theta_format)

warnings.simplefilter('ignore', AstropyDeprecationWarning)
warnings.simplefilter('error', UserWarning)

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'zs_catalog_list.json')) as f:
    catalog_list = json.load(f)

parser = argparse.ArgumentParser(description='Find the zenith star at a given time and location')
parser.add_argument('--datetime', help='Local date and time. Default is now.')
addr_or_latlong = parser.add_mutually_exclusive_group(required=True)
addr_or_latlong.add_argument('--address', help='Local address, e.g. "1600 Pennsylvania Ave NW, Washington, DC 20500"')
addr_or_latlong.add_argument('--latlong', nargs=2, help='Latitude and longitude, e.g. 38.9 -77.0')
addr_or_latlong.add_argument('--catalog-desc', choices = list(catalog_list.keys()).append('all'), help='Get catalog details.')
parser.add_argument('--limiting-mag', default=6, help='Limiting magnitude, i.e. minimum brightness for stars. Default is 6, which is the approximate naked-eye limit. Higher values correspond to dimmer stars.')
rad_or_box = parser.add_mutually_exclusive_group()
rad_or_box.add_argument('--search-radius', default=2, help='Radius in degrees of search query around zenith. Default is 2 degrees.')
rad_or_box.add_argument('--search-box', nargs=2, help='Height and width in degrees of search box around zenith, e.g. 0.01 0.5')
parser.add_argument('--catalog', choices = catalog_list.keys(), help='Star catalog to query. Default is determined dynamically based on limiting magnitude.')
parser.add_argument('--brightest', action='store_true', help='Return brightest star in search area rather than star closest to zenith.')
args = parser.parse_args()

if args.catalog_desc:
    #if --catalog-desc argument passed, print catalog description and then exit
    if not args.catalog_desc == 'all':
        catalog_list = {args.catalog_desc: catalog_list[args.catalog_desc]}
    
    for catalog in catalog_list:
        print()
        print(catalog)
        print(catalog_list[catalog]['catalog_name'])
        print('Object count: %s' % catalog_list[catalog]['obj_count'])
        print('Limiting magnitude: %s' % catalog_list[catalog]['visual_mag']['limit'])
    
    sys.exit()
elif args.address:
    #get latitude and longitude of street address
    print('\nCalculating location...')
    geolocator = Nominatim(user_agent='zenith_star')
    obs_loc = geolocator.geocode(args.address)
    obs_lat = obs_loc.latitude
    obs_lon = obs_loc.longitude
elif args.latlong:
    obs_lat, obs_lon = map(float, args.latlong)

obs_earth = EarthLocation(lat=obs_lat*u.deg, lon=obs_lon*u.deg)
print('Location: Latitude: %.2f\u00b0 %s, Longitude: %.2f\u00b0 %s' % (obs_lat, 'S' if obs_lat < 0 else 'N', abs(obs_lon), 'W' if obs_lon < 0 else 'E'))

#get timezone based on latitude and longitude
tzf = TimezoneFinder()
tz_str = tzf.timezone_at(lng=obs_lon, lat=obs_lat)
obs_timezone = ZoneInfo(tz_str)
print('Timezone: %s' % tz_str)

#get local datetime
if args.datetime:
    local_datetime = parse(args.datetime).replace(tzinfo=obs_timezone)
else:
    local_datetime = datetime.now(tzinfo=obs_timezone)

#convert local datetime to utc time
obs_time_utc = Time(local_datetime.astimezone(timezone.utc), scale='utc')
print('UTC Time: %s' % obs_time_utc)

#get ra/dec coordinates of zenith at observer location
altaz_frame = AltAz(obstime=obs_time_utc, location=obs_earth)
zenith_altaz = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=altaz_frame)
zenith_radec = zenith_altaz.transform_to('icrs')

#determine star catalog to query
limiting_mag = float(args.limiting_mag)
if args.catalog:
    catalog = args.catalog
else:
    if limiting_mag <= max(c['visual_mag']['limit'] for c in catalog_list.values()):
        #if the user-specified limiting magnitude is brighter than the maximum magnitude across all catalogs, return the catalog with the lowest magnitude limit above the limiting mag, with secondary priority for catalogs with higher object count.
        catalog = min((c for c in catalog_list if catalog_list[c]['visual_mag']['limit'] >= limiting_mag), key=lambda x: (catalog_list[x]['visual_mag']['limit'], -catalog_list[x]['obj_count']))
    else:
        #otherwise just return the largest catalog
        catalog = max((c for c in catalog_list), key=lambda x: catalog_list[x]['obj_count'])

catalog_specs = catalog_list[catalog]

#create Vizier catalog object with selected catalog's standard columns plus magnitude columns.
vizier = Vizier(columns=['*', *catalog_specs['visual_mag']['columns']], catalog=catalog)
vizier.ROW_LIMIT = -1

#querying catalog for stars within the specified search area around the zenith
print('\nQuerying VizieR Catalogue Service...')
if args.search_box:
    search_height, search_width = map(float, args.search_box)
    catalog_query = vizier.query_region(zenith_radec, height=search_height*u.deg, width=search_width*u.deg, column_filters={catalog_specs['visual_mag']['columns'][0]: '<%s' % args.limiting_mag})
else:
    search_radius =  float(args.search_radius)
    catalog_query = vizier.query_region(zenith_radec, radius=search_radius*u.deg, column_filters={catalog_specs['visual_mag']['columns'][0]: '<%s' % args.limiting_mag})

try:
    star_list = catalog_query[0]
except IndexError:
    print('No stars less than magnitude (%s) %s within zenith search area found in %s. Raise limiting magnitude, expand search area, or select a larger catalog.' % (catalog_specs['visual_mag']['columns'][0], args.limiting_mag, catalog_specs['catalog_name']))
    sys.exit()

#remove stars with calculated visual mag greater than limiting mag
if 'conversion' in catalog_specs['visual_mag']:
    se = simpleeval.SimpleEval()
    under_star_list = []
    for star in star_list:
        for var in catalog_specs['visual_mag']['columns']:
            se.names[var] = star[var]

        if all([is_float(star[var]) for var in catalog_specs['visual_mag']['columns']]):
            vmag = se.eval(catalog_specs['visual_mag']['conversion'])
        elif is_float(star[catalog_specs['visual_mag']['columns'][0]]):
            vmag = star[catalog_specs['visual_mag']['columns'][0]]
        
        try:
            if vmag <= limiting_mag:
                under_star_list.append(star)
        except:
            pass
    
    star_list = under_star_list

print('%s stars less than magnitude %s within zenith search area found in %s.' % (str(len(star_list)), args.limiting_mag, catalog_specs['catalog_name']))

if args.brightest:
    #get list of star magnitudes
    min_list = [(star[catalog_specs['visual_mag']['columns'][0]], n) for n, star in enumerate(star_list)]
    zenith_text = 'Brightest star near zenith identified.'
else:
    #get coordinates of candidate stars
    ra_list = [to_angle(star[catalog_specs['RA']['column']], catalog_specs['RA']['format']) for star in star_list]
    dec_list = [to_angle(star[catalog_specs['Dec']['column']], catalog_specs['Dec']['format']) for star in star_list]
    star_coords = SkyCoord(ra=ra_list*u.deg, dec=dec_list*u.deg, frame='icrs')

    #get angle separation between zenith and candidate stars
    min_list = [(zenith_radec.separation(star), n) for n, star in enumerate(star_coords)]
    zenith_text = 'Star closest to zenith identified.'

#find star with minimized parameter (either magnitude or separation), i.e. zenith star
min_idx = min(min_list, key=lambda x: x[0])[1]
zenith_star = star_list[min_idx]
zenith_star_ra = to_angle(zenith_star[catalog_specs['RA']['column']], catalog_specs['RA']['format']) if args.brightest else ra_list[min_idx]
zenith_star_dec = to_angle(zenith_star[catalog_specs['Dec']['column']], catalog_specs['Dec']['format']) if args.brightest else dec_list[min_idx]
zenith_star_coords = SkyCoord(ra=zenith_star_ra, dec=zenith_star_dec, unit='deg', frame='icrs') if args.brightest else star_coords[min_idx]
zenith_star_altaz = zenith_star_coords.transform_to(altaz_frame)

print('\n%s' % zenith_text)
print('\nCatalog ID: %s %s' % (catalog_specs['id'], zenith_star[catalog_specs['id']]))
print('Altitude: %.2f\u00b0' % Angle(zenith_star_altaz.alt).deg)
print('Azimuth: %.2f\u00b0' % Angle(zenith_star_altaz.az).deg)
print('RA: %.2f\u00b0' % zenith_star_ra.deg)
print('Dec: %.2f\u00b0' % zenith_star_dec.deg)
print('%s magnitude: %s' % (catalog_specs['visual_mag']['columns'][0], zenith_star[catalog_specs['visual_mag']['columns'][0]]))