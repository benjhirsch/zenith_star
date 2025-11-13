import os
import json
import argparse
import sys
from datetime import datetime
import warnings

import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle
from astropy.time import Time
from astropy.utils.exceptions import AstropyDeprecationWarning
from astroquery.vizier import Vizier

from zsparse import zsparse
from zsconvert import zsconvert

#todo: fixing delta_ra for large search areas near zenith at 0 degrees ra?

def to_angle(theta, theta_format):
    #handles hms/dms strings
    return Angle(theta if isinstance(theta, float) else tuple(float(a) for a in theta.split()), unit=theta_format)

def star_id(obj_id, star):
    #format star id string
    star_id_str = obj_id['format']
    for var in obj_id['columns']:
        star_id_str = star_id_str.replace('{%s}' % var, str(star[var]))
    return star_id_str

def make_float(x):
    try:
        return float(x)
    except:
        return x

warnings.simplefilter('ignore', AstropyDeprecationWarning)
warnings.simplefilter('error', UserWarning)

#load catalog list
with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'zs_catalog_list.json')) as f:
    catalog_list = json.load(f)

vizier_cache_path = os.path.abspath('./cache/vizier_cache.json')

parser = argparse.ArgumentParser(description='Find the zenith star at a given time and location')
parser.add_argument('--datetime', help='Local date and time. Default is now.')
addr_or_latlong = parser.add_mutually_exclusive_group(required=True)
addr_or_latlong.add_argument('--address', help='Local address, e.g. "1600 Pennsylvania Ave NW, Washington, DC 20500"')
addr_or_latlong.add_argument('--latlong', nargs=2, type=float, help='Latitude and longitude, e.g. 38.9 -77.0')
addr_or_latlong.add_argument('--catalog-desc', choices = list(catalog_list.keys()).append('all'), help='Get catalog details.')
addr_or_latlong.add_argument('--clear-cache', action='store_true', help='Clear cached queries.')
parser.add_argument('--limiting-mag', default=6, type=float, help='Limiting magnitude, i.e. minimum brightness for stars. Default is 6, which is the approximate naked-eye limit. Higher values correspond to dimmer stars.')
rad_or_box = parser.add_mutually_exclusive_group()
rad_or_box.add_argument('--search-radius', default=2, type=float, help='Radius in degrees of search query around zenith. Default is 2 degrees.')
rad_or_box.add_argument('--search-box', nargs=2, type=float, help='Height and width in degrees of search box around zenith, e.g. 0.01 0.5')
parser.add_argument('--catalog', choices = catalog_list.keys(), help='Star catalog to query. Default is determined dynamically based on limiting magnitude.')
parser.add_argument('--brightest', action='store_true', help='Return brightest star in search area rather than star closest to zenith.')
parser.add_argument('--display', nargs='?', const=10, type=int, help='Display a visualization of stars near the zenith. Defaults to the 10 closest/brightest.')
parser.add_argument('--disable-cache', action='store_true', help='Include to prevent caching of address, date, and catalog queries. Default stores these for later use.')

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
        print('Limiting magnitude: %s' % catalog_list[catalog]['magnitude']['limit'])
    
    sys.exit()
elif args.address:
    #get latitude and longitude of street address
    obs_lat, obs_lon = zsparse.addr_parse(args.address, not args.disable_cache)
elif args.latlong:
    obs_lat, obs_lon = args.latlong
elif args.clear_cache:
    #clear address, datetime, and vizier caches
    if os.path.exists(zsparse.addr_cache_path):
        os.remove(zsparse.addr_cache_path)
        print('Address cache cleared.')
    if os.path.exists(zsparse.datetime_cache_path):
        os.remove(zsparse.datetime_cache_path)
        print('Datetime cache cleared.')
    if os.path.exists(vizier_cache_path):
        os.remove(vizier_cache_path)
        print('Vizier cache cleared.')
    sys.exit()

#turn lat and long into coordinate object
try:
    obs_earth = EarthLocation(lat=obs_lat*u.deg, lon=obs_lon*u.deg)
except:
    print('Invalid latitude/longitude.')
    sys.exit()

print('Location: Latitude: %.2f\u00b0 %s, Longitude: %.2f\u00b0 %s' % (obs_lat, 'S' if obs_lat < 0 else 'N', abs(obs_lon), 'W' if obs_lon < 0 else 'E'))

#get timezone based on latitude and longitude
lon_time_offset = obs_lon / 15 * u.hour
print('Longitude-based time offset: %s' % lon_time_offset)

#get local datetime
if args.datetime:
    try:
        #if in a nice format, just use datetime module
        local_datetime = datetime.fromisoformat(args.datetime)
    except:
        #otherwise, parse with dateutils
        local_datetime = zsparse.dt_parse(args.datetime, not args.disable_cache)
else:
    #or just use the current time if none given
    local_datetime = datetime.now()

#convert local datetime to utc time
obs_time_utc = Time(local_datetime, scale='utc') - lon_time_offset
print('UTC Time: %s' % obs_time_utc)

#get ra/dec coordinates of zenith at observer location and time
altaz_frame = AltAz(obstime=obs_time_utc, location=obs_earth)
zenith_altaz = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=altaz_frame)
zenith_radec = zenith_altaz.transform_to('icrs')

#determine star catalog to query
if args.catalog:
    catalog = args.catalog
else:
    if args.limiting_mag <= max(c['magnitude']['limit'] for c in catalog_list.values()):
        #if the user-specified limiting magnitude is brighter than the maximum magnitude across all catalogs, return the catalog with the lowest magnitude limit above the limiting mag, with secondary priority for catalogs with higher object count.
        catalog = min((c for c in catalog_list if catalog_list[c]['magnitude']['limit'] >= args.limiting_mag), key=lambda x: (catalog_list[x]['magnitude']['limit'], -catalog_list[x]['obj_count']))
    else:
        #otherwise just return the largest catalog
        catalog = max((c for c in catalog_list), key=lambda x: catalog_list[x]['obj_count'])

if args.search_box:
    search_height, search_width = args.search_box
    search_kwargs = {'height': search_height*u.deg, 'width': search_width*u.deg}
else:
    search_radius = args.search_radius
    search_kwargs = {'radius': search_radius*u.deg}

catalog_specs = catalog_list[catalog]
catalog_mag = catalog_specs['magnitude']
mag_cols = catalog_mag['columns']
ra_col = catalog_specs['RA']['column']
ra_format = catalog_specs['RA']['format']
dec_col = catalog_specs['Dec']['column']
dec_format = catalog_specs['Dec']['format']

#create unique signature for this query to use for caching
vizier_sig = 'cat:%s,mag:%s,box:%s,rad:%s,ra:%s,dec:%s' % (catalog, args.limiting_mag, args.search_box if args.search_box else 'n/a', args.search_radius if not args.search_box else 'n/a', zenith_radec.ra.deg, zenith_radec.dec.deg)

#check cache for previous query results
vizier_cache = {}
if os.path.exists(vizier_cache_path):
    with open(vizier_cache_path) as f:
        vizier_cache = json.load(f)

if vizier_sig in vizier_cache:
    #load cached star list
    star_list = vizier_cache[vizier_sig]
    for star in star_list:
        for col in mag_cols:
            star[col] = make_float(star[col])
        
        if ra_format == 'deg':
            star[ra_col] = make_float(star[ra_col])

        if dec_format == 'deg':
            star[dec_col] = make_float(star[dec_col])
else:
    #create Vizier catalog object with selected catalog's standard columns plus magnitude columns.
    vizier = Vizier(columns=['*', *mag_cols], catalog=catalog)
    vizier.ROW_LIMIT = -1

    #querying catalog for stars within the specified search area around the zenith
    print('\nQuerying VizieR Catalogue Service...')
    catalog_query = vizier.query_region(zenith_radec, **search_kwargs, column_filters={mag_cols[0]: '<%s' % args.limiting_mag})

    try:
        star_list = catalog_query[0]
    except IndexError:
        print('No stars less than magnitude (%s) %s within zenith search area found in %s. Raise limiting magnitude, expand search area, or select a larger catalog.' % (mag_cols[0], args.limiting_mag, catalog_specs['catalog_name']))
        sys.exit()

    if not args.disable_cache:
        #cache star list
        star_cache = []
        for star in star_list:
            star_dic = {}
            for col in catalog_specs['obj_id']['columns']:
                star_dic[col] = str(star[col])
            star_dic[ra_col] = str(star[ra_col])
            star_dic[dec_col] = str(star[dec_col])
            for col in mag_cols:
                star_dic[col] = str(star[col])

            star_cache.append(star_dic)
            
        vizier_cache[vizier_sig] = star_cache
        os.makedirs(os.path.dirname(vizier_cache_path), exist_ok=True)
        with open(os.path.abspath(vizier_cache_path), 'w') as f:
            q = json.dump(vizier_cache, f, indent=4)

#remove stars with calculated visual mag greater than limiting mag
if 'visual_mag_conversion' in catalog_mag:
    under_star_list = []
    for star in star_list:
        vmag = zsconvert.to_vmag(star, mag_cols, catalog_mag['visual_mag_conversion'])
        try:
            if vmag <= args.limiting_mag:
                under_star_list.append(star)
        except:
            pass
    
    star_list = under_star_list

print('%s stars less than magnitude %s within zenith search area found in %s.' % (str(len(star_list)), args.limiting_mag, catalog_specs['catalog_name']))

if args.brightest:
    #get list of star magnitudes
    min_list = [(star[mag_cols[0]], n) for n, star in enumerate(star_list)]
    zenith_text = 'Brightest star near zenith identified.'
else:
    #get coordinates of candidate stars
    ra_list = [to_angle(star[ra_col], ra_format) for star in star_list]
    dec_list = [to_angle(star[dec_col], dec_format) for star in star_list]
    star_coords = SkyCoord(ra=ra_list*u.deg, dec=dec_list*u.deg, frame='icrs')

    #get angle separation between zenith and candidate stars
    min_list = [(zenith_radec.separation(star), n) for n, star in enumerate(star_coords)]
    zenith_text = 'Star closest to zenith identified.'

#find star with minimized parameter (either magnitude or separation), i.e. zenith star
min_idx = min(min_list, key=lambda x: x[0])[1]
zenith_star = star_list[min_idx]
zenith_star_ra = to_angle(zenith_star[ra_col], ra_format) if args.brightest else ra_list[min_idx]
zenith_star_dec = to_angle(zenith_star[dec_col], dec_format) if args.brightest else dec_list[min_idx]
zenith_star_coords = SkyCoord(ra=zenith_star_ra, dec=zenith_star_dec, unit='deg', frame='icrs') if args.brightest else star_coords[min_idx]
zenith_star_altaz = zenith_star_coords.transform_to(altaz_frame)

zenith_star_id = star_id(catalog_specs['obj_id'], zenith_star)

print('\n%s' % zenith_text)
print('\nCatalog ID: %s' % zenith_star_id)
print('Altitude: %.2f\u00b0' % Angle(zenith_star_altaz.alt).deg)
print('Azimuth: %.2f\u00b0' % Angle(zenith_star_altaz.az).deg)
print('RA: %.2f\u00b0' % zenith_star_ra.deg)
print('Dec: %.2f\u00b0' % zenith_star_dec.deg)
print('%s magnitude: %s' % (mag_cols[0], zenith_star[mag_cols[0]]))

if args.display:
    from zsplot import zsplot

    sorted_min_list = sorted(min_list, key=lambda x:x[0])
    display_stars_list = []

    #display lesser of number of found stars or user-given number
    for stat, n in sorted_min_list[:min([len(min_list), args.display])]:
        ra = to_angle(star_list[n][ra_col], ra_format) if args.brightest else ra_list[n]
        dec = to_angle(star_list[n][dec_col], dec_format) if args.brightest else dec_list[n]
        
        coords = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs') if args.brightest else star_coords[n]
        altaz = coords.transform_to(altaz_frame)
        
        mag = stat if args.brightest else star_list[n][mag_cols[0]]

        display_star_id = star_id(catalog_specs['obj_id'], star_list[n])
        
        delta_ra = zsplot.delta_ra(ra, zenith_radec.ra, zenith_radec.dec)
        delta_dec = zsplot.delta_dec(dec, zenith_radec.dec)

        #convert B-V color index or temperature to rgb for display purposes
        if 'temp' in catalog_mag:
            temp = star_list[n][catalog_mag['temp']]
        else:
            if 'b-v' in catalog_mag:
                bv = star_list[n][catalog_mag['b-v']]
            elif 'b-v_conversion' in catalog_mag:
                bv = zsconvert.to_bv(star_list[n], mag_cols, catalog_mag['b-v_conversion'])
            else:
                bv = None

            temp = zsplot.bv_to_temp(bv)

        rgb = zsplot.temp_to_rgb(temp)

        display_stars_list.append({'RA': ra.deg, 'Dec': dec.deg, 'Alt': altaz.alt.deg, 'Az': altaz.az.deg, 'Mag': mag, 'ID': display_star_id, 'delta_RA': -delta_ra.deg, 'delta_Dec': delta_dec.deg, 'rgb': '#%02x%02x%02x' % tuple(rgb), 'Temp': zsplot.round_to(temp, 2)})

    if args.search_box:
        search_x = search_width / 2
        search_y = search_height / 2
    else:
        search_x = search_radius
        search_y = search_radius

    fig = zsplot.zenith_field(display_stars_list, search_x, search_y, zenith_radec.ra.deg, zenith_radec.dec.deg, args.brightest)
    
    fig.show()
