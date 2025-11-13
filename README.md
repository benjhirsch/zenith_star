# Zenith Star
Python tool for finding the zenith star at a given time and location.

Zenith Star is a command line tool that takes as input a time and location, then queries star catalogs to find the star closest to the calculated zenith. It is a Python script run by the Python interpreter and requires Python (3.x) to use. Make sure to place the tool somewhere in your Python PATH.

The basic command for running the tool is:

python zenith_star.py --datetime "[date in any common format, but preferably ISO (YYYY-MM-DD HH:MM:SS)]" --latlong [two numbers, separated by a space, e.g. 38.9 -77.0]

With this command, Zenith Star will calculate the celestial coordinates of the point directly overhead, query the [Bright Star Catalogue](https://en.wikipedia.org/wiki/Bright_Star_Catalogue) for all stars brighter than 6th magnitude (the standard limit of naked-eye vision) within 2° of that point, and return the one closest to it.

## Optional Parameters
--address: Instead of --latlong, provide an address in any common format. Address will be parsed to return latitude and longitude.  
--limiting-mag: Limiting magnitude, i.e. the minimum brightness. Stars above this limit (higher=dimmer) will not be returned. The default is 6, the approximate naked-eye limit.  
--search-radius: Radius in degrees of the search area around the zenith. Default is 2 degrees.  
--search-box: Alternatively, define a search box (two numbers, height and width) in degrees centered on the zenith.  
--brightest: Include to return the brightest star in the search area rather than the one closest to it.  
--catalog: Star catalog to query. The default is determined dynamically based on the limiting magnitude. A list of catalogs can be found in [zs_catalog_list.json](zs_catalog_list.json), which can be expanded.  
--catalog-desc: This option will get the details of one or all the available star catalogs, print them, and exit.  
--display: Include to display a visualization of stars found in the search area, including size and color. Default is 10 if no number given.  
--disable-cache: Include to disable caching of address, datetime, Vizier queries.  
--clear-cache: Include to delete cache files, then exit.

## Usage Tips
* Increasing the limiting magnitude will return many more stars. Be sure to constrain the search area more tightly for higher magnitudes/larger catalogs.
* The coordinates of the zenith are very sensitive to time. 1 minute = 0.25 degrees (about half the width of a full moon) of movement across the sky. For approximate times, you may want to specify a --search-box with a width proportional to the time range.
* Precise geographic coordinates matter less. A 0.25 degree difference on the sky corresponds to ~28 km on the ground.
* When uncertainty is large and many stars are found near the zenith, you may want to return the --brightest star in the search area rather than the one that happens to be closest to the calculated zenith.

## Installation
Clone or otherwise download this repository, then run "pip install -r requirements.txt" from the zenith_star directory. Several modules are optional and not included in requirements.txt. They are:
* geopy: Only necessary to convert the --address parameter to latitude and longitude.
* python-dateutil: Only necessary if you don't use ISO date format ("YYYY-MM-DD HH:MM:SS") with --datetime.
* pandas, plotly, and numpy: Only necessary with the --display option.
* simpleeval: Only necessary for certain catalogs where magnitude values need to be converted from one system to another (those with a "conversion" attribute in [zs_catalog_list.json](zs_catalog_list.json)).
