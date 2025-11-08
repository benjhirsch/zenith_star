# Zenith Star
Python tool for finding the zenith star at a given time and location.

Zenith Star is a command line tool that takes as input a time and location, then queries star catalogs to find the star closest to the calculated zenith. It is a Python script run by the Python interpreter and requires Python (3.x) to use. Make sure to place the tool somewhere in your Python PATH.

The basic command for running the tool is:

python zenith_star.py --datetime "[date in any common format]" --address "[address in common format]"

With this command, Zenith Star will calculate the celestial coordinates of the point directly overhead, query the [Bright Star Catalogue](https://en.wikipedia.org/wiki/Bright_Star_Catalogue) for all stars brighter than 6th magnitude (the standard limit of naked-eye vision) within 2° of that point, and return the one closest to it.

## Optional Parameters
--latlong			Instead of an address, provide exact geographic coordinates, e.g. 38.9 -77.0. The format is two numbers separated by a space.  
--limiting-mag		Limiting magnitude, i.e. the minimum brightness. Stars above this limit (higher=dimmer) will not be returned. The default is 6, the approximate naked-eye limit.  
--search-radius		Radius in degrees of the search area around the zenith. Default is 2 degrees.  
--search-box		Alternatively, define a search box (two numbers, height and width) in degrees centered the zenith.  
--brightest			Include to return the brightest star in the search area rather than the one closest to it.  
--catalog			Star catalog to query. The default is determined dynamically based on the limiting magnitude. A list of catalogs can be found in [za_catalog_list.json](za_catalog_list.json), which can be expanded.  
--catalog-desc		This option will get the details of one or all the available star catalogs, print them, and exit.  
