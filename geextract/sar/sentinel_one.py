"""geextract.sar.sentinelone"""

import ee
import sqlite3
import pandas as pd
import re
from datetime import datetime
import warnings
from enum import Enum

from geextract import simplify, dictify

# Silence pandas warning
warnings.simplefilter(action='ignore')

ee.Initialize()

SATELLITES = ["S1A", "S1B"]

# Interferometric Wide swath (IW), Extra-Wide swath (EW), StripMap (SM)
MODES = ["IW", "EW", "SM"]

POLARIZATIONS = [
    "VV",
    "VH"
]

ORBITS = ["ASCENDING", "DESCENDING"]

class Orbit(Enum):
    ASCENDING = 0
    DESCENDING = 1


def get_date(filename, start_only=True):
    """Retrieve date information from typical Sentinel-1 GRD filenames

    Args:
        filename (str): Sentinel-1 GRD file name
        start_only (boolean): If True, return only a datetime object for acquisition start

    Returns:
        datetime.date : The corresponding acquisition starting date of the filename.
        datetime.date (optional) : The corresponding acquisition ending date of the filename


    Examples:
        >>> from geextract.sar import sentinel_one
        >>> sentinel_one.get_date('S1A_IW_GRDH_1SDV_20190105T000740_20190105T000805_025335_02CDC8_E4C2')
    """
    p0 = re.compile(r'''
        (?P<sensor>S1A|S1B)_
        (?P<acquisition_mode>IW|SW|SM)_GRD
        (?P<resolution>M|H)_
        (?P<processing_level>\d{1})SDV_
        (?P<start_date>\d{8})T
        (?P<start_time>\d{6})_
        (?P<end_date>\d{8})T
        (?P<end_time>\d{6})_
        (?P<absolute_orbit_number>\w{6})_
        (?P<mission_data_take_id>\w{6})_
        (?P<product_unique_id>\w{4})
        ''')
    if p0.search(filename):
        m = p0.search(filename)
        start_d = datetime.strptime(f"{m.group('start_date')}T{m.group('start_time')}", '%Y%m%dT%H%M%S').date()
        if not start_only:
            end_d = datetime.strptime(f"{m.group('end_date')}T{m.group('end_time')}", '%Y%m%dT%H%M%S').date()
            return (start_d, end_d)
        else:
            return start_d
    else:
        raise ValueError('Unknown pattern')


def ts_extract(start, mode = "IW", lon = None, lat = None,
               end = datetime.today(), radius = None, feature = None,
               polar = POLARIZATIONS, orbit = None, stats = 'mean'):
    """Perform a spatio temporal query to extract Sentinel 1 GRD data
        from gee

    Args:
        lon (float): Center longitude in decimal degree
        lat (float): Center latitude in decimal degree
        mode (str): Acquisition mode. Needs to be in ``['IW', 'EW', 'SM']``
            Default is ``'IW'``
        start (datetime.datetime): Start date of first acquisition
        end (datetime.datetime): Optional end date of last acquisition; automatically set as today if unset
        radius (float): Optional radius around center point in meters. If unset,
            time-series of a single pixel are queried. Otherwise a reducer is used
            to spatially aggregate the pixels intersecting the circular feature
            built.
        feature (dict): Optional dictionary representation of a polygon feature
            in longlat CRS. If unset, time-series of a single pixel are queried.
            Otherwise a reducer is used to spatially aggregate the pixels intersecting
            the given feature.
        polar (list of str or str): Chosen polarization. Can be multiple polarizations for PolSAR processing. 
            Default is ``['VV', 'VH]``
        orbit (int or None): Chosen orbit type (0: ASCENDING, 1: DESCENDING). If None, then data from both orbit types will be acquired. If int, we filter the ImageCollection by orbit type. 
            Default is ``None``.
        stats (str): Spatial aggregation function to use. Only relevant
            if a radius value is set.

    Returns:
        dict: A dictionary representation of the json data returned by the gee platform.

    Example:
        >>> from geextract.sar import sentinel_one
        >>> from pprint import pprint
        >>> from datetime import datetime

        >>> lon = -89.8107197
        >>> lat = 20.4159611

        >>> out = sentinel_one.ts_extract(lon=lon, lat=lat, mode = 'IW', start=datetime(2018, 1, 1, 0, 0),
        >>>                            radius=500)
        >>> pprint(out)

    """

    ################
    # CHECK INPUTS #
    ################

    # Mode
    if mode not in MODES:
        raise ValueError(f"Unknown mode '{mode}' (Must be one of {MODES})")

    # Polar
    if polar is None:
        polar = POLARIZATIONS
    elif (isinstance(polar, list)):
        for p in polar:
            if p not in POLARIZATIONS:
                raise ValueError(f"Unknown polarization '{p}' (Must be one of {POLARIZATIONS})")
    else:
        # To have cleaner code once building GEE query
        polar = [polar]

    # Orbit
    if orbit is not None:
        if orbit not in [0,1]:
            raise ValueError(f"Unknown orbit type '{orbit}' (Must be one of {ORBITS} where 0 -> ASCENDING, 1 -> DESCENDING)")
        else:
            orbit = Orbit(orbit).name

    ###########################
    # INITIALIZING COLLECTION #
    ###########################

    collection_name = 'COPERNICUS/S1_GRD'
    sentinel_ic = ee.ImageCollection(collection_name) \
            .filterDate(start=start, opt_end=end)
    ########################
    # FILTERING COLLECTION #
    ########################
    
    # Mode
    sentinel_ic = sentinel_ic \
        .filter(ee.Filter.eq('instrumentMode', mode))

    # Orbit
    if orbit != None:
        sentinel_ic = sentinel_ic \
            .filter(ee.Filter.eq('orbitProperties_pass', orbit))
    
    if radius is not None or feature is not None:
        # Define spatial aggregation function
        if stats == 'mean':
            fun = ee.Reducer.mean()
        elif stats == 'median':
            fun = ee.Reducer.median()
        elif stats == 'max':
            fun = ee.Reducer.max()
        elif stats == 'min':
            fun = ee.Reducer.min()
        else:
            raise ValueError('Unknown spatial aggregation function. Must be one of mean, median, max, or min')

        if feature is not None:
            geometry = ee.Geometry.Polygon(feature['geometry']['coordinates'])
        else: # Geometry defined by point and radius
            geometry = ee.Geometry.Point(lon, lat).buffer(radius)
        
        # Define function to map over imageCollection to perform spatial aggregation 
        def _reduce_region(image):
            """Spatial aggregation function for a single image and a polygon feature"""
            stat_dict = image.reduceRegion(fun, geometry, 10)
            # Feature needs to be rebuilt because the backend doesn't accept to map
            # functions that return dictionaries
            return ee.Feature(None, stat_dict)
        try:
            fc = sentinel_ic.filterBounds(geometry).select(*tuple(polar)).map(_reduce_region).getInfo()
        except ee.EEException as e:
            raise ValueError("No data found for these parameters. Please change the orbit type or the coordinates.") from None
        out = simplify(fc)
    else:
        # Extraction with a point, no spatial aggregation, etc
        geometry = ee.Geometry.Point(lon, lat)
        try:
            l = sentinel_ic.filterBounds(geometry).select(*tuple(polar)).getRegion(geometry, 10).getInfo()
        except ee.EEException as e:
            raise ValueError("No data found for these parameters. Please change the orbit type or the coordinates.") from None
        out = dictify(l)
        # pop longitude and latitude keys from dict collection so that band aliases can
        # be replaced by their color names
        [d.pop('longitude', None) for d in out]
        [d.pop('latitude', None) for d in out]
        [d.pop('time', None) for d in out]
    return out
