from geextract.sar.sentinel_one import ts_extract
import unittest
from datetime import datetime

class TestTsExtraction(unittest.TestCase):
    def test_point(self):
        a = ts_extract(
            start = datetime(2019,1,1), 
            end = datetime(2019,6,1),
            lon = -37.7, 
            lat=-4.6)

        self.assertTrue(len(a) == 13)
        self.assertTrue(isinstance(a[0], dict))
        self.assertEqual(set(a[0].keys()),
                         set(['id', 'VV', 'VH']))

    def test_point_radius(self):
        a = ts_extract(
            start = datetime(2019,1,1), 
            end = datetime(2019,6,1),
            lon = -37.7, 
            lat=-4.6,
            radius=500)
        self.assertTrue(len(a) == 13)
        self.assertTrue(isinstance(a[0], dict))
        self.assertEqual(set(a[0].keys()),
                         set(['id', 'VV', 'VH']))

    def test_point_radius_tm(self):
        a = ts_extract(
            start = datetime(2019,1,1), 
            end = datetime(2019,6,1),
            lon = -37.7, 
            lat=-4.6,
            radius=500,
            stats="median")
        self.assertTrue(len(a) == 13)
        self.assertTrue(isinstance(a[0], dict))
        self.assertEqual(set(a[0].keys()),
                         set(['id', 'VV', 'VH']))

    def test_point_radius_descending(self):
        a = ts_extract(
            start = datetime(2019,1,1), 
            end = datetime(2019,6,1),
            lon = -37.7, 
            lat=-4.6,
            radius=500,
            orbit=1,
            stats="mean")
        self.assertTrue(len(a) == 13)
        self.assertTrue(isinstance(a[0], dict))
        self.assertEqual(set(a[0].keys()),
                         set(['id', 'VV', 'VH']))

    def test_exceptions(self):

        # no data for this orbit
        kwargs_1 = {'lon': -37.7,
                  'lat': -4.6,
                  'orbit': 0,
                    'start' : datetime(2019,1,1), 
                    'end' : datetime(2019,6,1),
                    
                    }
        self.assertRaises(ValueError, ts_extract, **kwargs_1)

        # orbit does not exist
        kwargs_1 = {'lon': -3,
                  'lat': 44.7,
                  'orbit': 2,
                  'start': datetime(1999, 1, 1),
                  'end': datetime(2000, 6, 1),
                  'radius': 300}
        self.assertRaises(ValueError, ts_extract, **kwargs_1)

        # mode does not exist
        kwargs_2 = {'lon': -3,
                  'lat': 44.7,
                  'mode': '',
                  'start': datetime(1999, 1, 1),
                  'end': datetime(2000, 6, 1),
                  'radius': 300}
        self.assertRaises(ValueError, ts_extract, **kwargs_2)

        # polar does not exist
        kwargs_3 = {'lon': -3,
                  'lat': 44.7,
                  'polar': 'HH',
                  'start': datetime(1999, 1, 1),
                  'end': datetime(2000, 6, 1),
                  'radius': 300,
                  'stats': 'mode'}
        self.assertRaises(ValueError, ts_extract, **kwargs_3)
