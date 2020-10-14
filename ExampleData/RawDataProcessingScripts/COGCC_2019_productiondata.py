"""
COGCC_2019_productiondata.py loads well location and production data that was downloaded from the COGCC.
The data are used to calculate annual production per site, and the resulting array is saved.
The COGCC provides production and location data by well. In order to group wells into sites, any wells within
50 meters of each other are assumed to be part of the same site.
"""

import pandas as pd
import numpy as np
import pysal as ps
import pysal.lib
import geopy.distance
from os.path import dirname, abspath
import pickle
from feast.input_data_classes import ProductionData
import os

rsc_path, _ = os.path.split(dirname(abspath(__file__)))
prod_path = os.path.join(rsc_path, 'RawData', 'COGCC-2019-Production.xlsx')
loc_path = os.path.join(rsc_path, 'RawData', 'COGCC-2019-well-locations.dbf')
file_out = os.path.join(rsc_path, 'DataObjectInstances', 'COGCC_site_prod_2019.p')


# Load location data
db = pysal.lib.io.open(loc_path)
d = {str.upper(col): db.by_col(col) for col in db.header}
wells = pd.DataFrame(d)
db.close()

# load production data
dat = pd.read_excel(prod_path)

# extract relevant data from the production data frame
sqnum = np.array(dat['api_seq_num'], dtype=str)
ccode = np.array(dat['api_county_code'], dtype=str)
sdtrack = np.array(dat['sidetrack_num'], dtype=str)
formcode = np.array(dat['formation_code'], dtype=str)
opnum = np.array(dat['operator_num'], dtype=str)

# Store county and sequence codes to a list
cscode = []
for ind in range(len(sqnum)):
    c, s = ccode[ind], sqnum[ind]
    while len(c) < 3:
        c = '0' + c
    while len(s) < 5:
        s = '0' + s
    cscode.append(c + s)

# extract relevant data from the location data object
cond = wells.FACIL_STAT == 'PR'
llsqnum = np.array(wells.API_SEQ[cond], dtype='str')
llcounty = np.array(wells.API_COUNTY[cond], dtype='str')
lats, longs, APIs = list(wells.LATITUDE[cond]), list(wells.LONGITUDE[cond]), list(wells.API[cond])
sortargs = np.argsort(lats)
lats = [lats[l_index] for l_index in sortargs]
longs = [longs[l_index] for l_index in sortargs]
APIs = [APIs[l_index] for l_index in sortargs]

# Determine which wells in the production data set are not in the location data set
missing_codes = []
for ap in APIs:
    if ap not in cscode:
        missing_codes.append(ap)

cond = np.ones(len(APIs), dtype=bool)
for ind in range(len(APIs)):
    if APIs[ind] in missing_codes:
        cond[ind] = 0
print("The number of excluded indexes is {:0.0f} and should be 144.".
      format(len(cond) - np.sum(cond)))

api = np.array(APIs)[cond]
lat = np.array(lats)[cond]
long = np.array(longs)[cond]
sortargs = np.argsort(lat)
lat = list(lat[sortargs])
long = list(long[sortargs])
api = list(api[sortargs])


# noinspection PyShadowingNames
def distance(places, point):
    """
    calculates the distance between an array of places and a point.
    :input places: an array or list of latitudes and longitudes (deg)
    :input point: the latitude and longitude of a point (deg)
    :return dist: a list of the distances between the elements of "places"
                    and "point" (km)
    """
    dist = np.zeros(len(places))
    for ind in range(len(places)):
        dist[ind] = geopy.distance.distance(places[ind][:2],
                                            point[:2]).km
    return dist


def crude_lat_dist(l1, l2):
    """
    An approximation of North-South distance between latitudes
    Intended to give a fast rejection of points that don't satisfy
    a maximum distance criteria
    :input l1: latitude 1 (degrees)
    :input l2: latitude 2 (degrees)
    :return: the north south distance between the latitudes
    """

    earth_radius = 6373  # km
    return earth_radius * np.pi / 180 * np.abs(l1 - l2)  # km


def binary_find_min(a, ls):
    """
    Performs a binary search to find the minimum difference between "a"
    in the elements in the list "ls"
    :input a: a scalar
    :input ls: a sorted list of scalars
    :output: the index of the value in ls that is closest to a
    """
    top = len(ls) - 1
    bottom = 0
    while top - bottom > 1:
        mid = int((top - bottom) / 2 + bottom)
        if ls[mid] > a:
            top = mid
        elif ls[mid] < a:
            bottom = mid
        else:
            return mid
    if np.abs(ls[top] - a) < np.abs(ls[bottom] - a):
        return top
    else:
        return bottom


def finder(ll, lats, longs, APIs, pad, max_dist=0.05):
    """
    Determine which lats and longs are within the criteria
    distance of ll to be included in the same pad
    The function is designed for use in a recursive program
    :input ll: latitude and longitude of a point to be considered
    :input lats: a sorted array of latitudes
    :input longs: an array of longitudes associated with lats
    :input APIs: an array of API values of the wells
    :input pad: a list of wells grouped into pads.
    :input max_dist: maximum distance to the nearest well to be considered
    a pad (km)
    :return: none
    """
    pad.append(ll)
    if len(lats) < 1:
        return
    minind = binary_find_min(ll[0], lats)
    # preliminary fast culling of distant wells
    if crude_lat_dist(ll[0], lats[minind]) > max_dist * 10:
        return
    else:
        cond = []
        temp = minind
        while temp >= 0 and \
                crude_lat_dist(ll[0], lats[temp]) <= max_dist * 10:
            cond.append(temp)
            temp -= 1
        temp = minind + 1
        while temp < len(lats) and \
                crude_lat_dist(ll[0], lats[temp]) <= max_dist * 10:
            cond.append(temp)
            temp += 1
        cond.sort(reverse=True)
        # final selection of nearby wells to be included in the pad
        la = [lats[l_index] for l_index in cond]
        lo = [longs[l_index] for l_index in cond]
        distances = distance(np.array(list(zip(la, lo))), ll)
        cond2 = np.where(distances < max_dist)[0]
        if len(cond) == 0:
            return
        winners = []
        for c2 in cond2:
            winners.append([lats[cond[c2]], longs[cond[c2]],
                            APIs[cond[c2]]])
            del (lats[cond[c2]])
            del (longs[cond[c2]])
            del (APIs[cond[c2]])
        for w in winners:
            finder(w, lats, longs, APIs, pad)


pad = []
while len(lat) > 0:
    if len(lat) % 1000 == 0:
        print(len(lat))
    ll = [lat[0], long[0], api[0]]
    del (lat[0], long[0], api[0])
    pad.append([])
    finder(ll, lat, long, api, pad[-1])

site_prod = []
cscode = np.array(cscode)
for site in pad:
    lats = []
    longs = []
    site_prod.append([0, 0, 0, 0])
    for well in site:
        cond = np.where(cscode == well[2])[0]
        site_prod[-1][2] += np.sum(dat['gas_prod'][cond])
        site_prod[-1][3] += np.sum(dat['oil_prod'][cond])
        lats.append(well[0])
        longs.append(well[1])
    site_prod[-1][0] = np.mean(lats)
    site_prod[-1][1] = np.mean(longs)
site_prod = np.array(site_prod)

raw_data_file = ['COGCC-2019-Production.xlsx', 'COGCC-2019-well-locations.dbf']
notes = """"
This object stores well site location and production data developed using data from the COGCC.
Well location and production data that were downloaded from the COGCC.

Production and location data by well for 2019 were downloaded from https://cogcc.state.co.us/data2.html#/downloads. 
The downloads were completed on Jan. 13, 2020. In order to group wells into sites, any wells within
50 meters of each other were assumed to be part of the same site.

site_prod.site_prod contains 4 columns: latitude (degrees), longitude (degrees), gas production (mcf/year), oil 
production (bbl/year)
"""
site_prod = ProductionData(notes=notes, raw_file_name=raw_data_file, data_prep_file='COGCC_2019_productiondata.py',
                     site_prod=site_prod)

pickle.dump(site_prod, open(file_out, 'wb'))
