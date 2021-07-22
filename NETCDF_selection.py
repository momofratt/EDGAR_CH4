#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 13:25:51 2021

@author: Cosimo Fratticioli
"""
from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
import numpy as np

###############################################################################
###                         Conversion functions                            ###
###############################################################################
def lat_to_index(lat, lon):
    # returns indexes relative to input latitude and longitude
    # the round command is used to round lat and long to .05, avoiding .00 (emi_ch4 data are centered at .05).
    # it works like this: 0.05 is added to lat and long, then they are rounded to .1, then 0.05 is subtracted again. 
    # the value 0.05 is subtracted once more in order to obtain a .1 rounded number that can be converted to index multiplying by 10 (the grid points spacing is 0.1°)
    lat_index = int( (round((lat + 0.05)*10)*0.1 - 0.1 + 90)*10 )
    lon_index = int( (round((lon + 0.05)*10)*0.1 - 0.1     )*10 )
    return lat_index, lon_index

def index_to_lat(lat_index, lon_index):
    lat = round( (lat_index - 900) * 0.1 + 0.05, 2)
    lon = round( lon_index * 0.1 + 0.05, 2)
    return lat, lon

def from_to(start_lat, start_lon, distance, bearing):
    # gives destination coordinates from start coordinates, bearing (azimuth) and distace
    radius = 6371 
    lat1 = start_lat * (math.pi / 180)
    lon1 = start_lon * (math.pi / 180)
    brng = bearing   * (math.pi / 180)
    lat2 = math.asin(math.sin(lat1) * math.cos(distance / radius) + math.cos(lat1) * math.sin(distance / radius) * math.cos(brng))
    lon2 = lon1 + math.atan2(math.sin(brng) * math.sin(distance / radius) * math.cos(lat1), math.cos(distance / radius) - math.sin(lat1) * math.sin(lat2));
    return lat2*180/math.pi, lon2*180/math.pi
    
###############################################################################
###                         Read file and parameters                        ###
###############################################################################

#################################################
################ USER PARAMETERS ################
file_path = './v6.0_CH4_2018_TOTALS.0.1x0.1/'
file_name = 'v6.0_CH4_2018_TOTALS.0.1x0.1.nc'
WD_1 = 0 # Deg direction of the first line
WD_2 = 75 # Deg direction of the second line
r = 150 # [km] maximum distance from CMN

lat_CMN = 44.19433 # coordinates of the station
lon_CMN = 10.70111
#################################################
#################################################

#################################################
################ FIXED PARAMETERS ###############
########### DO NOT MODIFY THIS SECTION ##########
d_lon = 0.1 # separation of grid points in deg
d_lat = 0.1
data = Dataset(file_path+file_name, 'r')
earth_rad = 6371 # [km]
j_CMN, i_CMN = lat_to_index(lat_CMN, lon_CMN)
lat_CMN, lon_CMN = index_to_lat(j_CMN, i_CMN) # obtain approximated lat_CMN, lon_CMN cendered on grid points 
#################################################
#################################################

###############################################################################
###                         Select emission region                          ###
###############################################################################
# arrays lines. These arrays are used only to plot lines in the final graph
line1_lons = []
line1_lats = []
line2_lons = []
line2_lats = []
line3_lons = []
line3_lats = []

# array for line indexes. they correspond to the lines coordinates in the index reference system. These are used to select the data
line1_ind = []
line2_ind = []
line3_ind = []

# define the extension of the grid that is used to evaluate te data that are included in the region of interest
# evaluate lat and indexes of the four points facing N,E,S,W at distance r from Mt. CMN
lat1,lon1 = from_to(lat_CMN, lon_CMN, r, 0)
lat2,lon2 = from_to(lat_CMN, lon_CMN, r, 90)
lat3,lon3 = from_to(lat_CMN, lon_CMN, r, 180)
lat4,lon4 = from_to(lat_CMN, lon_CMN, r, 270)
j1,i1 = lat_to_index(lat1,lon1)
j2,i2 = lat_to_index(lat2,lon2)
j3,i3 = lat_to_index(lat3,lon3)
j4,i4 = lat_to_index(lat4,lon4)

nbin = max(i2-i4, j1-j3) #maximum extension in number of bins between N-S and E-W directions
if nbin % 2 != 0:
    nbin=nbin+1
    
nbin_2=int(0.5*nbin) # half width of nbin

#################################### line WD_1 #########################################
for j in range(0, int(r/10) + 1 ):
    # define points for the line at angle WD_1. These are sorted from 0 to r
    temp_lat1, temp_lon1 = from_to(lat_CMN, lon_CMN, r*j/(int(r/10)), WD_1)
    temp_j1, temp_i1 = lat_to_index(temp_lat1, temp_lon1)    
    line1_lons.append(temp_lon1)
    line1_lats.append(temp_lat1)
    line1_ind.append((temp_j1-j_CMN+nbin_2, temp_i1-i_CMN+nbin_2))
    
#################################### line WD_2 #########################################
for j in range(0, int(r/10) + 1 ):
    # define points for the line at angle WD_2. These are sorted from r to 0 to obtain continuity in the arc_sector array that is defined later
    temp_lat2, temp_lon2 = from_to(lat_CMN, lon_CMN, r-r*j/(int(r/10)), WD_2)
    temp_j2, temp_i2 = lat_to_index(temp_lat2, temp_lon2)
    line2_lons.append(temp_lon2)
    line2_lats.append(temp_lat2)
    line2_ind.append((temp_j2-j_CMN+nbin_2, temp_i2-i_CMN+nbin_2))
  
################################### arc line (3) #######################################
# define pointd for the arc line
# the if condition is used in order to solve the discontinuity at 360°-0°
if WD_2>WD_1:
    # in this case the arc does not pass across the discontinuity
    for j in range(0, abs(WD_2-WD_1) +1 ):
        temp_lat3, temp_lon3 = from_to(lat_CMN, lon_CMN, r, WD_1+(WD_2-WD_1)*j/(abs(WD_2-WD_1)))
        temp_j3, temp_i3 = lat_to_index(temp_lat3, temp_lon3)
        line3_lons.append(temp_lon3)
        line3_lats.append(temp_lat3)
        line3_ind.append((temp_j3-j_CMN+nbin_2, temp_i3-i_CMN+nbin_2))

else:
    # in this case the arc passes across the discontinuity and another algorithm is needed
    for j in range(0, abs(WD_2 + 360 - WD_1) +1 ):
        alpha = WD_1 + j
        if alpha > 360:
            alpha = alpha-360
        temp_lat3, temp_lon3 = from_to(lat_CMN, lon_CMN, r, alpha)
        temp_j3, temp_i3 = lat_to_index(temp_lat3, temp_lon3)
        line3_lons.append(temp_lon3)
        line3_lats.append(temp_lat3)
        line3_ind.append((temp_j3-j_CMN+nbin_2, temp_i3-i_CMN+nbin_2))
    
####################### find points in the region of interest ###########################

x, y = np.meshgrid(np.arange(nbin), np.arange(nbin)) # make a canvas with coordinates
x, y = x.flatten(), y.flatten()
points = np.vstack((x,y)).T 
arc_sector = line1_ind + line3_ind + line2_ind #defin arc_sector using the arrays of the three lines in the index ref. system 

p = Path(arc_sector) # make a polygon from arc_sector array
grid = p.contains_points(points) # make a grid
mask = grid.reshape(nbin,nbin) # mask with points inside a polygon

###############################################################################
###                      Write output on files                              ###
###############################################################################

# define arrays to store selected data with respective latitudes and longitudes
sel_data=np.zeros((nbin,nbin))
used_lat=np.zeros((nbin,nbin))
used_lon=np.zeros((nbin,nbin))

# output file with three columns (lat, lon, emi_ch4)
outfile = open("selected_emissions_R"+ str(r) + "_WD" + str(WD_1) + "-" + str(WD_2) + ".txt", 'w')
outfile.write("Lat Lon emi_CH4[kg/m2/s]\n")

for j in range(0, nbin):
    for i in range(0,nbin):
        if mask[i,j]:
            sel_data[j][i] = (( data['emi_ch4'][j_CMN-nbin_2 + j ][ i_CMN-nbin_2 + i] ))
            used_lat[j][i], used_lon[j][i] = index_to_lat(j_CMN-nbin_2 +j , i_CMN-nbin_2 +i) 
            outfile.write(str(used_lat[j][i]) + " " + str(used_lon[j][i]) + " " + str(sel_data[j][i]) + "\n")

outfile.close()

# output file with a lat-long matrix which elements are the emission values of CH4 at the respective coordinates
outfile2D = open("selected_emissions_2D_R" + str(r) + "_WD" + str(WD_1) + "-" + str(WD_2) + ".txt", 'w')

for i in range(0, nbin-1):
    _, temp_long = index_to_lat(j_CMN-nbin_2 , i_CMN-nbin_2+i)
    outfile2D.write(str(temp_long) + " ")
i=nbin-1
_, temp_long = index_to_lat(j_CMN-nbin_2 , i_CMN-nbin_2+i)
outfile2D.write(str(temp_long) + "\n")

for j in range(0, nbin):       # loop over latitudes 
    i=0
    temp_lat, _ = index_to_lat(j_CMN-nbin_2 + j , i_CMN-nbin_2)
    outfile2D.write(str(temp_lat) + " ")
    for i in range(1, nbin-1):        # loop over longitudes
        outfile2D.write(str(sel_data[j][i]) + " ")
    i=nbin-1
    outfile2D.write(str(sel_data[j][i]) + "\n")
      
outfile2D.close()

###############################################################################
###                      Evaluate total emission                            ###
###############################################################################
tot_emi_ch4 = 0
earth_rad2 = pow(earth_rad*1000,2) # square earth radius [m^2]
d_lat_rad = d_lat *math.pi/180
d_lon_rad = d_lon *math.pi/180
pi = math.pi

for j in range(0, nbin):    # loop over latitudes to evaluate total emission
    temp_lat, _ = index_to_lat(j_CMN-nbin_2 + j , _ ) # get latitude of jth iteration
    bin_surf = earth_rad2 * d_lat_rad * math.sin( 0.5*pi - abs(temp_lat) *pi/180) * d_lon_rad # [m^2] evaluate bin surface = r*Dtheta * r*sin(theta)*Dphi
    for i in range(1, nbin-1):      # loop over longitudes
        tot_emi_ch4 = tot_emi_ch4 + sel_data[j][i]*bin_surf
        
outfile_emi = open('total_emission_ch4_R' + str(r) + '_WD' + str(WD_1) + '-' + str(WD_2) + '.txt', 'w')
outfile_emi.write("total CH4 emission = " + str(tot_emi_ch4) + " [kg/s]")
outfile_emi.close()

###############################################################################
###                               PLOT                                      ###
###############################################################################

# Get some parameters for the Stereographic Projection
lon_0 = lon_CMN
lat_0 = lat_CMN

m = Basemap(width= 4*r*1000,height=4*r*1000,
            resolution='l',projection='stere',\
            lat_ts=40,lat_0=lat_0,lon_0=lon_0)

dlat = 0.5 + abs((max(lat1,lat3)-lat_CMN))
dlon =   1 + abs((max(lon2,lon4)-lon_CMN))

latinf, loninf = lat_to_index(lat_CMN-dlat, lon_CMN-dlon)
latsup, lonsup = lat_to_index(lat_CMN+dlat, lon_CMN+dlon)

### Plot selected region ###
lats = data['lat'][latinf:latsup]
lons = data['lon'][loninf:lonsup]
emis = data['emi_ch4'][latinf:latsup, loninf:lonsup] 
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)
cs = m.pcolor(xi,yi,np.squeeze(emis), cmap = 'viridis', norm=colors.LogNorm())

# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawrivers()  

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label('[kg m$^-2$ s$^-1$]')

# Add Title
plt.title('CH$_4$ emission [kg m$^{-2}$ s$^{-1}$]')

# Add lines
x_line1, y_line1 = m(line1_lons, line1_lats)
x_line2, y_line2 = m(line2_lons, line2_lats)
x_line3, y_line3 = m(line3_lons, line3_lats)

m.plot(x_line1, y_line1,marker=None,color='m', lw=.5)
m.plot(x_line2, y_line2,marker=None,color='m', lw=.5)
m.plot(x_line3, y_line3,marker=None,color='m', lw=.5)

# Add cities
x_CMN, y_CMN = m(lon_CMN, lat_CMN)
m.plot(x_CMN, y_CMN, marker='^',color='red', label='Mt. Cimone', markersize=2)
x_bolo, y_bolo = m(11.3435, 44.4947)
m.plot(x_bolo, y_bolo, marker='x',color='orange', label='Bologna', markersize=3)

for i in range(0,nbin):
    x_used, y_used = m(used_lon[i], used_lat[i])
    m.plot(x_used, y_used, 'o', color='m', markersize=.5)
    
plt.savefig('CH4_emission_R' + str(r) + '_WD' + str(WD_1) + '-' + str(WD_2) + '.pdf', format='pdf')
plt.show()