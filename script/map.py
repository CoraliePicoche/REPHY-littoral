import csv
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from math import floor,ceil

a_res='f' #low resolution for now, will switch to fine resolution once we get the right fig
fixed_lat=0.7
fixed_lon=1.4

min_lat_fr=42
min_lon_fr=-5.5
max_lat_fr=51.5
max_lon_fr=8.5

name_site=[]
lon_site=[]
lat_site=[]
f=open('data/local_sites.csv','rU')
reader=csv.reader(f,delimiter=';')
header=reader.next()
for r in reader:
	name_site.append(r[0])
	lon_site.append(r[3])
	lat_site.append(r[4])
f.close()

#BZ
min_lat=47.1
max_lat=min_lat+fixed_lat
min_lon=-3.4
max_lon=min_lon+fixed_lon
m_BZ=Basemap(llcrnrlon=min_lon,llcrnrlat=min_lat,urcrnrlon=max_lon,urcrnrlat=max_lat,resolution=a_res,projection='cyl')

fig=plt.figure(figsize=(10,5))
#ax1=fig.add_subplot(221)
ax1 = plt.subplot2grid((2,2), (0,0), colspan=1,rowspan=1) #Map
plt.annotate("A",xy=(.025,.875), xycoords='axes fraction',weight="bold")
m_BZ.drawcoastlines(color ='k', linewidth =0.2)
m_BZ.fillcontinents(color='lightgrey')
m_BZ.drawparallels(np.arange(floor(min_lat),ceil(max_lat),.5),labels=[1,0,0,0])
m_BZ.drawmeridians(np.arange(floor(min_lon),ceil(max_lon),.5),labels=[0,0,1,0])
for i in range(0,len(name_site)):
        xtmp,ytmp=m_BZ(lon_site[i],lat_site[i])
        m_BZ.plot(xtmp,ytmp,'ko')
        plt.annotate(name_site[i],xy=(xtmp,ytmp), xycoords='data',horizontalalignment='center',verticalalignment='top')

#MO
min_lat=45.7
max_lat=min_lat+fixed_lat
min_lon=-1.9
max_lon=min_lon+fixed_lon
m_MO=Basemap(llcrnrlon=min_lon,llcrnrlat=min_lat,urcrnrlon=max_lon,urcrnrlat=max_lat,resolution=a_res,projection='cyl')

#ax2=fig.add_subplot(222)
ax2 = plt.subplot2grid((2,2), (0,1), colspan=1,rowspan=1) #Map
plt.annotate("B",xy=(.025,.875), xycoords='axes fraction',weight="bold")
m_MO.drawcoastlines(color ='k', linewidth =0.2)
m_MO.fillcontinents(color='lightgrey')
m_MO.drawparallels(np.arange(floor(min_lat),ceil(max_lat),.5),labels=[0,1,0,0])
m_MO.drawmeridians(np.arange(floor(min_lon),ceil(max_lon),.5),labels=[0,0,1,0])
for i in range(0,len(name_site)):
        xtmp,ytmp=m_MO(lon_site[i],lat_site[i])
        m_MO.plot(xtmp,ytmp,'ko')
        plt.annotate(name_site[i],xy=(xtmp,ytmp), xycoords='data',horizontalalignment='center',verticalalignment='top')

#AR
min_lat=44.25
max_lat=min_lat+fixed_lat
min_lon=-1.75
max_lon=min_lon+fixed_lon
m_AR=Basemap(llcrnrlon=min_lon,llcrnrlat=min_lat,urcrnrlon=max_lon,urcrnrlat=max_lat,resolution=a_res,projection='cyl')

#ax3=fig.add_subplot(223)
ax3 = plt.subplot2grid((2,2), (1,0), colspan=1,rowspan=1) #Map
plt.annotate("C",xy=(.025,.875), xycoords='axes fraction',weight="bold")
m_AR.drawcoastlines(color ='k', linewidth =0.2)
m_AR.fillcontinents(color='lightgrey')
m_AR.drawparallels(np.arange(floor(min_lat),ceil(max_lat),.5),labels=[1,0,0,0])
m_AR.drawmeridians(np.arange(floor(min_lon),ceil(max_lon),.5),labels=[0,0,0,1])
for i in range(0,len(name_site)):
        xtmp,ytmp=m_AR(lon_site[i],lat_site[i])
        m_AR.plot(xtmp,ytmp,'ko')
        plt.annotate(name_site[i],xy=(xtmp,ytmp), xycoords='data',horizontalalignment='center',verticalalignment='top')


#SU
min_lat=42.95
min_lon=4.7
max_lat=min_lat+fixed_lat
max_lon=min_lon+fixed_lon

m_SU=Basemap(llcrnrlon=min_lon,llcrnrlat=min_lat,urcrnrlon=max_lon,urcrnrlat=max_lat,resolution=a_res,projection='cyl')

#ax4=fig.add_subplot(2,2,4)
ax4 = plt.subplot2grid((2,2), (1,1), colspan=1,rowspan=1) #Map
plt.annotate("D",xy=(.025,.875), xycoords='axes fraction',weight="bold")
m_SU.drawcoastlines(color ='k', linewidth =0.2)
m_SU.fillcontinents(color='lightgrey')
m_SU.drawparallels(np.arange(floor(min_lat),ceil(max_lat),.5),labels=[0,1,0,0])
m_SU.drawmeridians(np.arange(floor(min_lon),ceil(max_lon),.5),labels=[0,0,0,1])
for i in range(0,len(name_site)):
      	xtmp,ytmp=m_SU(lon_site[i],lat_site[i])
        m_SU.plot(xtmp,ytmp,'ko')
        plt.annotate(name_site[i],xy=(xtmp,ytmp), xycoords='data',horizontalalignment='center',verticalalignment='top')
#m_SU.llcrnrlon = min_lon
#m_SU.urcrnrlon = max_lon
#m_SU.llcrnrlat = min_lat
#m_SU.urcrnrlat = max_lat
m_SU.ax = ax4
#This code is not so nice, this is only a proxy
plt.annotate("10 km",xy=(min_lon+0.05,min_lat+0.05),xytext=(4.873935,42.985),arrowprops=dict(arrowstyle="-"))

plt.subplots_adjust(left=0.05,bottom=0.08,right=0.95,top=0.92,wspace=0, hspace=0.1)
#fig.tight_layout()

axin = inset_axes(m_SU.ax, width="50%", height="50%", loc=1)
        # Global inset map.
inmap = Basemap(projection='cyl', llcrnrlon=min_lon_fr,llcrnrlat=min_lat_fr,urcrnrlon=max_lon_fr,urcrnrlat=max_lat_fr,ax=axin)
inmap.drawcountries(color='white')
inmap.fillcontinents(color='gray')

bx, by = inmap(m_SU.boundarylons, m_SU.boundarylats)
xy = list(zip(bx, by))
mapboundary = Polygon(xy, edgecolor='k', linewidth=1, fill=False)
inmap.ax.add_patch(mapboundary)
plt.annotate("D",xy=(max(bx),max(by)))

bx, by = inmap(m_BZ.boundarylons, m_BZ.boundarylats)
xy = list(zip(bx, by))
mapboundary = Polygon(xy, edgecolor='k', linewidth=1, fill=False)
inmap.ax.add_patch(mapboundary)
plt.annotate("A",xy=(max(bx),max(by)))

bx, by = inmap(m_MO.boundarylons, m_MO.boundarylats)
xy = list(zip(bx, by))
mapboundary = Polygon(xy, edgecolor='k', linewidth=1, fill=False)
inmap.ax.add_patch(mapboundary)
plt.annotate("B",xy=(max(bx),max(by)))

bx, by = inmap(m_AR.boundarylons, m_AR.boundarylats)
xy = list(zip(bx, by))
mapboundary = Polygon(xy, edgecolor='k', linewidth=.75, fill=False)
inmap.ax.add_patch(mapboundary)
plt.annotate("C",xy=(max(bx),max(by)))

for i in range(0,len(name_site)):
        xtmp,ytmp=inmap(lon_site[i],lat_site[i])
        inmap.plot(xtmp,ytmp,'ko',markersize=1.5)

plt.savefig('article/graphe/placing_points.pdf')#,bbox_inches="tight")
