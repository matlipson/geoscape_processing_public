import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt

import cartopy.geodesic as cgeo

import os
oshome=os.getenv('HOME')

res = 'CCI'
domain = 'sydney'
version = 'v1.01'

################################################################################

projpath = f'.'
datapath = f'{projpath}/data'
plotpath = f'{projpath}/figures'

def main():

    ds = get_geoscape_data(domain,res,version)

    plot_fig5(ds,domain,res)
    plot_landcover(ds,domain,res)
    plot_morphology(ds,domain,res)

    return

def plot_fig5(ds,domain,res):

    plt.close('all')
    fig,ax = plt.subplots(ncols=2,nrows=4,figsize=(8,12),sharex=True,sharey=True)

    ds['total_built'].plot(ax=ax[0,0],cmap='Greys', vmin=0,vmax=1,cbar_kwargs={'label':''})
    ds['total_vegetation'].plot(ax=ax[1,0],cmap='Greens',vmin=0,vmax=1,cbar_kwargs={'label':''})
    ds['building_fraction'].plot(ax=ax[2,0],cmap='Reds',vmin=0,vmax=1,cbar_kwargs={'label':''})
    ds['roadpath_fraction'].plot(ax=ax[3,0],cmap='Blues',vmin=0,vmax=1,cbar_kwargs={'label':''})

    ax[0,0].set_title('Total built fraction',fontsize=11)
    ax[1,0].set_title('Total vegetation fraction',fontsize=11)
    ax[2,0].set_title('Building plan fraction',fontsize=11)
    ax[3,0].set_title('Road and paths fraction',fontsize=11)

    cmap = 'plasma'

    ds['building_height'].plot(ax=ax[0,1],cmap=cmap, vmin=0, vmax=25,extend='max',cbar_kwargs={'label':''})
    ds['frontal_density'].plot(ax=ax[1,1],cmap=cmap,vmin=0, vmax=0.5,extend='max',cbar_kwargs={'label':''})
    ds['height_to_width'].plot(ax=ax[2,1],cmap=cmap,vmin=0, vmax=2,extend='max',cbar_kwargs={'label':''})
    ds['roughness_mac'].plot(ax=ax[3,1],cmap=cmap,vmin=0, vmax=1,extend='max',cbar_kwargs={'label':''})

    ax[0,1].set_title('Building height (mean)',fontsize=11)
    ax[1,1].set_title('Frontal area density',fontsize=11)
    ax[2,1].set_title('Canyon aspect ratio',fontsize=11)
    ax[3,1].set_title('Roughness length (Mac1998)',fontsize=11)

    lat = ds.latitude.max().values - 0.09
    lon = ds.longitude.min().values + 0.06

    for a in ax.flatten():
        a.set_xlabel('')
        a.set_ylabel('')
        a = distance_bar(lat,lon,a1=a)
        
    for i in [0,1,2,3]:
        ax[i,0].set_ylabel('latitude')

    for i in [0,1]:
        ax[3,i].set_xlabel('longitude')

    fig.subplots_adjust(wspace=0.07)
    
    # plt.show()
    fig.savefig(f'{plotpath}/{domain}_fig5_{res}_{version}.png',bbox_inches='tight',dpi=300)

    return

def plot_landcover(ds,domain,res):

    # cmap = 'plasma'

    cmap = create_cmap('Greys',1000,start_col='white')

    plt.close('all')
    fig,ax = plt.subplots(ncols=2,nrows=4,figsize=(8,12),sharex=True,sharey=True)

    ds['total_built'].plot(ax=ax[0,0],cmap=cmap, vmin=0,vmax=1,cbar_kwargs={'label':''})
    ds['building_fraction'].plot(ax=ax[1,0],cmap='Reds',vmin=0,vmax=1,cbar_kwargs={'label':''})
    ds['roadpath_fraction'].plot(ax=ax[2,0],cmap='Reds',vmin=0,vmax=1,cbar_kwargs={'label':''})
    ds['bareearth_fraction'].plot(ax=ax[3,0],cmap='Oranges',vmin=0, vmax=1,cbar_kwargs={'label':''})

    ax[0,0].set_title('Total built fraction',fontsize=11)
    ax[1,0].set_title('Building plan fraction',fontsize=11)
    ax[2,0].set_title('Road and paths fraction',fontsize=11)
    ax[3,0].set_title('Bare earth fraction',fontsize=11)

    ds['total_pervious'].plot(ax=ax[0,1],cmap=cmap, vmin=0, vmax=1,cbar_kwargs={'label':''})
    ds['tree_fraction'].plot(ax=ax[1,1],cmap='Greens',vmin=0, vmax=1,cbar_kwargs={'label':''})
    (ds['grass_fraction']+ds['shrub_fraction']).plot(ax=ax[2,1],cmap='Greens',vmin=0, vmax=1,cbar_kwargs={'label':''})
    ds['water_fraction'].plot(ax=ax[3,1],cmap='Blues',vmin=0,vmax=1,cbar_kwargs={'label':''})

    ax[0,1].set_title('Total pervious fraction',fontsize=11)
    ax[1,1].set_title('Tree canopy fraction',fontsize=11)
    ax[2,1].set_title('Grass and shrub fraction',fontsize=11)
    ax[3,1].set_title('Water fraction',fontsize=11)

    lat = ds.latitude.max().values - 0.09
    lon = ds.longitude.min().values + 0.06

    for a in ax.flatten():
        a.set_xlabel('')
        a.set_ylabel('')
        a = distance_bar(lat,lon,a1=a)
        
    for i in [0,1,2,3]:
        ax[i,0].set_ylabel('latitude')

    for i in [0,1]:
        ax[3,i].set_xlabel('longitude')

    fig.subplots_adjust(wspace=0.07)
    
    # plt.show()
    fig.savefig(f'{plotpath}/{domain}_landcover_{res}_{version}.png',bbox_inches='tight',dpi=300)

    return

def plot_morphology(ds,domain,res):

    # cmap = 'plasma'

    cmap = create_cmap('viridis',1000,start_col='black')

    plt.close('all')
    fig,ax = plt.subplots(ncols=2,nrows=4,figsize=(8,12),sharex=True,sharey=True)

    ds['building_height'].plot(ax=ax[0,0],cmap=cmap, vmin=0,vmax=25,cbar_kwargs={'label':''})
    ds['building_height_max'].plot(ax=ax[1,0],cmap=cmap,vmin=0,vmax=25,cbar_kwargs={'label':''})
    ds['wall_density'].plot(ax=ax[2,0],cmap=cmap,vmin=0,vmax=1,cbar_kwargs={'label':''})
    ds['roughness_kanda'].plot(ax=ax[3,0],cmap=cmap,vmin=0,vmax=1,cbar_kwargs={'label':''})

    ax[0,0].set_title('Building height mean (m)',fontsize=11)
    ax[1,0].set_title('Building height max (m)',fontsize=11)
    ax[2,0].set_title('Wall area density (-)',fontsize=11)
    ax[3,0].set_title('Roughness length: Kanda (m)',fontsize=11)

    ds['tree_height'].plot(ax=ax[0,1],cmap=cmap, vmin=0, vmax=25,cbar_kwargs={'label':''})
    ds['skyview_factor'].plot(ax=ax[1,1],cmap=cmap,vmin=0, vmax=1,cbar_kwargs={'label':''})
    ds['height_to_width'].plot(ax=ax[2,1],cmap=cmap,vmin=0, vmax=1,cbar_kwargs={'label':''})
    ds['displacement_kanda'].plot(ax=ax[3,1],cmap=cmap,vmin=0, vmax=10,cbar_kwargs={'label':''})

    ax[0,1].set_title('Tree height mean (m)',fontsize=11)
    ax[1,1].set_title('Sky view factor: buildings (-)',fontsize=11)
    ax[2,1].set_title('Canyon height to width ratio (-)',fontsize=11)
    ax[3,1].set_title('Displacement height: Kanda (m)',fontsize=11)

    lat = ds.latitude.max().values - 0.09
    lon = ds.longitude.min().values + 0.06

    for a in ax.flatten():
        a.set_xlabel('')
        a.set_ylabel('')
        a = distance_bar(lat,lon,a1=a)
        
    for i in [0,1,2,3]:
        ax[i,0].set_ylabel('latitude')

    for i in [0,1]:
        ax[3,i].set_xlabel('longitude')

    fig.subplots_adjust(wspace=0.07)
    
    # plt.show()
    fig.savefig(f'{plotpath}/{domain}_morphology_{res}_{version}.png',bbox_inches='tight',dpi=300)

    return

def get_geoscape_data(domain,res,version):

    ds = xr.open_dataset(f'{datapath}/{domain.capitalize()}_full_surface_description_{res}_{version}.nc')
    ds['total_vegetation'] = ds['tree_fraction'] + ds['grass_fraction'] + ds['shrub_fraction']
    ds['total_vegetation'].values = np.where(ds['total_vegetation']>1,1,ds['total_vegetation'])
    ds['total_built'].values = np.where(ds['total_built']>1,1,ds['total_built'])
    
    return ds

def create_cmap(name,num,start_col='black'):

    from matplotlib.colors import LinearSegmentedColormap

    cmap = plt.get_cmap(name,num)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # cmaplist[0] = (1.0,1.0,1.0,1.0)
    if start_col=='black':
        cmaplist[0] = (0.,0.,0.,1.)
    elif start_col=='white':
        cmaplist[0] = (1.,1.,1.,1.)
    cmap = LinearSegmentedColormap.from_list('mcm',cmaplist, cmap.N)

    return cmap

def distance_bar(lat,lon,a1):

    # plot distance bar
    start = (lon,lat)
    end = cgeo.Geodesic().direct(points=start,azimuths=90,distances=20000).base[:,0:2][0]
    a1.plot([start[0],end[0]],[start[1],end[1]], color='k', linewidth=1.5, mew=1)
    a1.text(start[0]-0.002,start[1]+0.005, '20 km', color='black', 
        fontsize=8, ha='left',va='bottom')

    return a1

################################################################################

if __name__ == "__main__":
    main()