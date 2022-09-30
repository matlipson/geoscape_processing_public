__title__ = "Regrid Geoscape onto template"
__version__ = "2022-05-25"
__author__ = "Mathew Lipson"
__email__ = "mathew.lipson@unsw.edu.au"

'''
Associated with the manuscript: A transformation in city-descriptive input data for urban climate model
Developed using Buildings © Geoscape Australia 2020: https://geoscape.com.au/legal/data-copyright-and-disclaimer/

Copyright 2022 Mathew Lipson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import os
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from rioxarray import merge

import shapely
import geopandas as gpd
import cartopy.geodesic as cgeo

import warnings
import matplotlib.pyplot as plt

oshome=os.getenv('HOME')

'''
Instructions: 
1. Install necessary python dependancies (see environment.yml)
2. Get sample data from Geoscape: https://geoscape.com.au/get-sample/ and extract into  ./data folder
3. Provide template file for grid or use supplied 300 m grid from 'CCI' (based on http://maps.elie.ucl.ac.be/CCI/viewer/index.php)
4. Update project and data paths if required
5. Run script

Note: using "do_shape_to_raster" and "replace_height_with_raster" will improve accuracy of height statistics by replacing
the centroid calculation (which applies statistics to a single grid only) with a 2 m resolution rasterised height map 
(for which parts of buildings can be applied to multiple grids). HOWEVER, processing takes orders of magnitude longer.
'''

################################################################################
# user input

domain = 'sample' # use 'sample' for data available at https://geoscape.com.au/get-sample/
grid   = 'CCI'
version = 'v1.02'
projpath = '.'

# if first run set to True
do_buildings = True  # process Geoscape Buildings shapefile into a GeoPackage
do_landcover = True  # process Geoscape Land Cover to a template grid
do_morphology = True # process buildings GeoPackage, landcover and Trees data into final dataset

# optional: WARNING: processing takes orders of magnitude longer.
do_shape_to_raster = False # burn vector to raster for building heights (long to process, uses GeoCube python package)
replace_height_with_raster = False # use pre-processed raster for building heights (accounts for buildings over multiple grids)
with_30m = False     # use 30 m land cover data as well as 2 m (mismatched)

missing_float = -999.

################################################################################
# v0.92 : for peer review in Frontiers in Environmental Science
# v1.0  : for publication in Frontiers in Environmental Science
# v1.01 : lowveg_fraction seperated into grass_fraction and shrub_fraction
# v1.02 : more accurate building perimeter calculation using geodesic distance

################################################################################

def main():

    for mkpath in [plotpath,outpath]:
        if not os.path.exists(mkpath):
            os.makedirs(mkpath)
            print("created folder : ", mkpath)

    if do_buildings:
        buildings = main_create_buildings_geopackage()
    if do_landcover:
        derived = main_resample_cover(grid=grid)
    if do_morphology:
        complete = main_calculate_morphology(grid=grid)

    return complete

def main_create_buildings_geopackage():
    '''
    create geopandas geopackage of geoscape buildings after calculating necessary building-level data
    saves as geopackage
    '''

    ################################################################################
    # format individual building information

    print('reading Geoscape building shapefile')
    orig = gpd.read_file(shp_fpath)

    raw = orig.copy()
    raw.columns = orig.columns.str.lower()

    if domain == 'melbourne':
        # drop buildings without height info
        raw = raw.rename(columns={'vic_heig_4':'bld_hgt'})
        raw = raw[raw.bld_hgt.notna()]
    else:
        print('assuming building height is average of roof and eave height')
        raw['bld_hgt'] = (raw['roof_hgt']+raw['eave_hgt'])/2
        raw['bld_max_hgt'] = raw['roof_hgt']
    
    buildings = raw[['bld_hgt','bld_max_hgt','area','geometry']]

    print('calculating each building perimeter')
    buildings['perimeter'] = geodesic_area_and_perimeter(buildings,mode='perim')
    # buildings['area'] = geodesic_area_and_perimeter(buildings,mode='area')

    print('calculating wall area per building')
    buildings['wall_area'] = buildings['perimeter']*buildings['bld_hgt']

    print('calculating frontal area per building')
    x_frontal = buildings.apply(get_dist_x,axis=1)*buildings['bld_hgt']
    y_frontal = buildings.apply(get_dist_y,axis=1)*buildings['bld_hgt']
    buildings['avg_frontal'] = (x_frontal + y_frontal)/2

    print('saving geometry as geopackage')
    buildings.to_file(geopkg_outpath, driver='GPKG')

    return buildings

################################################################################

def main_resample_cover(grid):
    '''main function for resampling geoscape surface cover data to a lower resolution and recategorise'''

    ################################################################################
    # opening geoscape raster

    print(f'opening {domain} geotif')

    # read tif
    tif_list = []
    for tif_fpath in tif_fpaths:
        tif_list.append(rxr.open_rasterio(tif_fpath,masked=True))
    print('merging tifs')
    tif_orig =  merge.merge_arrays(tif_list)

    # drop band 1
    tif_orig = tif_orig.drop(['band']).squeeze()

    print('preprojecting geotif to lat/lon')
    geo_orig = tif_orig.rio.reproject('epsg:4326')

    # reconfigure 
    geo_orig = geo_orig.rename({'y':'latitude','x':'longitude'})

    ################################################################################
    # rescaling geoscape raster

    print(f'loading {grid} for template')
    template_full = xr.open_dataset(template_fpath)

    # use buildings to define extent
    xmin,ymin,xmax,ymax = get_bld_extent(shp_fpath,template_full)

    # use land cover to define extent
    # xmin,ymin,xmax,ymax = geo_orig.rio.bounds()

    try: 
        ykey,xkey = 'latitude','longitude'
        template = template_full.sel({'latitude':slice(ymin,ymax),'longitude':slice(xmin,xmax)})
    except Exception:
        print('trying different coordinate name')
        ykey,xkey = 'lat','lon'
        template = template_full.sel({'lat':slice(ymin,ymax),'lon':slice(xmin,xmax)}).rename({ykey:'latitude',xkey:'longitude'})

    if len(template['latitude'])==0:  # check if any data
        print('no data found with those extents, trying reverse lat bounds')
        template = template_full.sel({ykey:slice(ymax,ymin),xkey:slice(xmin,xmax)}).rename({ykey:'latitude',xkey:'longitude'})

    template = template[['latitude','longitude']]

    print(f'creating resampled dataset with {grid} grid')
    geo_frac = create_fractions(template,geo_orig)
    geo_frac = set_geoscape_names(geo_frac)
    geo_frac = set_attributes(geo_frac,grid=grid)

    ################################################################################

    print('saving derived surface cover')
    derived = calc_derived_fractions(geo_frac)
    derived = remove_partial_sums(derived)
    derived = set_attributes(derived,grid=grid)

    write_netcdf(derived,geoscape_outpath,tif=False)

    ################################################################################

    assert_tests(geo_frac,derived)

    return derived

################################################################################

def main_calculate_morphology(grid):

    '''main function for converting building geopackage and tree data into final dataset'''

    print('reading Geoscape tree raster')
    tree_orig = rxr.open_rasterio(tree_fpath,masked=True).squeeze(['band'],drop=True)
    tree_orig = tree_orig.where(tree_orig>0)

    print('reprojecting Geoscape trees')
    tree_orig.rio.write_nodata(np.nan,inplace=True)
    tree_ds = tree_orig.rio.reproject('epsg:4326').rename({'y':'latitude','x':'longitude'})

    print('loading template for regridding')
    template = xr.open_dataset(geoscape_outpath)

    print('regridding tree data')
    tree_gridded = regrid_raster_for_hgt(tree_ds,template,'tree')

    write_netcdf(tree_gridded,f'{outpath}/{domain.capitalize()}_tree_gridded_{grid}.nc',tif=False,dtype='float32')

    print('calculating building gridded values from template')
    buildings = gpd.read_file(geopkg_outpath)
    building_gridded,cell = calc_all_grids_from_template(buildings,template)

    # burn height shapefile to netcdf (one-off)
    if do_shape_to_raster:
        '''This takes a longer time to process and requires additional modules (GeoCube)'''
        rasterise_height_shapefile(buildings,tree_ds,key='bld_hgt')
        rasterise_height_shapefile(buildings,tree_ds,key='bld_max_hgt')

    if replace_height_with_raster:
        ''' vector point heights have the problem that large buildings can not be split across grids.
        Converting vector (shapefile) data to raster for building heights is an advantage because
        larger buildings can affect multiple grids. This method is not appropriate for wall statistics,
        as it is not clear how to split a wall or height-to-width statistics between grids'''
        
        fname_method = 'byraster'

        print('replacing vector centroid point heights with raster heights')

        full_raster_hgt_fpath = f'{raster_hgt_fpath}_bld_hgt.nc'
        bld_hgt1 = xr.open_dataset(full_raster_hgt_fpath)
        bld_gridded1 = regrid_raster_for_hgt(bld_hgt1['bld_hgt'],template,'building')
        building_gridded.update(bld_gridded1[['building_height','building_height_std']])

        print('replacing vector centroid point heights with raster heights max')
        
        full_raster_hgt_fpath = f'{raster_hgt_fpath}_bld_max_hgt.nc'
        bld_hgt2 = xr.open_dataset(full_raster_hgt_fpath)
        bld_gridded2 = regrid_raster_for_hgt(bld_hgt2['bld_max_hgt'],template,'building')
        building_gridded.update(bld_gridded2[['building_height_max']])

    else:
        fname_method = 'bypoint'
        

    ds = xr.merge([building_gridded,tree_gridded])

    ds['building_fraction']  = template['building_fraction']
    ds['tree_fraction']      = template['tree_fraction']
    ds['grass_fraction']     = template['grass_fraction']
    ds['shrub_fraction']     = template['shrub_fraction']
    ds['water_fraction']     = template['water_fraction']
    ds['bareearth_fraction'] = template['bareearth_fraction']
    ds['roadpath_fraction']  = template['roadpath_fraction']
    ds['total_built']        = template['total_built']
    ds['total_pervious']     = template['total_pervious']
    ds['total_vegetation']   = template['total_vegetation']

    print('calculating height to width ratio')
    # from masson 2020: h/w = lambda_w/(2*(1-lambda_p))
    ds['height_to_width'] = ds['wall_density']/(2*(1.-ds['building_fraction'])).round(4)

    # replace undefined with nan (where building fraction = 1)
    ds['height_to_width'].values = np.where(ds['height_to_width'].values==np.inf,np.nan,ds['height_to_width'].values)

    print('calculating skyview factor')
    # from masson 2020: [(h/w)^2 + 1]^(1/2) - h/w
    ds['skyview_factor'] = ((ds['height_to_width']**2 + 1.)**(1./2) - ds['height_to_width']).round(4)

    ds['displacement_mac'], ds['roughness_mac'] = calc_site_roughness(
        bldhgt_ave=ds['building_height'],
        sigma_h=ds['building_height_std'],
        lambda_p=ds['building_fraction'],
        lambda_pv=ds['tree_fraction'],
        lambda_f=ds['frontal_density'],
        bldhgt_max=ds['building_height_max'],
        mode=0) # macdonald et al 1998 staggered grid

    ds['displacement_kanda'], ds['roughness_kanda'] = calc_site_roughness(
        bldhgt_ave=ds['building_height'],
        sigma_h=ds['building_height_std'],
        lambda_p=ds['building_fraction'],
        lambda_pv=ds['tree_fraction'],
        lambda_f=ds['frontal_density'],
        bldhgt_max=ds['building_height_max'],
        mode=3) # Kanda et al 2013

    ds['displacement_mac'] = ds['displacement_mac'].round(4)
    ds['roughness_mac'] = ds['roughness_mac'].round(4)
    ds['displacement_mac'] = ds['displacement_mac'].round(4)
    ds['roughness_mac'] = ds['roughness_mac'].round(4)
    ds['displacement_kanda'] = ds['displacement_kanda'].round(4)
    ds['roughness_kanda'] = ds['roughness_kanda'].round(4)

    ds = set_attributes(ds,grid=grid)

    ################################################################################

    complete_outpath = f'{outpath}/{domain.capitalize()}_full_surface_description_{fname_method}_{grid}_{version}.nc'

    complete_outpath = f'{outpath}/{domain.capitalize()}_full_surface_description_{grid}_{version}.nc'

    write_netcdf(ds,complete_outpath,tif=True,zlib=True,dtype='float32')

    return ds

################################################################################

def get_bld_extent(shp_fpath,template_full,border_pxl=1):
    '''
    calculates the extent of the template domain based on building shapefile 
    '''

    print('reading Geoscape building shapefile')
    bld = gpd.read_file(shp_fpath)

    print('subsetting to buildings extent')
    xmin,ymin,xmax,ymax = bld.total_bounds
    # add one pixel border around buildings extent to capture all
    xmin2 = xmin-border_pxl*abs(template_full.rio.resolution()[0])
    xmax2 = xmax+border_pxl*abs(template_full.rio.resolution()[0])
    ymin2 = ymin-border_pxl*abs(template_full.rio.resolution()[1])
    ymax2 = ymax+border_pxl*abs(template_full.rio.resolution()[1])

    return xmin2,ymin2,xmax2,ymax2

def create_fractions(template,hires,xkey='longitude',ykey='latitude'):
    '''
    Description:
        Calculates the sub-tile surface fractions from big and small tiled data, 
        ensuring grid edges align between large and small datasets
    Arguments:
        template: xarray dataset, land cover fraction dataset with big tiles
        hires: xarray dataset, land cover fraction dataset with small tiles

    Returns:
        frac_ds: xarray dataset, fractional land cover at big resolution
    '''

    xbounds,ybounds = f'{xkey}_bounds',f'{ykey}_bounds'

    # match extents
    s1,s2,s3,s4 = hires.rio.bounds()
    b1,b2,b3,b4 = template.rio.bounds()
    big = template.sel({ykey:slice(s2,s4),xkey:slice(s1,s3)})
    small = hires.sel({ykey:slice(b2,b4),xkey:slice(b1,b3)})

    # check if problematic latitude order
    if len(small.latitude) == 0:
        small = hires.sel({ykey:slice(b4,b2),xkey:slice(b1,b3)})
    if len(big.latitude) == 0:
        big = template.sel({ykey:slice(s4,s2),xkey:slice(b1,b3)})

    # add cell bounds
    big = add_bounds(big)

    frac_list = []

    for lat in big[ykey].values:
        for lon in big[xkey].values:

            big_cell = big.sel({ykey:lat,xkey:lon})

            ybound = big_cell[ybounds].values
            xbound = big_cell[xbounds].values

            # get cells bounded by larger cell
            small_cell = small.sel({ykey:slice(ybound[0],ybound[1]),xkey:slice(xbound[0],xbound[1])})
            if len(small_cell.latitude) == 0:
                small_cell = small.sel({ykey:slice(ybound[1],ybound[0]),xkey:slice(xbound[0],xbound[1])})
            # if float(small_cell.count()) == 0:
            #     continue

            counts = small_cell.to_series().value_counts()
            counts.index = counts.index.map(str)

            frac = pd.Series({ykey:lat,xkey:lon})
            frac = frac.append(counts/counts.sum())

            frac_list.append(frac)

    df = pd.DataFrame(frac_list)

    # order columns
    cols = [ykey,xkey] + sorted(df.drop(columns=[ykey,xkey]).columns)
    df = df[cols]

    frac_ds = df.set_index([ykey,xkey]).to_xarray()

    # latitude go from 90 to -90
    frac_ds = frac_ds.sortby('latitude',ascending=False)

    return frac_ds

def set_geoscape_names(ds):

    # set geoscape names
    all_keys = {
        '2.0': 'bare_earth',
        '3.0': 'road_path',
        '4.0': 'grass',
        '5.0': 'trees',
        '6.0': 'other_veg',
        '7.0': 'built_area',
        '8.0': 'water',
        '9.0': 'buildings',
        '10.0':'cloud',
        '11.0':'shadow',
        '12.0':'swimming_pool'
        }

    inc_keys = {key: all_keys[key] for key in [k for k in ds.keys()]}

    ds = ds.rename_vars(inc_keys)

    return ds

def regrid_raster_for_hgt(da,template,varkey,xkey='longitude',ykey='latitude'):
    '''
    Description:
        Calculates tile height statistics from small tiled data to big template,
        ensuring grid edges align between large and small datasets
    Arguments:
        da: xarray data array to be regridded
        template: xarray dataset, dataset with big tiles

    Returns:
        ds: xarray dataset, gridded height stats at template grid
    '''

    # add cell bounds
    big = add_bounds(template)

    xbounds,ybounds = f'{xkey}_bounds',f'{ykey}_bounds'

    # select small extent to edge of big
    ymin,ymax = big[ybounds].values[0,0], big[ybounds].values[1,-1]
    xmin,xmax = big[xbounds].values[0,0], big[xbounds].values[1,-1]
    
    small = da.sel({xkey:slice(xmin,xmax), ykey:slice(ymin,ymax)})
    if len(small.values.flatten()) == 0:
        small = da.sel({xkey:slice(xmin,xmax), ykey:slice(ymax,ymin)})
        if len(small.values.flatten()) == 0:
            small = da.sel({xkey:slice(xmax,xmin), ykey:slice(ymin,ymax)})

    ser_list = []

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')

        for lat in big[ykey].values:
            for lon in big[xkey].values:

                # big_cell = big.sel(lat=lat,lon=lon)
                big_cell = big.sel({ykey:lat,xkey:lon})

                ybound = big_cell[ybounds].values
                xbound = big_cell[xbounds].values

                # get cells bounded by larger cell
                small_cell = small.sel({xkey:slice(xbound.min(),xbound.max()), ykey:slice(ybound.max(),ybound.min())})
                if len(small_cell.values.flatten()) == 0:
                    small_cell = small.sel({xkey:slice(xbound.min(),xbound.max()), ykey:slice(ybound.min(),ybound.max())})
                    if len(small_cell.values.flatten()) == 0:
                        small_cell = small.sel({xkey:slice(xbound.max(),xbound.min()), ykey:slice(ybound.max(),ybound.min())})

                if len(small_cell.values.flatten()) == 0:
                    ser = pd.Series({
                        ykey: lat,
                        xkey: lon,
                        f'{varkey}_height'    :np.nan,
                        f'{varkey}_height_max':np.nan,
                        f'{varkey}_height_std':np.nan
                        })
                else:
                    ser = pd.Series({
                        ykey: lat,
                        xkey: lon,
                        f'{varkey}_height'    :float(small_cell.mean()),
                        f'{varkey}_height_max':float(small_cell.max()),
                        f'{varkey}_height_std':float(small_cell.std())
                        })

                ser_list.append(ser)

    df = pd.DataFrame(ser_list)

    ds = df.set_index([ykey,xkey]).to_xarray()

    # latitude go from 90 to -90
    ds = ds.sortby('latitude',ascending=False)

    return ds

def write_netcdf(ds,fpath,tif=True,crs='epsg:4326',dtype='float32',zlib=False):

    # ensure output directory exists
    outpath = os.path.dirname(fpath)
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # remove wayward attributes
    for key in ds.data_vars:
        try:
            del ds[key].attrs['_FillValue']
            print('removed _FillValue')
        except Exception:
            pass
        try:
            del ds[key].attrs['grid_mapping']
            print('removed grid_mapping')
        except Exception:
            pass

    # simplify spatial_ref for netcdf
    ds.rio.write_crs(crs, inplace=True)
    ds.spatial_ref.attrs = {'grid_mapping_name': 'latitude_longitude'}

    for key in ds.data_vars:
        ds[key].encoding.update({'zlib':zlib, '_FillValue':missing_float, 'dtype':dtype})

    ds.to_netcdf(fpath,format='NETCDF4')

    if tif:
        ds2 = xr.open_dataset(fpath)
        ds2.rio.write_crs(crs, inplace=True)
        fpath2 = fpath.split('.nc')[0] + '.tiff'
        try:
            ds2.rio.to_raster(fpath2, compress='LZW')
        except Exception as e:
            print('raster tiff save error')
            print(e)

    return

################################################################################

def create_grid_from_nc_template(template,ykey='latitude',xkey='longitude'):

    '''calculate shapely polygon for each grid cell using netcdf template'''

    tmp = add_bounds(template)

    # get lat/lons of template grid bounds
    xbound,ybound = f'{xkey}_bounds', f'{ykey}_bounds'

    aa = np.meshgrid( tmp[xbound][0].values, tmp[ybound][0].values )[0].flatten()
    bb = np.meshgrid( tmp[xbound][0].values, tmp[ybound][0].values )[1].flatten()
    cc = np.meshgrid( tmp[xbound][1].values, tmp[ybound][1].values )[0].flatten()
    dd = np.meshgrid( tmp[xbound][1].values, tmp[ybound][1].values )[1].flatten()

    grid_cells = [shapely.geometry.box(a,b,c,d) for a,b,c,d in zip(aa,bb,cc,dd)]

    # get lat/lons of template grid centres
    lons,lats = np.meshgrid( tmp[xkey].values, tmp[ykey].values )
    lons,lats = lons.flatten(),lats.flatten()

    df = pd.DataFrame({
            'latitude':lats,
            'longitude':lons,
            'geometry':grid_cells
            })

    cell = gpd.GeoDataFrame(df,crs='epsg:4326')

    # calculate cell area from geodesic distances
    cell['cell_area'] = cell.apply(get_dist_x,axis=1)*cell.apply(get_dist_y,axis=1)

    return cell

def get_dist_x(row):
    '''building frontal wall area in x direction'''

    bounds = row['geometry'].bounds

    spoint = (bounds[0],bounds[1])
    epoint = (bounds[2],bounds[1])

    try:
        x = cgeo.Geodesic().inverse(spoint,epoint).base[0][0]  # old version
    except Exception:
        x = cgeo.Geodesic().inverse(spoint,epoint)[0][0]

    return x

def get_dist_y(row):
    '''building frontal wall area in y direction'''

    bounds = row['geometry'].bounds

    spoint = (bounds[0],bounds[1])
    epoint = (bounds[0],bounds[3])

    try: 
        y = cgeo.Geodesic().inverse(spoint,epoint).base[0][0]   # old version
    except Exception:
        y = cgeo.Geodesic().inverse(spoint,epoint)[0][0]

    return y

def get_midpoint(gdf):

    xmin, ymin, xmax, ymax= gdf.total_bounds
    xmid = (xmin+xmax)/2
    ymid = (ymin+ymax)/2

    return xmid,ymid

def get_perimeter(geom):

    # convert shapely polygon from 3D to 2D (ignore height)
    geom_2d = shapely.wkb.loads(shapely.wkb.dumps(geom, output_dimension=2))
        
    # calculate perimeter along the geodesic Earth
    perim = cgeo.Geodesic().geometry_length(geom_2d)

    return perim

def geodesic_area_and_perimeter(geodf,mode='perim'):
    '''calculate area and perimeter of shapely polygon using geodesic distance
    See https://gis.stackexchange.com/questions/413349/calculating-area-of-lat-lon-polygons-without-transformation-using-geopandas
    See https://pyproj4.github.io/pyproj/stable/api/geod.html'''
    if not geodf.crs and geodf.crs.is_geographic:
        raise TypeError('geodataframe should have geographic coordinate system')

    geod = geodf.crs.get_geod()

    # not needed
    def calc_perim(geom):
        '''use geopandas get_geod to calculate geometry length'''
        perim = geodf.crs.get_geod().geometry_length(geom)
        return perim

    def calc_area(geom):
        '''use geopandas get_geod to calculate geometry perimeter area'''
        area, perim = geodf.crs.get_geod().geometry_area_perimeter(geom)
        return abs(area)

    if mode=='perim':
        print('returning geodesic perimeter (m)')
        result = geodf.geometry.apply(calc_perim)
    if mode=='area':
        print('returning geodesic area (m2)')
        result = geodf.geometry.apply(calc_area)
    
    return result

def calc_all_grids_from_template(buildings,template):

    '''calculates grid-level morphology from building data using template grid '''

    print('creating grid from template')

    # new method accounts for uneven grids (as in CCI)
    cell = create_grid_from_nc_template(template)

    cell = cell.sort_values(['longitude','latitude'], ascending=[True,False]).reset_index(drop=True)

    # ## plot ###
    # plt.close('all')
    # fig, ax = plt.subplots(figsize=(20,10))
    # buildings.plot(ax=ax)
    # template['building_fraction'].plot(ax=ax,edgecolor='k',linewidth=0.3,zorder=-1)
    # # cell_orig.plot(ax=ax,column='cell_area',edgecolor='r',linewidth=0.3,alpha=0.5)
    # cell.plot(ax=ax,column='cell_area',edgecolor='cyan',linewidth=0.1,alpha=0.5)
    # plt.show()
    # # plt.savefig('figures/grid_new.pdf',dpi=600)

    print('calculate building centroids')
    # NOTE: Warning regarding centroid being incorrect can be ignored as the small scales being processed will have negligible error
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        geometry = gpd.points_from_xy(buildings.geometry.centroid.x, buildings.geometry.centroid.y)

    df = pd.DataFrame(buildings.drop(columns='geometry'))
    points = gpd.GeoDataFrame(df, crs='epsg:4326', geometry=geometry)
    # points2 = gpd.GeoDataFrame(buildings, crs='epsg:4283', geometry=geometry)

    print(f'calculating all grids using template')
    # calculate building height weighted by buidling area
    points['bld_hgt_wgt'] = points['bld_hgt'] * points['area']

    # spatial join building points within grid cells
    building_points = gpd.sjoin(points, cell, how='left', op='within')
    grp = building_points.groupby('index_right')


    gridded = pd.DataFrame(index=cell.index,data=cell[['latitude','longitude','cell_area']])

    print('calculating building heights')
    gridded['building_height'] = grp['bld_hgt_wgt'].sum()/grp['area'].sum() # use weighted calculation so larger buildings in grid have greater effect
    gridded['building_height_max'] = grp['bld_max_hgt'].max()
    gridded['building_height_std'] = grp['bld_hgt'].std()

    gridded['wall_density'] = grp['wall_area'].sum()/cell['cell_area']
    gridded['frontal_density'] = grp['avg_frontal'].sum()/cell['cell_area']

    gridded = gridded.set_index(['latitude','longitude'])

    ds = gridded.to_xarray()

    # latitude go from 90 to -90
    ds = ds.sortby('latitude',ascending=False)

    return ds,cell

def add_bounds(ds):
    """
    Returns a new dataset with bounds variables. The bounds values are guessed assuming
    equal spacing on either side of a coordinate label.
    derived from cf-xarray: https://github.com/xarray-contrib/cf-xarray/blob/533ed5c5ceefbea94a01dee3d6e1b62ab26e0866/cf_xarray/accessor.py
    """

    obj = ds.copy(deep=True) # don't affect original ds
    dimensions = [key for key in ds.dims.keys()]

    for dim in dimensions:
        if dim == 'time':
            continue
        if '_bounds' in dim:
            continue
        bname = f"{dim}_bounds"
        if bname in obj.variables:
            raise ValueError(f"Bounds variable name {bname!r} will conflict!")
        obj.coords[bname] = _guess_bounds_dim(obj[dim].reset_coords(drop=True))
        obj[dim].attrs["bounds"] = bname

    return obj

def _guess_bounds_dim(da):
    """
    Guess bounds values given a 1D coordinate variable.
    Assumes equal spacing on either side of the coordinate label.
    from cf-xarray: https://github.com/xarray-contrib/cf-xarray/blob/533ed5c5ceefbea94a01dee3d6e1b62ab26e0866/cf_xarray/accessor.py
        """
    assert da.ndim == 1

    dim = da.dims[0]
    diff = da.diff(dim)
    lower = da - diff / 2
    upper = da + diff / 2
    bounds = xr.concat([lower, upper], dim="bounds")

    first = (bounds.isel({dim: 0}) - diff[0]).assign_coords({dim: da[dim][0]})
    result = xr.concat([first, bounds], dim=dim)

    return result

def set_attributes(ds_orig,grid):

    # don't overight original

    ds = ds_orig.copy(deep=True)

    '''set dataset attributes '''

    ### coordinates ###
    for key in ['latitude','lat']:
        if key in ds.keys():
            ds[key].attrs = {
                    'long_name'    : 'latitude',
                    'units'        : 'degrees_north',
                    'standard_name': 'latitude',
                    'long_name'    : 'latitude',
                    'axis'         : 'Y',
                }

    for key in ['longitude','lon']:
        if key in ds.keys():
            ds[key].attrs = {
                    'long_name'    : 'longitude',
                    'units'        : 'degrees_east',
                    'standard_name': 'longitude',
                    'long_name'    : 'longitude',
                    'axis'         : 'X',
                }

    ### geoscape names

    key = 'bare_earth'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Bare earth',
                'description'   : 'Includes sand dunes, desert, rock outcrops, bare soil other than bare agricultural land, and sparsely vegetated areas of grass and shrub. Non-vegetated strip mines and quarries except where covered by development or water.',
                'units'         : '-',
        }
    key = 'road_path'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Road and path',
                'description'   : 'Roads and parking lots covered in a man-made material excluding hard packed dirt trails.',
                'units'         : '-',
        }
    key = 'grass'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Grass',
                'description'   : 'Grass and herbaceous areas. The category may include herbaceous wetlands if images are collected during dry season or periods of drought.',
                'units'         : '-',
        }
    key = 'trees'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Trees',
                'description'   : 'All trees including deciduous and evergreen woody vegetation.',
                'units'         : '-',
        }
    key = 'other_veg'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Unspecified Vegetation',
                'description'   : 'Any other vegetative material not included within the Grass or Tree class. This may include, but is not limited to, shrub, scrub, agriculture, and aquatic plants.',
                'units'         : '-',
        }
    key = 'built_area'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Built-up Areas',
                'description'   : 'Any areas of man-made environments and infrastructure excluding road and paths and buildings.',
                'units'         : '-',
        }
    key = 'water'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Water',
                'description'   : 'Depending on the resolution quality of the imagery used, natural water will include streams, canals, ponds, lakes, reservoirs, estuaries and bays.',
                'units'         : '-',
        }
    key = 'buildings'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Buildings',
                'description'   : 'Where the majority of a pixel intersects a Building, vector building polygon representation.',
                'units'         : '-',
        }
    key = 'cloud'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Cloud',
                'description'   : 'The area covered with cloud on Date of collection.',
                'units'         : '-',
        }
    key = 'shadow'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Shadow',
                'description'   : 'The area covered with shadow on Date/time of collection.',
                'units'         : '-',
        }
    key = 'swimming_pool'
    if key in ds.keys():
        ds[key].attrs = {
                'long_name'     : 'Swimming Pool',
                'description'   : 'An area identified as a swimming pool.',
                'units'         : '-',
        }


    ### surface cover ###
    if 'building_fraction' in ds.keys():
        ds['building_fraction'].attrs = {
                'long_name'    : 'Building footprint fraction (lambda_p)',
                'description'  : 'building footprint area as fraction of grid area, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'tree_fraction' in ds.keys():
        ds['tree_fraction'].attrs = {
                'long_name'    : 'Tree canopy fraction (lambda_vt)',
                'description'  : 'tree canopy plan area as fraction of grid area, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'lowveg_fraction' in ds.keys():
        ds['lowveg_fraction'].attrs = {
                'long_name'    : 'Grass and low vegetation fraction (lambda_vl)',
                'description'  : 'low vegetation (grass, shrubs, other vegetation) as fraction of grid area, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'grass_fraction' in ds.keys():
        ds['grass_fraction'].attrs = {
                'long_name'    : 'Grass fraction (lambda_vg)',
                'description'  : 'grass as fraction of grid area, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'shrub_fraction' in ds.keys():
        ds['shrub_fraction'].attrs = {
                'long_name'    : 'Shrub fraction (lambda_vs)',
                'description'  : 'shrub and other vegetetaion as fraction of grid area, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'water_fraction' in ds.keys():
        ds['water_fraction'].attrs = {
                'long_name'    : 'Open water fraction (lambda_wa)',
                'description'  : 'all open water (ocean, lakes, pools) as fraction of grid area, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'bareearth_fraction' in ds.keys():
        ds['bareearth_fraction'].attrs = {
                'long_name'    : 'Bare earth, sand and rock fraction (lambda_be)',
                'description'  : 'bare earth including construction sites, rock, sand and sparsely vegetated areas as fraction of grid area, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'roadpath_fraction' in ds.keys():
        ds['roadpath_fraction'].attrs = {
                'long_name'    : 'Road, path and other hard surfaces on ground (lambda_if)',
                'description'  : 'all hard surfaces on ground excluding buildings, defined as "impervious surface fraction" in Stewart and Oke, 2012, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'total_built' in ds.keys():
        ds['total_built'].attrs = {
                'long_name'    : 'Total built fraction (lambda_tb)',
                'description'  : 'all impervious surfaces including buildings, roads, paths and other hard surfaces, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'total_pervious' in ds.keys():
        ds['total_pervious'].attrs = {
                'long_name'    : 'Total pervious fraction (lambda_tp)',
                'description'  : 'all pervious surfaces including vegetation, water and bare earth, corrected for cloud and shadow fractions',
                'units'        : '-',
            }
    if 'total_vegetation' in ds.keys():
        ds['total_vegetation'].attrs = {
                'long_name'    : 'Total vegetation fraction (lambda_v)',
                'description'  : 'all vegetation, corrected for cloud and shadow fractions',
                'units'        : '-',
            }

    # morphology
    if 'frontal_density' in ds.keys():
        ds['frontal_density'].attrs = {
                'long_name'    : 'Frontal area density (lambda_f)',
                'description'  : 'sum of cardinally averaged building frontal area as fraction of grid area',
                'units'        : '-',
            }
    if 'wall_density' in ds.keys():
        ds['wall_density'].attrs = {
                'long_name'    : 'Wall area density (lambda_w)',
                'description'  : 'sum of building wall area as fraction of grid area',
                'units'        : '-',
            }

    if 'building_height' in ds.keys():
        ds['building_height'].attrs = {
                'long_name'    : 'Building height mean (H_avg)',
                'description'  : 'mean building height in grid cell (avg. of geoscape roof and eave height)',
                'units'        : 'm',
            }
    if 'building_height_max' in ds.keys():
        ds['building_height_max'].attrs = {
                'long_name'    : 'Building height maximum (H_max)',
                'description'  : 'maximum building height in grid cell (avg. of geoscape roof and eave height)',
                'units'        : 'm',
            }
    if 'building_height_std' in ds.keys():
        ds['building_height_std'].attrs = {
                'long_name'    : 'Building height standard deviation (H_std)',
                'description'  : 'standard deviation of building height in grid cell (avg. of geoscape roof and eave height)',
                'units'        : 'm',
            }

    if 'tree_height' in ds.keys():
        ds['tree_height'].attrs = {
                'long_name'    : 'Tree height average (HT_ave)',
                'description'  : 'average vegetation canopy height in grid',
                'units'        : 'm',
            }
    if 'tree_height_std' in ds.keys():
        ds['tree_height_std'].attrs = {
                'long_name'    : 'Tree height standard deviation (HT_ave)',
                'description'  : 'standard deviation in vegetation canopy height in grid',
                'units'        : 'm',
            }

    if 'roughness_mac' in ds.keys():
        ds['roughness_mac'].attrs = {
                'long_name'    : 'Aerodynamic roughness length (Z_0,mac)',
                'description'  : 'roughness length for staggered arrays from Macdonald et al., 1998: https://doi.org/10.1016/S1352-2310(97)00403-2',
                'units'        : 'm',
            }
    if 'displacement_mac' in ds.keys():
        ds['displacement_mac'].attrs = {
                'long_name'    : 'Displacement height (Z_d,mac)',
                'description'  : 'zero-plane displacement height from Macdonald et al., 1998: https://doi.org/10.1016/S1352-2310(97)00403-2',
                'units'        : 'm',
            }
    if 'roughness_kanda' in ds.keys():
        ds['roughness_kanda'].attrs = {
                'long_name'    : 'Aerodynamic roughness length (Z_0,kan)',
                'description'  : 'roughness length for staggered arrays from Kanda et al., 2013: https://doi.org/10.1007/s10546-013-9818-x',
                'units'        : 'm',
            }
    if 'displacement_kanda' in ds.keys():
        ds['displacement_kanda'].attrs = {
                'long_name'    : 'Displacement height (Z_d,kan)',
                'description'  : 'zero-plane displacement height from Kanda et al., 2013: https://doi.org/10.1007/s10546-013-9818-x',
                'units'        : 'm',
            }

    ### canyon attributes ###
    if 'height_to_width' in ds.keys():
        ds['height_to_width'].attrs = {
                'long_name'    : 'Canyon height to width aspect ratio (h/w)',
                'description'  : 'average aspect ratio assuming street canyon geometry using Eq 1 of Masson et al. 2020: https://doi.org/10.1016/j.uclim.2019.100536',
                'units'        : '-',
            }
    if 'skyview_factor' in ds.keys():
        ds['skyview_factor'].attrs = {
                'long_name'    : 'Skyview factor (psi)',
                'description'  : 'average skyview factor assuming street canyon geometry using Eq 2 of Masson et al. 2020: https://doi.org/10.1016/j.uclim.2019.100536',
                'units'        : '-',
            }

    # other
    if 'cell_area' in ds.keys():
        ds['cell_area'].attrs = {
                'long_name'    : 'Cell area',
                'description'  : f'area of grid cell, at {grid} resolution',
                'units'        : 'm2',
            }
    if 'bld_hgt' in ds.keys():
        ds['bld_hgt'].attrs = {
                'long_name'    : 'Individual building height',
                'description'  : 'average for midpoint of roof and eave height, per building',
                'units'        : 'm',
            }
    if 'bld_max_hgt' in ds.keys():
        ds['bld_max_hgt'].attrs = {
                'long_name'    : 'Individual maximum building height',
                'description'  : 'Maximum of roof height, per building',
                'units'        : 'm',
            }
    if 'bld_hgt_wgt' in ds.keys():
        ds['bld_hgt_wgt'].attrs = {
                'long_name'    : 'Building height weighted average',
                'description'  : 'area weighted average for midpoint of roof and eave height in grid',
                'units'        : 'm',
            }
    if 'roof_hgt' in ds.keys():
        ds['roof_hgt'].attrs = {
                'long_name'    : 'Roof height average',
                'description'  : 'average roof height in grid (no buiding area weighting)',
                'units'        : 'm',
            }
    if 'roof_hgt_wgt' in ds.keys():
        ds['roof_hgt_wgt'].attrs = {
                'long_name'    : 'Roof height weighted average',
                'description'  : 'area weighted average roof height in grid',
                'units'        : 'm',
            }
    if 'eave_hgt' in ds.keys():
        ds['eave_hgt'].attrs = {
                'long_name'    : 'Eave height average',
                'description'  : 'average eave height in grid (no buiding area weighting)',
                'units'        : 'm',
            }
    if 'eave_hgt_wgt' in ds.keys():
        ds['eave_hgt_wgt'].attrs = {
                'long_name'    : 'Eave height weighted average',
                'description'  : 'area weighted average eave height in grid',
                'units'        : 'm',
            }

    ds.attrs['title']       = f'Urban characteristics for {domain.capitalize()} at {grid} resolution derived from Geoscape datasets'
    ds.attrs['version']     = version
    ds.attrs['author']      = 'Mathew Lipson <m.lipson@unsw.edu.au>'
    ds.attrs['institution'] = 'ARC Centre of Excellence for Climate Extremes, UNSW Sydney, Australia'
    ds.attrs['source']      = 'Developed using Geoscape Buildings v2.0, Trees v1.6 and Surface cover v1.6 (c) Geoscape Australia 2020. https://geoscape.com.au/legal/data-copyright-and-disclaimer/'
    ds.attrs['licence']     = 'Data in this file is available under Creative Commons Attribution 4.0 International (CC-BY) with attribution: https://creativecommons.org/licenses/by/4.0/legalcode'
    ds.attrs['publication']   = 'A transformation in city-descriptive input data for urban climate models: Frontiers in Environmental Science 2022'
    ds.attrs['publication_authors'] = 'Mathew Lipson, Negin Nazarian, Melissa hart, Kerry Nice, Brooke Conroy'

    return ds

def calc_site_roughness(bldhgt_ave,sigma_h,lambda_p,lambda_pv,lambda_f,bldhgt_max=None,mode=0):
    '''estimate urban site roughness length for momentum and zero-plane displacement based on various methods 
    as described in Grimmond and Oke (1999): Aerodynamic Properties of Urban Areas Derived from Analysis of Surface Form .

    modes: four different morphometric methods to calculate roughness and displacement:
        - 0: Macdonald et al. 1998: An improved method for the estimation of surface roughness of obstacle arrays
        - 1: Kent et al 2017: Aerodynamic roughness parameters in cities: Inclusion of vegetation
        - 2: Millward-Hopkins et al., (2011) per Kent et al. 2017 eq 11 and 12 # TYPO IN EQ 12 OF KENT PAPER
        - 3: Kanda et al., 2013

    Inputs
    ------
    bldhgt_ave [m] : building mean height
    sigma_h [m]    : building height standard deviation
    lambda_p [0-1] : building plan area fraction
    lambda_pv [0-1]: tree plan area fraction
    lambda_f [1]   : wall frontal area fraction
    bldhgt_max [m] : maximum building height (if none assume 1.25 times SD)

    Returns
    -------
    zd [m] zero-plane displacement
    z0 [m] rougness length for momentum

    Warning
    -------
    in mode=0 wall frontal area (assuming random canyon orientation) calculated from Porsen et al 2010 (DOI:10.1002/qj.668)
    in mode=1 (including vegetation), vegetation frontal area index assumed to be equal to tree plan fraction (from Table 1 in Kent et al 2017)
    in mode=3 (Kanda) maximum building height assumed to be 2 standard deviations
    '''

    vonkarm = 0.4
    dragbld = 1.2

    if mode==0:
        print('zd,z0 from Macdonald et al. 1998')

        alpha   = 4.43 # for staggered arrays
        beta    = 1.0  # for staggered arrays

        zd = (1. + alpha**(-lambda_p)*(lambda_p - 1.))*bldhgt_ave    # eq 23 from Mac1998
        z0 = ((1. - zd/bldhgt_ave) * np.exp( -1.*(0.5*beta*(dragbld/vonkarm**2)*(1.-zd/bldhgt_ave)*lambda_f)**(-0.5) ))*bldhgt_ave  # eq 26 from Mac1998

    if mode == 1: # with veg from Kent et al (2018) https://doi.org/10.1007/s11252-017-0710-1
        print('zd,z0 from Macdonald with vegetation from Kent et al. 2017 and 2018')
        # Aerodynamic roughness variation with vegetation: analysis in a suburban neighbourhood and a city park

        alpha   = 4.43 # for staggered arrays
        beta    = 1.0  # for staggered arrays

        lambda_fv = lambda_pv # estimated equality based on Table 1 of Kent et al 2017 for different sites
        P_3D = 0.4   # leaf-on porosity

        lambda_tot = (lambda_p + lambda_pv*(1. - P_3D)) # eq 4 from Kent et al (2018)
        Pv  = (-1.251*P_3D**2 + 0.489*P_3D + 0.803)/dragbld  # eq 7 from Kent et al (2018)
        weighted_frontal_area = lambda_f + lambda_fv*Pv

        zd = (1. + alpha**(-lambda_tot)*(lambda_tot-1.))*bldhgt_ave # eq 5
        z0 = (1. - zd/bldhgt_ave)*np.exp( -1.*(0.5*beta*(dragbld/vonkarm**2)*(1.-zd/bldhgt_ave)*weighted_frontal_area)**(-0.5) )*bldhgt_ave # eq 6

    if mode == 2: # Millward-Hopkins et al. (2011), per Kent et al. 2017 eq 11 and 12
        print('zd,z0 from Millward-Hopkins et al. 2011, per Kent et al. 2017 eq 11 and 12')

        MhoUzd_on_Hav_A = (19.2*lambda_p - 1. + np.exp(-19.2*lambda_p))/(19.2*lambda_p*(1.- np.exp(-19.2*lambda_p)))
        MhoUzd_on_Hav_B = (117*lambda_p + (187.2*lambda_p**3 - 6.1)*(1.- np.exp(-19.2*lambda_p)) )/( (1.+114*lambda_p + 187*lambda_p**3)*(1.-np.exp(-19.2*lambda_p)) )
        MhoUzd_on_Hav = np.where( lambda_p >= 0.19, MhoUzd_on_Hav_A,  MhoUzd_on_Hav_B )

        Mho_exp = np.exp( -1.*(0.5*dragbld*vonkarm**(-2)*lambda_f)**(-0.5) )
        MhoUz0_on_Hav = ( (1. - MhoUzd_on_Hav)*Mho_exp)

        zd = bldhgt_ave*(MhoUzd_on_Hav + ((0.2375 * np.log(lambda_p) + 1.1738)*(sigma_h/bldhgt_ave)) )

        # z0 = bldhgt_ave*(MhoUz0_on_Hav + (np.exp(0.8867*lambda_f) - 1.)*(sigma_h/bldhgt_ave)**(np.exp(2.3271*lambda_f))) # ERROR IN KENT ET AL PAPER
        z0 = bldhgt_ave*(MhoUz0_on_Hav + np.exp(0.8867*lambda_f - 1.)*(sigma_h/bldhgt_ave)**(np.exp(2.3271*lambda_f))) # correction based on UMEP

    if mode == 3:
        print('zd,z0 Kanda et al 2013')

        alpha   = 4.43 # for staggered arrays
        beta    = 1.0  # for staggered arrays

        if np.all(bldhgt_max)==None:
            bldhgt_max = 1.5*sigma_h + bldhgt_ave # assume bldhgt_max
            # bldhgt_max = 12.51*sigma_h**0.77 # eq 3 from Kanda et al 2013

        # first calculate macdonald et al 1998 (for later scaling)
        mac_zd = (1. + (lambda_p - 1.)*alpha**(-lambda_p))*bldhgt_ave
        mac_z0 = ((1. - mac_zd/bldhgt_ave) * np.exp( -1.*(0.5*beta*(dragbld/vonkarm**2)*(1.-mac_zd/bldhgt_ave)*lambda_f)**(-0.5) ))*bldhgt_ave

        # then Kanda et al 2013 scaling parameters
        a0,b0,c0 = 1.29, 0.36, -0.17
        a1,b1,c1 = 0.71, 20.21, -0.77

        X = (sigma_h + bldhgt_ave)/bldhgt_max  # eq. 10b
        Y = (lambda_p*sigma_h)/bldhgt_ave        # eq. 12b

        zd = bldhgt_max*(c0*X**2 + (a0*lambda_p**b0 - c0)*X) # eq 10a 
        z0 = mac_z0*(a1 + b1*Y**2 + c1*Y) # eq 12a

    zd = np.round(zd,4)
    zd = xr.where(zd<0,0,zd)
    z0 = np.round(z0,4)
    z0 = xr.where(z0<0,0,z0)

    return zd,z0

def calc_derived_fractions(ds):

    '''create derived classes and attributes (normalising for cloud and shadow fractions)'''

    all_keys = [str(key) for key in ds.data_vars.keys()]

    # test for cloud and shadow to rescale variables

    if set(['cloud','shadow']).issubset(all_keys):
        denominator = 1. - ds['cloud'].fillna(0) - ds['shadow'].fillna(0)
    if set(['cloud']).issubset(all_keys):
        denominator = 1. - ds['cloud'].fillna(0)
    if set(['shadow']).issubset(all_keys):
        denominator = 1. - ds['shadow'].fillna(0)
    else:
        denominator = 1.

    total = ds.to_array().sum(axis=0)

    # fill with 0, then mask from total (counting any value)
    dsf = ds.fillna(0).where(total>0)

    derived = xr.Dataset()

    derived['building_fraction']  = (dsf['buildings'])/denominator
    derived['tree_fraction']      = (dsf['trees'])/denominator
    derived['grass_fraction']     = (dsf['grass'])/denominator
    derived['shrub_fraction']     = (dsf['other_veg'])/denominator
    derived['bareearth_fraction'] = (dsf['bare_earth'])/denominator
    derived['roadpath_fraction']  = (dsf['road_path'] + dsf['built_area'])/denominator
    derived['water_fraction']     = (dsf['water'] + dsf['swimming_pool'])/denominator
    derived['total_built']        = (dsf['road_path'] + dsf['built_area'] + dsf['buildings'])/denominator
    derived['total_vegetation']   = (dsf['grass'] + dsf['trees'] + dsf['other_veg'])/denominator
    derived['total_pervious']     = (dsf['grass'] + dsf['trees'] + dsf['other_veg'] + dsf['bare_earth'] + dsf['water'] + dsf['swimming_pool'])/denominator

    return derived

def remove_partial_sums(ds):

    keys = ['building_fraction','tree_fraction','grass_fraction','shrub_fraction','bareearth_fraction','roadpath_fraction','water_fraction']

    total = ds[keys].to_array().sum(axis=0)
    mask = total.where((total>0.99999) & (total<1.00001)).notnull()

    ds_masked = ds.where(mask)

    print('number of masked cells:')
    print(ds.count() - ds_masked.count())

    return ds_masked

def assert_tests(geoscape,derived,rtol=1E-6):

    '''test fractions sum to 1'''

    # test geoscape classess sum to 1
    total = np.full_like(geoscape.buildings.values,0)
    for key in geoscape.keys():
        total = total + geoscape[key]
    np.testing.assert_allclose(total.to_series().dropna().values, 1, rtol)

    # test derived classes sum to 1
    total = np.full_like(derived.building_fraction.values,0)
    for key in ['building_fraction','tree_fraction','grass_fraction','shrub_fraction','bareearth_fraction','roadpath_fraction','water_fraction']:
        total = total + derived[key]
    np.testing.assert_allclose(total.to_series().dropna().values, 1, rtol)

    # test total classes sum to 1
    total = np.full_like(derived.building_fraction.values,0)
    for key in ['total_built','total_pervious']:
        total = total + derived[key]
    np.testing.assert_allclose(total.to_series().dropna().values, 1, rtol)

    # test total classes sum to 1
    total = np.full_like(derived.building_fraction.values,0)
    for key in ['roadpath_fraction','total_pervious','building_fraction']:
        total = total + derived[key]
    np.testing.assert_allclose(total.to_series().dropna().values, 1, rtol)

    print('tests passed (all cells add to 1)')

    return

def rasterise_height_shapefile(buildings,tree_ds,key='bld_hgt'):

    from geocube.api.core import make_geocube
    # from geocube.rasterize import rasterize_image
    # from functools import partial

    print('making cube')
    ds = make_geocube(
        vector_data=buildings[[key,'geometry']],
        resolution=tree_ds.rio.resolution(),
        output_crs='epsg:4326',
        # rasterize_function=partial(rasterize_image, all_touched=False)
        )

    # latitude go from 90 to -90
    ds = ds.rename({'y':'latitude','x':'longitude' }).sortby('latitude',ascending=False)

    ds['spatial_ref'].attrs = {'grid_mapping_name':'latitude_longitude'}
    encoding = {var:{'zlib':True, '_FillValue':missing_float} for var in ds.data_vars}

    for var in ds.data_vars:
        ds[var].attrs = {}

    full_raster_hgt_fpath = f'{raster_hgt_fpath}_{key}.nc'

    print('saving netcdf')
    ds.to_netcdf(full_raster_hgt_fpath,format='NETCDF4',encoding=encoding)

    return

################################################################################

if __name__ == "__main__":

    ################################################################################
    # path setup

    datapath = f'{projpath}/data'
    outpath = f'{projpath}/outputs'   # netcdf/tif outputs
    plotpath = f'{projpath}/figures'  # figures output

    if grid in ['CCI','cci']: 
        template_fpath = f'{datapath}/grid_template_cci.nc' # template for regridding

    if domain == 'sample':
        tif_dir = f'{datapath}/Geoscape_Sample_data_2021/06_Surface_Cover_Trees'
        # place 30 m resolution tifs last for merge
        tif_fpaths = [
            f'{tif_dir}/NSW_SURFACECOVER_2M_Z56_16837.tif']
        if with_30m:
            tif_fpaths = tif_fpaths+[
            f'{tif_dir}/NSW_SURFACECOVER_30M_Z56.tif']
        tree_fpath = f'{tif_dir}/NSW_TREES_Z56_16837.tif'
        shp_fpath =  f'{datapath}/Geoscape_Sample_data_2021/01_Buildings/SHP/buildings.shp'

    if domain == 'sydney':
        # place 30 m resolution tifs last for merge
        tif_fpaths = [f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16824.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16825.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16827.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16828.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16829.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16830.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16831.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16832.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16833.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16834.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16835.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16836.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16837.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_2M_Z56_16838.tif']
        if with_30m:
            tif_fpaths = tif_fpaths+[
                         f'{datapath}/sydney/ACTNSW_SURFACECOVER_30M_Z55.tif',
                         f'{datapath}/sydney/NSW_SURFACECOVER_30M_Z56.tif']
        shp_fpath = f'{datapath}/sydney/geoscape_buildings_gsyd_2020.shp'
        tree_fpath = f'{datapath}/sydney/NSW_TREES_Z56_16837.tif'

    geoscape_outpath = f'{outpath}/{domain.capitalize()}_geoscape_derived_surface_cover_{grid}_{version}.nc'
    geopkg_outpath = f'{outpath}/{domain.capitalize()}_geoscape_buildings.gpkg'
    raster_hgt_fpath = f'{outpath}/{domain.capitalize()}_geoscape_building_height_rasterised'

    ################################################################################
    if not os.path.exists(tif_fpaths[0]):
        print(f'ERROR: geoscape tif not found at {tif_fpaths[0]}\n \
            Try downloading sample from https://geoscape.com.au/get-sample/ and extract into data folder')
        exit()

    ds = main()

    print(f'done! see output in {outpath}')
