import ee
import matplotlib.pyplot as plt
import numpy as np
import geemap as gm
# import geemap.foliumap as emap
import pandas as pd
import geopandas as gpd
import datetime
import sys
from matplotlib.colors import LinearSegmentedColormap

# initialize Earth Engine
ee.Initialize()
# suppress warnings from geopandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


class plot_props(object):
    def __init__(self):
        '''
        returns plotting properties such as colormaps
        '''
        
        def get_et_cmap():
            et_palette = ['#DEC29B', '#E6CDA1', '#EDD9A6', '#F5E4A9', '#FFF4AD', '#C3E683', 
                          '#6BCC5C', '#3BB369', '#20998F', '#1C8691', '#16678A', '#114982', '#0B2C7A']
            et_cmap = LinearSegmentedColormap.from_list('etcmap', et_palette, N=100)
            return et_cmap
        def get_esi_cmap():
            esi_cmap = 'coolwarm_r'
            return esi_cmap
        self.ET_cmap = get_et_cmap()
        self.ESI_cmap= get_esi_cmap()

class data_review(object):

    def __init__(self, 
                 model, 
                 state_list, 
                 out_path, 
                 process_scale=1000,
                 start_date = '1990',
                 end_date = '2021',
                 cr_mask = True
                ):
        """
        Args:
            model: ptjpl, ssebop, eemetric, geesebal, disalexi
            
        """
        
        self.model = model
        self.band = self.get_band()
        self.state_list = state_list
        self.out_path = out_path
        self.process_scale = process_scale
        self.s_date = start_date
        self.e_date = end_date
        self.cr_mask = cr_mask
        ## Features
        self.huc8_all = ee.FeatureCollection('USGS/WBD/2017/HUC08')
        self.states = ee.FeatureCollection('TIGER/2018/States')
        self.state = self.get_state()
        self.state_gdf = gm.ee_to_gdf(self.state)
        self.state_hucs = self.get_huc8s()
        self.state_hucs_gdf = gm.ee_to_gdf(self.state_hucs)
        self.plot_bounds = self.get_plot_bounds()
        ## ET relevant
        self.ET = self.get_ET()
        self.ET_clim_m = self.get_month_clim()
        self.ET_clim_std = self.get_month_clim_std()
        self.ET_anomalies = self.get_ET_anoms()
        self.cloud_free_count = self.get_cloudfree_count()
        self.EnsRange = self.get_ensemble_range()

    ## Image & Image Collection Operations
    def get_state(self):
        ''' returns state geometry from list'''
        return self.states.filter(ee.Filter.inList('STUSPS',self.state_list));

    def get_huc8s(self):
        ''' returns huc 8 basins that intersect with a state
        if cr_mask = True, only huc 8s that are part of Colorado River Basin are returned'''
        h8 = ee.FeatureCollection('USGS/WBD/2017/HUC08').filter(ee.Filter.bounds(self.state.first().geometry()));
        if self.cr_mask == True:
            def add_hucnum(feature):
                num = ee.Number.parse(feature.get('huc8'));
                return feature.set('huc8num', num);
            sheds = h8.map(add_hucnum)
            huc8 = sheds.filter(ee.Filter.gte('huc8num', 14000000)).filter(ee.Filter.lte('huc8num', 16000000));
        else:
            huc8 = h8
        return huc8.select(ee.List(['huc8','name','areasqkm']))

    def get_plot_bounds(self):
        ''' returns the bounding box for the state'''
        coords  = self.state_hucs.geometry().bounds().coordinates().getInfo()[0]
        lon_min = np.array(coords)[:,0].min()-0.15
        lon_max = np.array(coords)[:,0].max()+0.15
        
        lat_min = np.array(coords)[:,1].min()-0.15
        lat_max = np.array(coords)[:,1].max()+0.15
        bbox = [lon_max, lat_min, lon_min, lat_max]
    
        return bbox
    
    def get_band(self):
        '''
        returns band associated with evapotranspiration model
        '''
        band = {'disalexi':'et',
                'eemetric':'et',
                'geesebal':'et',
                'ptjpl':'et',
                'ssebop':'et',
                'sims':'et',
                'ensemble':'et_ensemble_mad'
               }
        return ee.String(band[self.model])
    
    def get_ET(self):
        '''
        returns GEE asset_id associated with evapotranspiration model
            * long term convert this to a csv file to make easier to update & edit.
        '''
        model_loc={'disalexi':'projects/openet/assets/disalexi/conus/gridmet/monthly/provisional',
                   'eemetric':'projects/openet/assets/eemetric/conus/gridmet/monthly/provisional',
                   'geesebal':'projects/openet/assets/geesebal/conus/gridmet/monthly/provisional',
                   'ptjpl'   :'projects/openet/assets/ptjpl/conus/gridmet/monthly/provisional',
                   'ssebop'  :'projects/openet/assets/ssebop/conus/gridmet/monthly/provisional',
                   'sims'    :'projects/openet/assets/sims/conus/gridmet/monthly/provisional',
                   'ensemble':'projects/openet/assets/ensemble/conus/gridmet/monthly/provisional'
                  }
        model_loc2={'disalexi':'projects/openet/assets/disalexi/conus/gridmet/monthly/v2_0',
                   'eemetric':'projects/openet/assets/eemetric/conus/gridmet/monthly/v2_0',
                   'geesebal':'projects/openet/assets/geesebal/conus/gridmet/monthly/v2_0',
                   'ptjpl'   :'projects/openet/assets/ptjpl/conus/gridmet/monthly/v2_0',
                   'ssebop'  :'projects/openet/assets/ssebop/conus/gridmet/monthly/v2_0',
                   'sims'    :'projects/openet/assets/sims/conus/gridmet/monthly/v2_0',
                   'ensemble':'projects/openet/assets/ensemble/conus/gridmet/monthly/v2_0'
                  }
        ET  = ee.ImageCollection(model_loc[self.model])
        ET2 = ee.ImageCollection(model_loc2[self.model])
        combined = ET.merge(ET2)
        ETall = combined.filterDate(ee.Date(self.s_date),ee.Date(self.e_date)).select(self.band).filterBounds(self.state_hucs.geometry())
        return ETall

    def get_ensemble_range(self):
        '''
        returns ensemble range to 
        '''

        ET  = ee.ImageCollection('projects/openet/assets/ensemble/conus/gridmet/monthly/provisional')
        ET2 = ee.ImageCollection('projects/openet/assets/ensemble/conus/gridmet/monthly/v2_0')
        combined = ET.merge(ET2)
        ETall = combined.filterDate(ee.Date(self.s_date),ee.Date(self.e_date)).filterBounds(self.state_hucs.geometry()).select(['et_ensemble_mad','et_ensemble_mad_min','et_ensemble_mad_max'])
        Ensval = ETall.select('et_ensemble_mad').mean()
        Ensmin = ETall.select('et_ensemble_mad_min').mean()
        Ensmax = ETall.select('et_ensemble_mad_max').mean()
        per_range = ee.Image(Ensmax.subtract(Ensmin)).divide(Ensval).rename('percent_dif');
        # UT_mask = ee.Image("projects/openet/data_review/UT/common/mosaiced_nlcd_agricultural_mask_UT");
        UCRB_mask = ee.Image("projects/openet/data_review/UCRB/common/mosaiced_nlcd_agricultural_mask_UCRB");
        return per_range#.updateMask(UCRB_mask)
    
    def get_cloudfree_count(self):
        '''
        returns cloud free count as Image Collection by model
        '''
        model_loc={'disalexi':'projects/openet/assets/disalexi/conus/gridmet/monthly/provisional',
                   'eemetric':'projects/openet/assets/eemetric/conus/gridmet/monthly/provisional',
                   'geesebal':'projects/openet/assets/geesebal/conus/gridmet/monthly/provisional',
                   'ptjpl'   :'projects/openet/assets/ptjpl/conus/gridmet/monthly/provisional',
                   'ssebop'  :'projects/openet/assets/ssebop/conus/gridmet/monthly/provisional',
                   'sims'    :'projects/openet/assets/sims/conus/gridmet/monthly/provisional',
                   'ensemble':'projects/openet/assets/ensemble/conus/gridmet/monthly/provisional'
                  }
        model_loc2={'disalexi':'projects/openet/assets/disalexi/conus/gridmet/monthly/v2_0',
                   'eemetric':'projects/openet/assets/eemetric/conus/gridmet/monthly/v2_0',
                   'geesebal':'projects/openet/assets/geesebal/conus/gridmet/monthly/v2_0',
                   'ptjpl'   :'projects/openet/assets/ptjpl/conus/gridmet/monthly/v2_0',
                   'ssebop'  :'projects/openet/assets/ssebop/conus/gridmet/monthly/v2_0',
                   'sims'    :'projects/openet/assets/sims/conus/gridmet/monthly/v2_0',
                   'ensemble':'projects/openet/assets/ensemble/conus/gridmet/monthly/v2_0'
                  }
        CF1  = ee.ImageCollection(model_loc[self.model])
        CF2 = ee.ImageCollection(model_loc2[self.model])
        combined = CF1.merge(CF2)
        CFall = combined.filterDate(ee.Date(self.s_date),ee.Date(self.e_date)).select('count').filterBounds(self.state_hucs.geometry())
        return CFall
    
    def get_month_clim(self):
        '''
        returns the monthly climatology for ET
        '''
        def month_mean(m):
            '''
            returns the monthly climatology for 
            '''
            out_image = self.ET.filter(ee.Filter.calendarRange(m, m, 'month')).mean().set('month', m).clip(self.state.geometry())
            return out_image
        
        months = ee.List.sequence(1, 12);
        ET_clim = ee.ImageCollection.fromImages(months.map(month_mean))
        return ET_clim

    def get_month_clim_std(self):
        '''
        returns the monthly climatology for ET
        '''
        def month_stdev(m):
            '''
            returns the monthly climatology for 
            '''
            out_image = self.ET.filter(ee.Filter.calendarRange(m, m, 'month')).reduce(ee.Reducer.stdDev()).set('month', m).clip(self.state.geometry())
            return out_image
        
        months = ee.List.sequence(1, 12);
        ET_clim_std = ee.ImageCollection.fromImages(months.map(month_stdev))
        return ET_clim_std

    def get_ET_anoms(self):
        '''
        returns the ET anomalies less than self.thresh
        '''
        def get_month_num(img):
            '''
            returns the month of system:time_start from an image
            '''
            return ee.Date(img.get('system:time_start')).get('month');
        
        def get_mean_monthly_ET(mm):
            '''
            returns the climatology of ET_clim for given month mm
            '''
            return self.ET_clim_m.filter(ee.Filter.eq('month', mm)).first();

        def get_std_monthly_ET(mm):
            '''
            returns the climatology of ET_clim for given month mm
            '''
            return self.ET_clim_std.filter(ee.Filter.eq('month', mm)).first();
        
        def calc_esi(img):
            '''
            returns the pixels with ET anoms less than X% as 1 and others as 0
            '''
            img_mm = get_month_num(img)
            clim_m_et = get_mean_monthly_ET(img_mm)
            clim_std_et = get_std_monthly_ET(img_mm)
            esi = img.subtract(clim_m_et).divide(clim_std_et)
            out_band_name = ee.String('esi')
            out_esi = esi.rename(out_band_name).copyProperties(img, ["system:time_start"])#.unmask(0)

            return out_esi
        
        return self.ET.map(calc_esi)
    
def get_month_from_IC(IC, date_str):
    '''
    IC is image
    date_str is a date time
    '''
    date_of_interest = ee.Date(date_str);
    end_doi = date_of_interest.advance(1, 'month')
    return IC.filterDate(date_of_interest, end_doi).mean()
    
## Get Statistics by HUC
def get_zonal_stats(self):
    '''
    returns a dataframe of the zonal statistics for the county huc boundaries
    '''
    out_table_stats = self.out_path+str(self.county_name)+'_et_anomalies_by_huc.csv'
    gm.zonal_statistics(self.ET_anomalies, self.state_hucs, out_table_stats, statistics_type='MEAN', scale=self.process_scale)
    raw_df = pd.read_csv(out_table_stats)
    dfT = raw_df.T[:-3]
    new_index=[]
    for i in dfT.index:
        new_index.append(i.split('_')[1])
    dfT['date_time'] = pd.to_datetime(new_index)
    dfT.set_index('date_time',inplace=True)
    new_names = list(map(str,raw_df.huc12))
    dfT.columns=new_names
    return dfT

