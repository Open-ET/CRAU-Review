import ee
import matplotlib.pyplot as plt
import numpy as np
import geemap as gm
import geemap.foliumap as emap
import pandas as pd
import geopandas as gpd
import datetime
import sys
# initialize Earth Engine
ee.Initialize()
# suppress warnings from geopandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


class fire_meta(object):
    
    def __init__(self,
                 state_list=ee.List(['CA','OR'])):
        """
            Args:
                model: ptjpl, ssebop, eemetric, geesebal, disalexi

        """                  
                 
        self.state_list = state_list
        self.counties = ee.FeatureCollection('TIGER/2018/Counties')
        self.states = ee.FeatureCollection('TIGER/2018/States')

    def get_county_list(self):
        out_list = self.counties.filter(ee.Filter.bounds(self.states.filter(ee.Filter.inList('STUSPS',self.state_list))))
        return out_list.aggregate_array('NAME').getInfo()
        

class fire_risk(object):

    def __init__(self, 
                 model, 
                 county_name, 
                 out_path, 
                 process_scale=1000,
                 start_date = '2016-07',
                 end_date = '2021'
                ):
        """
        Args:
            model: ptjpl, ssebop, eemetric, geesebal, disalexi
            
        """
        
        self.model = model
        self.band = self.get_band()
        self.thresh = -0.2
        self.county_name = county_name
        self.out_path = out_path
        self.process_scale = process_scale
        self.s_date = start_date
        self.e_date = end_date
        ## Features
        self.huc12_all = ee.FeatureCollection('USGS/WBD/2017/HUC12')
        self.counties = ee.FeatureCollection('TIGER/2018/Counties')
        self.states = ee.FeatureCollection('TIGER/2018/States')
        self.westernUS = self.get_western_us_States()
        self.county = self.get_county_boundary()
        self.county_hucs = self.huc12_all.filterBounds(self.county).select(ee.List(['huc12','areaacres']))
        self.county_hucs_gdf = gm.ee_to_gdf(self.county_hucs)
        self.all_fires = ee.FeatureCollection('USFS/GTAC/MTBS/burned_area_boundaries/v1').select(ee.List(['BurnBndAc','Ig_Date']))
        self.county_fires = self.all_fires.filterBounds(self.county_hucs)
        ## ET relevant
        self.ET = self.get_ET()
        self.ET_clim_m = self.get_month_clim()
        self.ET_anomalies = self.get_ET_anoms()

    ## Feature & Feature Collection Operations
    
    def get_western_us_States(self):
        """
        returns western US states for analysis
        """
        return self.states.filter(ee.Filter.inList('STUSPS',['CA','OR']))
    
    def get_county_boundary(self):
        """
        returns county boundary for county_name
        """
        countyWUS = self.counties.filterBounds(self.westernUS)
        county = countyWUS.filter(ee.Filter.eq('NAME',self.county_name))
        return county
        
    def assign_fire_times(self):
        '''
        assigns datetime to 
        '''
        fires_gdf = gm.ee_to_gdf(self.county_fires)
        fires_gdf['date_time'] = pd.to_datetime(fires_gdf['Ig_Date'], unit='ms')
        fires_gdf.set_index('date_time',inplace=True)
        return fires_gdf
    
    def filter_fires_by_year(self, yyyy):
        ''' 
        returns fire gdf filtered by year 
        inputs: yyyy [int]
        '''
        return self.assign_fire_times().loc[str(yyyy)+'-05-01':str(yyyy)+'-10-31']
    
    def co_huc_fires_by_year(self, yyyy):
        '''
        returns a geodataframe of huc12s with fire occurences
        '''
        huc_gdf = gm.ee_to_gdf(self.county_hucs)
        fire_gdf = self.filter_fires_by_year(str(yyyy))
        return huc_gdf.sjoin(fire_gdf)
    
    ## Image & Image Collection Operations
    
    def get_band(self):
        '''
        returns band associated with evapotranspiration model
        '''
        band = {'disalexi':'et',
                'eemetric':'et',
                'geesebal':'et',
                'ptjpl':'et',
                'ssebop':'et'
               }
        return ee.String(band[self.model])
    
    def get_ET(self):
        '''
        returns GEE asset_id associated with evapotranspiration model
            * long term convert this to a csv file to make easier to update & edit.
        '''
        model_loc={'disalexi':'projects/openet/disalexi/conus_gridmet/monthly_provisional',
                   'eemetric':'projects/openet/eemetric/conus_gridmet/monthly_provisional',
                   'geesebal':'projects/openet/geesebal/conus_gridmet/monthly_provisional',
                   'ptjpl':'projects/openet/ptjpl/conus_gridmet/monthly_provisional',
                   'ssebop':'projects/openet/ssebop/conus_gridmet/monthly_provisional'
                  }
        
        ET = ee.ImageCollection(model_loc[self.model]).filterDate(ee.Date(self.s_date), ee.Date(self.e_date)).select(self.band)
        return ET
            
    def get_month_clim(self):
        '''
        returns the monthly climatology for ET
        '''
        def month_mean(m):
            '''
            returns the monthly climatology for 
            '''
            out_image = self.ET.filter(ee.Filter.calendarRange(m, m, 'month')).mean().set('month', m).clip(self.county_hucs.geometry())
            return out_image
        
        months = ee.List.sequence(1, 12);
        ET_clim = ee.ImageCollection.fromImages(months.map(month_mean))
        return ET_clim
    
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
        
        def find_stress(img):
            '''
            returns the pixels with ET anoms less than X% as 1 and others as 0
            '''
            img_mm = get_month_num(img)
            clim_m_et = get_mean_monthly_ET(img_mm)
            per_dif =  img.subtract(clim_m_et).divide(clim_m_et)
            et_thresh = per_dif.select(self.band).lte(self.thresh);
            out_band_name = ee.String('stress')
            et_conditional = et_thresh.rename(out_band_name).unmask(0)

            return et_conditional
        
        return self.ET.map(find_stress)
    
    ## Get Statistics by HUC
    def get_zonal_stats(self):
        '''
        returns a dataframe of the zonal statistics for the county huc boundaries
        '''
        out_table_stats = self.out_path+str(self.county_name)+'_et_anomalies_by_huc.csv'
        gm.zonal_statistics(self.ET_anomalies, self.county_hucs, out_table_stats, statistics_type='MEAN', scale=self.process_scale)
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
    
    ## Mapping & Plotting Operations
