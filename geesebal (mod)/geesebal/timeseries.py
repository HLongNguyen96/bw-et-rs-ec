#----------------------------------------------------------------------------------------#
#---------------------------------------//GEESEBAL//-------------------------------------#
#GEESEBAL - GOOGLE EARTH ENGINE APP FOR SURFACE ENERGY BALANCE ALGORITHM FOR LAND (SEBAL)
#CREATED BY: LEONARD LAIPELT, RAFAEL KAYSER, ANDERSON RUHOFF, AND AYAN FLEICHMANN
#PROJECT - ET BRASIL https://etbrasil.org/
#LAB - HIDROLOGIA DE GRANDE ESCALA [HGE] website: https://www.ufrgs.br/hge/author/hge/
#UNIVERSITY - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL - UFRGS
#RIO GRANDE DO SUL, BRAZIL

#DOI
#VERSION 0.1.1
#CONTACT US: leonardo.laipelt@ufrgs.br

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#MODIFIED BY: HOANG LONG NGUYEN
#LAB - ECOHYDROLOGY website: https://www.sallyethompson.com/
#UNIVERSITY OF WESTERN AUSTRALIA - UWA
#PERTH, AUSTRALIA

#CONTACT: hoanglong.nguyen@research.uwa.edu.au

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

#PYTHON PACKAGES
#Call EE
import ee
ee.Initialize()
from datetime import date
import datetime

#FOLDERS
from .landsatcollection import (fexp_landsat_5Coordinate, fexp_landsat_7Coordinate, 
                                fexp_landsat_8Coordinate, fexp_landsat_9Coordinate)
from .masks import (f_cloudMaskL57_SR, f_cloudMaskL89_SR, f_albedoL5L7, f_albedoL8L9)
from .meteorology import get_meteorology
from .tools import (fexp_spec_ind, LST_DEM_correction, fexp_radlong_up, fexp_radshort_down, 
                    fexp_radlong_down, fexp_radbalance, fexp_soil_heat, fexp_sensible_heat_flux)
from .endmembers import (fexp_cold_pixel, fexp_hot_pixel)
from .evapotranspiration import fexp_et

#TIMESRIES FUNCTION
class TimeSeries2():
    
    #ENDMEMBERS DEFAULT
    #ALLEN ET AL. (2013)
    def __init__(self, year_i, month_i, day_i, year_e, month_e, day_e,
                 cloud_cover, coordinate, whole, north, south, NDVI_cold=5,Ts_cold=20,NDVI_hot=10,Ts_hot=20):
        
        #INFORMATION
        self.coordinate = coordinate
        self.whole = whole
        self.north = north
        self.south = south
        self.cloud_cover = cloud_cover
        self.start_date = ee.Date.fromYMD(year_i,month_i,day_i)
        self.i_date = date(year_i,month_i,day_i)
        self.end_date = date(year_e,month_e,day_e)
        self.n_search_days = self.end_date - self.i_date
        self.n_search_days = self.n_search_days.days
        self.end_date = self.start_date.advance(self.n_search_days,'day')
        
        #COLLECTIONS
        self.collection_l9 = fexp_landsat_9Coordinate(self.start_date, self.end_date, self.coordinate, self.cloud_cover)
        self.collection_l8 = fexp_landsat_8Coordinate(self.start_date, self.end_date, self.coordinate, self.cloud_cover)
        self.collection_l7 = fexp_landsat_7Coordinate(self.start_date, self.end_date, self.coordinate, self.cloud_cover)
        self.collection_l5 = fexp_landsat_5Coordinate(self.start_date, self.end_date, self.coordinate, self.cloud_cover)
        
        #LIST OF IMAGES
        self.sceneListL9 = self.collection_l9.aggregate_array('system:index').getInfo()
        self.sceneListL8 = self.collection_l8.aggregate_array('system:index').getInfo()
        self.sceneListL7 = self.collection_l7.aggregate_array('system:index').getInfo()
        self.sceneListL5 = self.collection_l5.aggregate_array('system:index').getInfo()
        
        self.collection = self.collection_l9.merge(self.collection_l8).merge(self.collection_l7).merge(self.collection_l5)
        self.CollectionList = self.collection.sort("system:time_start").aggregate_array('system:index').getInfo()
        self.CollectionList_image = self.collection.aggregate_array('system:index').getInfo()
        self.count = self.collection.size().getInfo()
        
        #PRINT NUMBER OF SCENE
        print("Number of scenes: ",self.count)
        n = 0
        
        #LIST FOR ET AND DATES
        self.List_ET_coordinate = []
        self.List_ET_whole = []
        self.List_ET_north = []
        self.List_ET_south = []
        self.List_Date = []
        
        #====== ITERATIVE PROCESS ======#
        #FOR EACH IMAGE ON THE LIST
        #ESTIMATE ET DAILY IMAGE
        while n < self.count:
            
            #GET IMAGE
            self.image = self.collection.filterMetadata('system:index','equals',self.CollectionList[n]).first()
            self.image = ee.Image(self.image)
            
            #PRINT ID
            print(self.image.get('LANDSAT_PRODUCT_ID').getInfo())
            
            #GET INFORMATIONS FROM IMAGE
            self._index = self.image.get('system:index')
            self.cloud_cover = self.image.get('CLOUD_COVER')
            self.LANDSAT_ID = self.image.get('LANDSAT_PRODUCT_ID').getInfo()
            self.landsat_version = self.image.get('SPACECRAFT_ID').getInfo()
            self.sun_elevation = self.image.get('SUN_ELEVATION')
            self.azimuth_angle = self.image.get('SUN_AZIMUTH')
            self.zenith_angle = ee.Number(90).subtract(self.sun_elevation)
            self.time_start = self.image.get('system:time_start')
            self._date = ee.Date(self.time_start)
            self._year = ee.Number(self._date.get('year'))
            self._month = ee.Number(self._date.get('month'))
            self._day = ee.Number(self._date.get('day'))
            self._hour = ee.Number(self._date.get('hour'))
            self._minutes = ee.Number(self._date.get('minutes'))
            self.crs = self.image.projection().crs()
            self.transform = ee.List(ee.Dictionary(ee.Algorithms.Describe(self.image.projection())).get('transform'))
            self.date_string = self._date.format('YYYY-MM-dd').getInfo()
            
            #ENDMEMBERS
            self.p_top_NDVI = ee.Number(NDVI_cold)
            self.p_coldest_Ts = ee.Number(Ts_cold)
            self.p_lowest_NDVI = ee.Number(NDVI_hot)
            self.p_hottest_Ts = ee.Number(Ts_hot)
            
            #MASKS
            if self.landsat_version == 'LANDSAT_5':
               self.image_toa = ee.Image('LANDSAT/LT05/C02/T1_TOA/'+ self.CollectionList[n][4:])
               
               #GET CALIBRATED RADIANCE
               self.col_rad = ee.Algorithms.Landsat.calibratedRadiance(self.image_toa)
               self.col_rad = self.image.addBands(self.col_rad.select([5],["T_RAD"]))
               
               #CLOUD REMOTION
               self.image = ee.ImageCollection(self.image).map(f_cloudMaskL57_SR)
               
               #ALBEDO TASUMI ET AL. (2008)
               self.image = self.image.map(f_albedoL5L7)
            
            elif self.landsat_version == 'LANDSAT_7':
               self.image_toa = ee.Image('LANDSAT/LE07/C02/T1_TOA/'+ self.CollectionList[n][4:])
               
               #GET CALIBRATED RADIANCE
               self.col_rad = ee.Algorithms.Landsat.calibratedRadiance(self.image_toa)
               self.col_rad = self.image.addBands(self.col_rad.select([5],["T_RAD"]))
               
               #CLOUD REMOTION
               self.image = ee.ImageCollection(self.image).map(f_cloudMaskL57_SR)
               
               #ALBEDO TASUMI ET AL. (2008)
               self.image = self.image.map(f_albedoL5L7)
               
            elif self.landsat_version == 'LANDSAT_8':
               self.image_toa = ee.Image('LANDSAT/LC08/C02/T1_TOA/'+self.CollectionList[n][2:])
                
               #GET CALIBRATED RADIANCE
               self.col_rad = ee.Algorithms.Landsat.calibratedRadiance(self.image_toa)
               self.col_rad = self.image.addBands(self.col_rad.select([9],["T_RAD"]))
               
               #CLOUD REMOTION
               self.image = ee.ImageCollection(self.image).map(f_cloudMaskL89_SR)
               
               #ALBEDO TASUMI ET AL. (2008)
               self.image = self.image.map(f_albedoL8L9)
               
            else:
               self.image_toa = ee.Image('LANDSAT/LC09/C02/T1_TOA/'+self.CollectionList[n][2:])
               
               #GET CALIBRATED RADIANCE
               self.col_rad = ee.Algorithms.Landsat.calibratedRadiance(self.image_toa)
               self.col_rad = self.image.addBands(self.col_rad.select([9],["T_RAD"]))
               
               #CLOUD REMOTION
               self.image = ee.ImageCollection(self.image).map(f_cloudMaskL89_SR)
               
               #ALBEDO TASUMI ET AL. (2008)
               self.image = self.image.map(f_albedoL8L9)
               
            #GEOMETRY
            self.geometryReducer = self.image.geometry().bounds().getInfo()
            self.geometry_download = self.geometryReducer['coordinates']
            self.camada_clip = self.image.select('BRT').first()
            
            #METEOROLOGY PARAMETERS
            col_meteorology= get_meteorology(self.image,self.time_start)
            
            #AIR TEMPERATURE [C]
            self.T_air = col_meteorology.select('AirT_G')
                
            #WIND SPEED [M/S]
            self.ux = col_meteorology.select('ux_G')
            
            #RELATIVE HUMIDITY [%]
            self.UR = col_meteorology.select('RH_G')
            
            #POTENTIAL ET [W/M2]
            self.Pet = col_meteorology.select('Pet')
            
            #NET RADIATION 24H [W/M2]
            self.Rn24hobs = col_meteorology.select('Rn24h_G')
            
            #SRTM DATA ELEVATION
            SRTM_ELEVATION = 'USGS/SRTMGL1_003'
            self.srtm = ee.Image(SRTM_ELEVATION).clip(self.geometryReducer)
            self.z_alt = self.srtm.select('elevation')
            
            #GET IMAGE
            self.image = self.image.first()
        
            #SPECTRAL IMAGES (NDVI, EVI, SAVI, LAI, T_LST, e_0, e_NB, long, lat)
            self.image = fexp_spec_ind(self.image)

            #LAND SURFACE TEMPERATURE
            self.image = LST_DEM_correction(self.image, self.z_alt, self.T_air, self.UR, self.sun_elevation, self._hour, self._minutes)
            
            #COLD PIXEL 
            self.d_cold_pixel = fexp_cold_pixel(self.image, self.geometryReducer, self.p_top_NDVI, self.p_coldest_Ts)

            #COLD PIXEL NUMBER
            self.n_Ts_cold = ee.Number(self.d_cold_pixel.get('temp').getInfo())
            
            #HOT PIXEL
            self.d_hot_pixel = fexp_hot_pixel(self.image, self.geometryReducer, self.p_lowest_NDVI, self.p_hottest_Ts)
    
            #INSTANTANEOUS OUTGOING LONG-WAVE RADIATION [W/M2]
            self.image = fexp_radlong_up(self.image)
            
            #INSTANTANEOUS INCOMING SHORT-WAVE RADIATION [W/M2]
            self.image = fexp_radshort_down(self.image, self.z_alt, self.T_air, self.UR, self.sun_elevation)
            
            #INSTANTANEOUS INCOMING LONGWAVE RADIATION [W/M2]
            self.image = fexp_radlong_down(self.image, self.n_Ts_cold)
            
            #INSTANTANEOUS NET RADIATON BALANCE (Rn) [W/M2]
            self.image = fexp_radbalance(self.image)
            
            #SOIL HEAT FLUX (G) [W/M2]
            self.image = fexp_soil_heat(self.image)
            
            #SENSIBLE HEAT FLUX (H) [W/M2]
            self.image = fexp_sensible_heat_flux(self.image, self.ux, self.UR, self.Rn24hobs, self.Pet, 
                                                 self.d_cold_pixel, self.d_hot_pixel, self.date_string, self.geometryReducer)
            
            #DAILY EVAPOTRANSPIRATION (ET_24H) [MM/DAY]
            self.image = fexp_et(self.image,self.Rn24hobs)
            
            self.NAME_FINAL=self.LANDSAT_ID[:5]+self.LANDSAT_ID[10:17]+self.LANDSAT_ID[17:25]
            self.image=self.image.addBands([self.image.select('ET_24h').rename(self.NAME_FINAL)])
            
            #EXTRACT ET VALUE
            self.ET_coordinate = self.ET_daily.reduceRegion(reducer = ee.Reducer.first() ,geometry = self.coordinate, scale = 30, maxPixels = 1e14)
            self.ET_whole = self.ET_daily.reduceRegion(reducer = ee.Reducer.mean() ,geometry = self.whole, scale = 30, maxPixels = 1e14)
            self.ET_north = self.ET_daily.reduceRegion(reducer = ee.Reducer.mean() ,geometry = self.north, scale = 30, maxPixels = 1e14)
            self.ET_south = self.ET_daily.reduceRegion(reducer = ee.Reducer.mean() ,geometry = self.south, scale = 30, maxPixels = 1e14)
            
            #GET DATE AND DAILY ET
            self._Date = datetime.datetime.strptime(self.date_string,'%Y-%m-%d')
            self.ET_coordinate_get = ee.Number(self.ET_coordinate.get(self.NAME_FINAL)).getInfo()
            self.ET_whole_get = ee.Number(self.ET_whole.get(self.NAME_FINAL)).getInfo()
            self.ET_north_get = ee.Number(self.ET_north.get(self.NAME_FINAL)).getInfo()
            self.ET_south_get = ee.Number(self.ET_south.get(self.NAME_FINAL)).getInfo()
            
            #ADD LIST
            self.List_Date.append(self._Date)
            self.List_ET_coordinate.append(self.ET_coordinate_get)
            self.List_ET_whole.append(self.ET_whole_get)
            self.List_ET_north.append(self.ET_north_get)
            self.List_ET_south.append(self.ET_south_get)
                        
            #MOVE TO THE NEXT ITERATION
            n=n+1
        