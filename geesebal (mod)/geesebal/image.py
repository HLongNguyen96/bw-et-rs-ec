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

#FOLDERS
from .masks import (f_cloudMaskL57_SR, f_cloudMaskL89_SR, f_albedoL5L7, f_albedoL8L9)
from .meteorology import get_meteorology
from .tools import (fexp_spec_ind, LST_DEM_correction, fexp_radlong_up, fexp_radshort_down, 
                    fexp_radlong_down, fexp_radbalance, fexp_soil_heat, fexp_sensible_heat_flux)
from .endmembers import (fexp_cold_pixel, fexp_hot_pixel)
from .evapotranspiration import fexp_et

#IMAGE FUNCTION
class Image2():
    
    #ENDMEMBERS DEFAULT
    #ALLEN ET AL. (2013)
    def __init__(self,image,NDVI_cold=5,Ts_cold=20,NDVI_hot=10,Ts_hot=20):
        
        #INFORMATION
        self.image = ee.Image(image)
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
            self.image = self.image.select([0,1,2,3,4,8,5,17],["B","GR","R","NIR","SWIR_1","BRT","SWIR_2","pixel_qa"])
            self.image_toa=ee.Image('LANDSAT/LT05/C02/T1_TOA/'+ self._index.getInfo())
            
            #GET CALIBRATED RADIANCE
            self.col_rad = ee.Algorithms.Landsat.calibratedRadiance(self.image_toa)
            self.col_rad = self.image.addBands(self.col_rad.select([5],["T_RAD"]))
            
            #CLOUD REMOTION
            self.image = ee.ImageCollection(self.image).map(f_cloudMaskL57_SR)
            
            #ALBEDO TASUMI ET AL. (2008)
            self.image = self.image.map(f_albedoL5L7)
            
        elif self.landsat_version == 'LANDSAT_7':
            self.image = self.image.select([0,1,2,3,4,8,5,17],["B","GR","R","NIR","SWIR_1","BRT","SWIR_2","pixel_qa"])
            self.image_toa = ee.Image('LANDSAT/LE07/C02/T1_TOA/'+ self._index.getInfo())
            
            #GET CALIBRATED RADIANCE
            self.col_rad = ee.Algorithms.Landsat.calibratedRadiance(self.image_toa)
            self.col_rad = self.image.addBands(self.col_rad.select([5],["T_RAD"]))
            
            #CLOUD REMOTION
            self.image = ee.ImageCollection(self.image).map(f_cloudMaskL57_SR)
            
            #ALBEDO TASUMI ET AL. (2008)
            self.image = self.image.map(f_albedoL5L7)
            
        elif self.landsat_version == 'LANDSAT_8':
            self.image = self.image.select([0,1,2,3,4,5,6,8,17],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
            self.image_toa = ee.Image('LANDSAT/LC08/C02/T1_TOA/'+ self._index.getInfo())
                
            #GET CALIBRATED RADIANCE
            self.col_rad = ee.Algorithms.Landsat.calibratedRadiance(self.image_toa)
            self.col_rad = self.image.addBands(self.col_rad.select([9],["T_RAD"]))
            
            #CLOUD REMOTION
            self.image = ee.ImageCollection(self.image).map(f_cloudMaskL89_SR)
            
            #ALBEDO TASUMI ET AL. (2008)
            self.image = self.image.map(f_albedoL8L9)
            
        else:
            self.image = self.image.select([0,1,2,3,4,5,6,8,17],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
            self.image_toa = ee.Image('LANDSAT/LC09/C02/T1_TOA/'+ self._index.getInfo())
                
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
         
