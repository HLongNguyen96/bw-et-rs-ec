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
# ee.Initialize()

#SURFACE REFLECTANCE
#ATMOSPHERICALLY CORRECTED

#GET LANDSAT 9 COLLECTIONS BY PATH ROW
def fexp_landsat_9PathRow(start_date,end_date,n_path,n_row,th_cloud_cover):
    col_SR_L9 = (ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
                        .filterDate(start_date, end_date)
                        .filterMetadata('WRS_PATH', 'equals', n_path)
                        .filterMetadata('WRS_ROW', 'equals', n_row)
                        .select([0,1,2,3,4,5,6,8,17],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));
    return col_SR_L9

#GET LANDSAT 8 COLLECTIONS BY PATH ROW
def fexp_landsat_8PathRow(start_date,end_date,n_path,n_row,th_cloud_cover):
    col_SR_L8 = (ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                        .filterDate(start_date, end_date)
                        .filterMetadata('WRS_PATH', 'equals', n_path)
                        .filterMetadata('WRS_ROW', 'equals', n_row)
                        .select([0,1,2,3,4,5,6,8,17],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));
    return col_SR_L8

#GET LANDSAT 7 COLLECTIONS BY PATH ROW
def fexp_landsat_7PathRow(start_date,end_date,n_path,n_row,th_cloud_cover):
    col_SR_L7 = (ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
                        .filterDate(start_date, end_date)
                        .filterMetadata('WRS_PATH', 'equals', n_path)
                        .filterMetadata('WRS_ROW', 'equals', n_row)
                        .select([0,1,2,3,4,8,5,17],["B","GR","R","NIR","SWIR_1","BRT","SWIR_2","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover))
    return col_SR_L7

#GET LANDSAT 5 COLLECTIONS BY PATH ROW
def fexp_landsat_5PathRow(start_date,end_date,n_path,n_row,th_cloud_cover):
    col_SR_L5 = (ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
                        .filterDate(start_date, end_date)
                        .filterMetadata('WRS_PATH', 'equals', n_path)
                        .filterMetadata('WRS_ROW', 'equals', n_row)
                        .select([0,1,2,3,4,8,5,17],["B","GR","R","NIR","SWIR_1","BRT","SWIR_2","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover))
    return col_SR_L5

#GET LANDSAT 9 COLLECTIONS BY COORDINATE
def fexp_landsat_9Coordinate(start_date,end_date,coordinate,th_cloud_cover):
    col_SR_L9 = (ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
                        .filterDate(start_date, end_date)
                        .filterBounds(coordinate)
                        .select([0,1,2,3,4,5,6,8,17],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));
    return col_SR_L9

#GET LANDSAT 8 COLLECTIONS BY COORDINATE
def fexp_landsat_8Coordinate(start_date,end_date,coordinate,th_cloud_cover):
    col_SR_L8 = (ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                        .filterDate(start_date, end_date)
                        .filterBounds(coordinate)
                        .select([0,1,2,3,4,5,6,8,17],["UB","B","GR","R","NIR","SWIR_1","SWIR_2","BRT","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover));
    return col_SR_L8

#GET LANDSAT 7 COLLECTIONS BY COORDINATE
def fexp_landsat_7Coordinate(start_date,end_date,coordinate,th_cloud_cover):
    col_SR_L7 = (ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
                        .filterDate(start_date, end_date)
                        .filterBounds(coordinate)
                        .select([0,1,2,3,4,8,5,17],["B","GR","R","NIR","SWIR_1","BRT","SWIR_2","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover))
    return col_SR_L7

#GET LANDSAT 5 COLLECTIONS BY COORDINATE
def fexp_landsat_5Coordinate(start_date,end_date,coordinate,th_cloud_cover):
    col_SR_L5 = (ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
                        .filterDate(start_date, end_date)
                        .filterBounds(coordinate)
                        .select([0,1,2,3,4,8,5,17],["B","GR","R","NIR","SWIR_1","BRT","SWIR_2","pixel_qa"])
                        .filterMetadata('CLOUD_COVER', 'less_than', th_cloud_cover))
    return col_SR_L5