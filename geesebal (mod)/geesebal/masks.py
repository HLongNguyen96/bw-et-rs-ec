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

#CLOUD REMOVAL

#FUNCTION TO MASK CLOUDS IN LANDSAT 5 AND 7 FOR SURFACE REFLECTANCE
def f_cloudMaskL57_SR(image):
    l57_qa_var = image.select('pixel_qa')
    l57_value1_var = l57_qa_var.eq(5440) #CLEAR, LOW SET
    l57_value2_var = l57_qa_var.eq(5504) #WATER, LOW SET
    l57_mask_var = l57_value1_var.Or(l57_value2_var)
    return image.updateMask(l57_mask_var)

#FUNCTION FO MASK CLOUD IN LANDSAT 8 AND 9 FOR SURFACE REFELCTANCE
def f_cloudMaskL89_SR(image):
    l89_qa_var = image.select('pixel_qa')
    l89_value1_var = l89_qa_var.eq(21824) #CLEAR, LOW SET
    l89_value2_var = l89_qa_var.eq(21888) #WATER, LOW SET
    l89_mask_var = l89_value1_var.Or(l89_value2_var)
    return image.updateMask(l89_mask_var)

#ALBEDO
#TASUMI ET AL(2008) FOR LANDSAT 5 AND 7
def f_albedoL5L7(image):
    alfa = image.expression(
      '(0.254*B1) + (0.149*B2) + (0.147*B3) + (0.311*B4) + (0.103*B5) + (0.036*B7)',{
        'B1' : image.select(['B']).multiply(0.0000275).subtract(0.2),
        'B2' : image.select(['GR']).multiply(0.0000275).subtract(0.2),
        'B3' : image.select(['R']).multiply(0.0000275).subtract(0.2),
        'B4' : image.select(['NIR']).multiply(0.0000275).subtract(0.2),
        'B5' : image.select(['SWIR_1']).multiply(0.0000275).subtract(0.2),
        'B7' : image.select(['SWIR_2']).multiply(0.0000275).subtract(0.2)
      }).rename('ALFA')
    #ADD BANDS
    return image.addBands(alfa)

#ALBEDO
#USING TASUMI ET AL. (2008) METHOD FOR LANDSAT 8 AND 9
#COEFFICIENTS FROM KE ET AL. (2016)
def f_albedoL8L9(image):
    alfa = image.expression(
      '(0.130*B1) + (0.115*B2) + (0.143*B3) + (0.180*B4) + (0.281*B5) + (0.108*B6) + (0.042*B7)',{  #// (Ke, Im  et al 2016)
        'B1' : image.select(['UB']).multiply(0.0000275).subtract(0.2),
        'B2' : image.select(['B']).multiply(0.0000275).subtract(0.2),
        'B3' : image.select(['GR']).multiply(0.0000275).subtract(0.2),
        'B4' : image.select(['R']).multiply(0.0000275).subtract(0.2),
        'B5' : image.select(['NIR']).multiply(0.0000275).subtract(0.2),
        'B6' : image.select(['SWIR_1']).multiply(0.0000275).subtract(0.2),
        'B7' : image.select(['SWIR_2']).multiply(0.0000275).subtract(0.2)
      }).rename('ALFA')
    #ADD BANDS
    return image.addBands(alfa)