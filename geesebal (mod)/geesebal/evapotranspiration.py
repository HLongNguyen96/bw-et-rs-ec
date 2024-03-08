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

def fexp_et(image, Rn24hobs):
    
    #NET DAILY RADIATION (Rn24h) [W/M2]
    #BRUIN (1982)
    Rn24hobs = Rn24hobs.select('Rn24h_G').multiply(ee.Number(1))
    
    #GET ENERGY FLUXES VARIABLES AND LST
    i_Rn = image.select('Rn')
    i_G = image.select('G')
    i_lst = image.select('T_LST_DEM')
    i_H_final = image.select('H')
    
    #FILTER VALUES
    i_H_final = i_H_final.where(i_H_final.lt(0),0)
    
    #INSTANSTANEOUS LATENT HEAT FLUX (LE) [W/M2]
    #BATIAANSSEN ET AL. (1998)
    i_lambda_ET = i_H_final.expression(
        '(i_Rn - i_G - i_H_final)',{
            'i_Rn': i_Rn,
            'i_G': i_G,
            'i_H_final': i_H_final}).rename('LE')
    
    #FILTER
    i_lambda_ET = i_lambda_ET.where(i_lambda_ET.lt(0),0)
    
    #LATENT HEAT OF VAPORISATION (LAMBDA) [J/KG]
    #BISHT ET AL. (2005)
    #LAGOUARDE AND BURNET (1983)
    i_lambda = i_H_final.expression(
        '(2.501-0.002361*(Ts-273.15))', {'Ts' : i_lst })
    
    #INSTANSTANEOUS ET (ET_inst) [MM/H]
    i_ET_inst = i_H_final.expression(
            '0.0036 * (i_lambda_ET/i_lambda)', {
            'i_lambda_ET' : i_lambda_ET,
            'i_lambda' : i_lambda  }).rename('ET_inst')

    #EVAPORATIVE FRACTION (EF)
    #CRAGO (1996)
    i_EF = i_H_final.expression(
            'i_lambda_ET/(i_Rn-i_G)',
            {'i_lambda_ET' : i_lambda_ET,
             'i_Rn' : i_Rn,
             'i_G' : i_G }).rename('EF')

    #DAILY EVAPOTRANSPIRATION (ET_24h) [MM DAY-1]
    i_ET24h_calc = i_H_final.expression(
            '(0.0864 *i_EF * Rn24hobs)/(i_lambda)', {
          'i_EF' : i_EF,
          'i_lambda' : i_lambda,
          'Rn24hobs' : Rn24hobs
          }).rename('ET_24h')

    #ADD BANDS
    image = image.addBands([i_ET_inst, i_ET24h_calc, i_lambda_ET, i_EF])
    return image    
    