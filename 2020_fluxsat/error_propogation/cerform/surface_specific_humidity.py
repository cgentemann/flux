from saturation_vapor_pressure import vapor_pressure
from numpy import power,exp,log
def sea_humidity_method_qsat26sea(sea_surface_temperature,surface_air_pressure):
    """
    come from coare4
    # computes surface saturation specific humidity [g/kg]
    Args:
        sea_surface_temperature (float): temperature sst [degC]
        surface_air_pressure (float): pressure [mb]
    Returns:
        sea_surface_specific_humidity (float): sea surface saturation specific humidity [g/kg]
    """
    ex = vapor_pressure(sea_surface_temperature,surface_air_pressure)
    es = 0.98*ex # reduction at sea surface
    sea_surface_specific_humidity = 622*es/(surface_air_pressure-0.378*es)
    return sea_surface_specific_humidity

def sea_humidity_method_bentamy(sea_surface_temperature):
    """
    method from Bentamy:
    Args:
        sea_surface_temperature (float): temperature sst [degC]
    Returns:
        sea_surface_specific_humidity (float): sea surface saturation specific humidity [g/kg]
    """
    TSM = sea_surface_temperature + 273.15
    a = -4.928
    b = 23.55
    c = -2937.0
    Psurf = 1013.25
    esurf = (power(TSM,a))*(exp((b + (c/TSM))*log(10)))
    Qsurf = (0.622*esurf)/(Psurf - esurf)
    sea_surface_specific_humidity = 1000.0*Qsurf
#     qsint  = Qsurf; 
    return sea_surface_specific_humidity

def sea_humidity_method_qsee(sea_surface_temperature,surface_air_pressure):
    """
    come from coare3
    # computes sea surface specific humidity [mb]
    Args:
        sea_surface_temperature (float): temperature [Celsius]
        surface_air_pressure (float): pressure [mb] 
    Returns:
        sea_surface_specific_humidity (float): sea surface saturation specific humidity [g/kg]
    """
#     x = y[:,0]
#     p = y[:,1]
#     x,p = y
    es = 6.112*exp(17.502*sea_surface_temperature/(sea_surface_temperature+240.97))*0.98*(1.0007+3.46e-6*surface_air_pressure)
    sea_surface_specific_humidity = es*621.97/(surface_air_pressure-0.378*es)
    return sea_surface_specific_humidity