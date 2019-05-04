from numpy import exp,power,log
from saturation_vapor_pressure import vapor_pressure
def air_humidity_method_qsat26air(air_temperature,surface_air_pressure,relative_humdity):
    """
    # computes saturation specific humidity
    
    Args:
        air_temperature (float): temperature air [degC]
        surface_air_pressure (float): pressure [mb]
        relative_humdity (float): relative humidity [0-100]
    Returns:
        air_humidity (float): surface air saturation specific humidity [g/kg]
    """
    es = vapor_pressure(air_temperature,surface_air_pressure)
    em = 0.01*relative_humdity*es
    air_humidity = 622.*em/(surface_air_pressure-0.378*em)
    return air_humidity

def air_humidity_method_bentamy(relative_humidity,air_temperature):
    """
    method from Bentamy
    
    Args:
        relative_humidity (float): in % [0-100]
        air_temperature (float): [degC]
    Returns:
        air_humidity (float): surface air saturation specific humidity [g/kg]
    
    """
    presint = 1013.25
    eps = 0.622
    a = -4.928
    b = 23.55
    c = -2937
    TA = air_temperature+ 273.0
    esatur = (power(TA,a))*(exp((b + (c/TA))*log(10)))
    Qsatur = (eps*esatur)/(presint - (1-eps)*esatur)
    Qsatur = 1000*Qsatur
    air_humidity = relative_humidity*Qsatur
    air_humidity = air_humidity/100.
    return air_humidity