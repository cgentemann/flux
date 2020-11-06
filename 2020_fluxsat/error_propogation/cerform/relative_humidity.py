from numpy import exp,power,log
def relative_humidity_method_fairall(air_temperature,air_pressure,air_specific_humidity):
    """
    method Fairall
    # computes relative humidity given T,P, & Q
    Args:
        air_temperature (float): in Celius
        air_pressure (float):  [mb]
        air_specific_humidity (float): [kg/kg]
    Returns:
        relative_humidity (float): [%] [0-100]
    """
    es = 6.1121*exp(17.502*air_temperature/(air_temperature+240.97))*(1.0007+3.46e-6*air_pressure)
    em = air_specific_humidity*air_pressure/(0.378*air_specific_humidity+0.622)
    relative_humidity = 100.*em/es
    return relative_humidity

def relative_humidity_method_bentamy(air_temperature,air_specific_humdity):
    """
    method Bentamy 
    Args:
        air_temperature (float): in Celius
        air_specific_humdity (float):  [kg/kg]
    Returns:
        relative_humidity (float) :[%] [0-0.1]
    """
    presint = 1013.25
    eps = 0.622
    a = -4.928
    b = 23.55
    c = -2937
    TA = air_temperature+ 273.0
    esatur = (power(TA,a))*(exp((b + (c/TA))*log(10)))
    Qsatur = (eps*esatur)/(presint - (1-eps)*esatur)
    relative_humidity = air_specific_humdity/Qsatur
    return relative_humidity