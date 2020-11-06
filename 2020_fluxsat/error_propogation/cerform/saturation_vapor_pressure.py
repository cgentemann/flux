from numpy import exp
def vapor_pressure(temperature,pressure):
    """
    come from coare40vn.m
    computes saturation vapor pressure
    Args:
        pressure (float):  [Celsius]
        pressure (float):  [mb] 
    Returns:
        saturation_vapor_pressure (float):  [mb]
    """
    saturation_vapor_pressure = 6.1121*exp(17.502*temperature/(temperature+240.97))*(1.0007+3.46e-6*pressure)
    return saturation_vapor_pressure