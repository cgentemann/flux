import logging
from numpy import random,meshgrid,power,sin,pi
def grv(lat,lon=None,shape_wanted=None):
    """
    purpose:give gravity cste depending of the latitude
    note: grv module returns vector if lat is a vector, to force the output to be a matrix of the same shape
    than over input fields, one can provide the lon vector as well or directly the shape of the matrix wanted.
    Args:
        lat (float ndarray): latitudes
        lon (float ncarray): longitudes [optional]
        shape_wanted (tuple): shape of the matrix wanted [optional]
    Returns:
        gg (float ndarray): gravity constant
        
        2/12/2020 removed long which was a python 2.7 code
    """
    gamma = 9.7803267715
    c1 = 0.0052790414
    c2 = 0.0000232718
    c3 = 0.0000001262
    c4 = 0.0000000007
    
    logging.debug('grv | lat:%s',lat.shape)
    if lon is not None:
        lon_m,lat_m = meshgrid(lon,lat)
    elif shape_wanted is not None and isinstance(lat, (int, float, complex))==False and len(shape_wanted)>1:
        randy = random.rand(shape_wanted[1])
#         lat_m = tile(lat,shape_wanted[1])
        lon_m,lat_m = meshgrid(randy,lat)
#         lat_m = lat_m.T
        logging.debug('grv | shape_wanted :%s lat_m:%s',shape_wanted,lat_m.shape)
    else:
        lat_m = lat
    phi = lat_m*pi/180.
    xx = sin(phi)
    gg = gamma*(1+c1*power(xx,2)+c2*power(xx,4)+c3*power(xx,6)+c4*power(xx,8))
    return gg