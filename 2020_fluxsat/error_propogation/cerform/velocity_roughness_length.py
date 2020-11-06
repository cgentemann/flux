from numpy import power,exp,sqrt 
def get_velocity_roughness_length(Cdn_10):
    """
    author: Antoine Grouazel
    Args:
        Cdn_10 (float):10-m velocity drag coeeficient, including gustiness
    Returns:
        z0 (float) : velocity_roughness_length [m]
    """
    von = 0.4 #cste copy paste from coare4
    z0 = 10./exp(sqrt((1000.*power(von,2))/Cdn_10))
    return z0