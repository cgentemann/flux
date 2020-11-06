from numpy import minimum,exp,arctan,sqrt,log,float64,ma,real
import logging
import pdb
import copy
def velocity_structure_method_psiu_26(zet):
    """
    # computes velocity structure function
    Args:
        zet (float ndarray):
    Returns:
        psi (float ndarray):
    """
    dzet = minimum(50,0.35*zet) # stable
    psi = -((1+zet)+0.6667*(zet-14.28)*exp(-dzet)+8.525)
    k = zet<0 # unstable
    zet0 = zet
    xou = (1-15*zet[k])**0.25
    xou = real(xou)
    psik = 2.*log((1.+xou)/2.)+log((1.+xou*xou)/2.)-2.*arctan(xou)+2.*arctan(1)
    xou = (1-10.15*zet[k])**0.3333
    xou = real(xou)
    zet = zet0
    psic = 1.5*log((1+xou+xou**2)/3.)-sqrt(3.)*arctan((1.+2.*xou)/sqrt(3.))+4.*arctan(1.)/sqrt(3.)
    f = zet[k]**2/(1.+zet[k]**2)
#     if hasattr(psi,'__size__'):
#         nb_elem = psi.size
#     else:
#         nb_elem = 1
#     if nb_elem>1:
    psi[k] = (1.-f)*psik+f*psic
#     else:
#         if zet<0:
#             psi = (1.-f)*psik+f*psic
    return psi

def velocity_structure_method_psiu_30(zet):
    """
    # computes velocity structure function
    #validated transcoding
    Args:
        zet (float ndarray):
    Returns:
        psi (float ndarray):
    """
    zet0 = zet.copy()
    zet = ma.array(zet,dtype=complex) #to avoid problem with power and decimal exposant (masked values for ma.array)
    xou = (1.-15.*zet)**0.25
    xou = real(xou)
    
    psik = 2.0*log((1.+xou)/2.)+log((1.+xou*xou)/2.)-2.*arctan(xou)+2.*arctan(1.)
    xou = (1.0-10.15*zet)**0.3333
    xou = real(xou)
    psic = 1.5*log((1.+xou+xou*xou)/3.)-sqrt(3.)*arctan((1.+2.*xou)/sqrt(3.))+4.*arctan(1.)/sqrt(3.)
#     zet = zet0
    fff = zet0*zet0/(1+zet0*zet0)
    psi = (1.-fff)*psik+fff*psic
    ii = zet0>0
    logging.debug('psiu_30  | ii:%s zet:%s psi:%s ',ii.shape,zet0.shape,psi.shape)
    ccc = minimum(50.,0.35*zet0)
    psi[ii] = -((1.+1.0*zet0[ii])**1.0+0.667*(zet0[ii]-14.28)/exp(ccc[ii])+8.525)
    return psi