from numpy import minimum,arctan,sqrt,log,exp,float64,real,ma
import pdb
def temperature_structure_method_psit_26(zet):
    """
    # computes temperature structure function
    Args:
        zet (float ndarray):
    Returns:
        psi (float ndarray):
    """
    zet0 = zet
    zet = ma.array(zet,dtype=complex)
    dzet = minimum(50,0.35*zet0) # stable
    psi = -((1+0.6667*zet)**1.5+0.6667*(zet0-14.28)*exp(-dzet)+8.525)
    psi = real(psi)
    k = zet0<0 # unstable
    x = (1-15*zet[k])**0.5
    x = real(x)
    psik = 2.*log((1.+x)/2.)
    x = (1-34.15*zet[k])**0.3333
    x = real(x)
    psic = 1.5*log((1+x+x**2)/3)-sqrt(3)*arctan((1+2*x)/sqrt(3))+4*arctan(1)/sqrt(3)
    f = zet0[k]**2/(1.+zet0[k]**2)
    psi[k] = (1-f)*psik+f*psic
    return psi


def temperature_structure_method_psit_30(zet):
    """
    # computes temperature structure function
    Args:
        zet (float ndarray):
    Returns:
        psi (float ndarray):
    """
    zet0 = zet.copy()
    zet = ma.array(zet,dtype=complex)
    x = (1.0-15.0*zet)**0.5
    x = real(x)
    psik = 2.*log((1.+x)/2.)
    x = (1-34.15*zet)**0.3333
    x= real(x)
    psic = 1.5*log((1+x+x*x)/3.)-sqrt(3.)*arctan((1.+2.*x)/sqrt(3.))+4.*arctan(1)/sqrt(3)
    f = zet0*zet0/(1+zet0*zet0)
    psi = (1-f)*psik+f*psic
    ii = zet0>0
    ccc = minimum(50.,0.35*zet0)
    psi[ii] = -((1+2/3*zet0[ii])**1.5+0.6667*(zet0[ii]-14.28)/exp(ccc[ii])+8.525)
    return psi