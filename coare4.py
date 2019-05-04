"""
.. module::cerform.flux.coare4

@author: Antoine Grouazel
@creation: 20/08/2015
@purpose: use coare4 in python for OHF project (compare buoys flux with L3 flux products)
@validation: tbd
#test

"""
from numpy import log,exp,sqrt,where,power,ones,minimum,maximum,ma,ndarray
import pdb
from gravity_constant import grv
from relative_humidity import relative_humidity_method_fairall
from surface_specific_humidity import sea_humidity_method_qsat26sea
from air_specific_humidity import air_humidity_method_qsat26air
from velocity_structure import velocity_structure_method_psiu_26
from temperature_structure import temperature_structure_method_psit_26
import logging
from default_values import complete_missing_args

def coare4(inputs):
    """
    Purpose: 
        transcoding from coare40vn.m coming from vectorized version of COARE3 code (Fairall et al, 2003) with modification based on the CLIMODE, MBL and CBLAST experiments (Edson et al., 2011). The cool skin option is retained but warm layer and surface wave options removed. 
        An important component of this code is whether the inputed ts 
        represents the skin temperature of a near surface temperature.  
        How this variable is treated is determined by the jcool parameter:
        set jcool=1 if Ts is bulk ocean temperature (default),
            jcool=0 if Ts is true ocean skin temperature. 
        The code assumes u,t,rh,ts are vectors; 
        sensor heights zu,zt,zl, latitude lat, and PBL height zi are constants;
        air pressure P and radiation Rs,Rl may be vectors or constants. 
        Default values are assigned for P,Rs,Rl,lat,and zi if these data are not 
        available.  Input NaNs to indicate no data. Defaults should be set to 
        representative regional values if possible.
    
    Args:  
        u (float): relative wind speed (m/s) at height zu(m)
        t (float): bulk air temperature (degC) at height zt(m)
        rh (float): relative humidity (#) at height zq(m) [0-100]
        P (float): surface air pressure (mb) (default = 1015)
        ts (float): water temperature (degC) see jcool below
        Rs (float): downward shortwave radiation (W/m^2) (default = 150) 
        Rl (float): downward longwave radiation (W/m^2) (default = 370)
        lat (float): latitude (default = +45 N)
        zi (float): PBL height (m) (default = 600m)
        
    Returns:
        usr (float): friction velocity that includes gustiness (m/s)
        tau (float): wind stress (N/m^2)
        hsb (float) : sensible heat flux into ocean (W/m^2)
        hlb (float) : latent heat flux into ocean (W/m^2)
        hbb (float) : buoyany flux into ocean (W/m^2)
        hsbb (float) : "sonic" buoyancy flux measured directly by sonic anemometer 
        tsr (float) : temperature scaling parameter (K)
        qsr (float) : specific humidity scaling parameter (g/Kg)
        zot (float) : thermal roughness length (m)
        zoq (float) : moisture roughness length (m)
        Cd (float) : wind stress transfer (drag) coefficient at height zu   
        Ch (float) : sensible heat transfer coefficient (Stanton number) at height zu   
        Ce (float) : latent heat transfer coefficient (Dalton number) at height zu
        L (float) : Obukhov length scale (m) 
        zet (float) : Monin-Obukhov stability parameter zu/L 
        dter (float) : cool-skin temperature depression (degC)
        dqer (float) : cool-skin humidity depression (degC)
        tkt (float) : cool-skin thickness (m)
        Urf (float) : wind speed at reference height (user can select height below)
        Tfr (float) : temperature at reference height
        Qfr (float) : specific humidity at reference height
        RHfr (float) : relative humidity at reference height
        UrfN (float) : neutral value of wind speed at reference height
        Rnl (float) : Upwelling IR radiation computed by COARE
        Le (float) : latent heat of vaporization
        rhoa (float) : density of air
        UN (float) : neutral value of wind speed at zu
        U10 (float) : wind speed adjusted to 10 m
        UN10 (float) : neutral value of wind speed at 10m
        Cdn_10 (float) : neutral value of drag coefficient at 10m    
        Chn_10 (float) : neutral value of Stanton number at 10m    
        Cen_10 (float) : neutral value of Dalton number at 10m    
    
    
    Notes:
        1) u is the relative wind speed, i.e., the magnitude of the difference between the wind (at zu) and ocean surface current vectors.
        2) Set jcool=0 in code if ts is true surface skin temperature,otherwise ts is assumed the bulk temperature and jcool=1.
        3) Set P=NaN to assign default value if no air pressure data available. 
        4) Set Rs=NaN, Rl=NaN if no radiation data available.  This assigns default values to Rs, Rl so that cool skin option can be applied. 
        5) Set lat=NaN and/or zi=NaN to assign default values if latitude and/or PBL height not given. 
        6) The code to compute the heat flux caused by precipitation is included if rain data is available (default is no rain).
        7) Code updates the cool-skin temperature depression dter and thickness tkt during iteration loop for consistency.
        8) Number of iterations set to nits = 6.
    
    Reference:
    
      Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
      Bulk parameterization of air sea fluxes: updates and verification for the 
      COARE algorithm, J. Climate, 16, 571-590.
    
    History:
     
     1. 12/14/05 - created based on scalar version coare26sn.m with input
        on vectorization from C. Moffat.  
     2. 12/21/05 - sign error in velocity_structure_method_psiu_26 corrected, and code added to use variable
        values from the first pass through the iteration loop for the stable case
        with very thin M-O length relative to zu (zetu>50) (as is done in the 
        scalar coare26sn and COARE3 codes).
     3. 7/26/11 - S = dt was corrected to read S = ut.
     4. 7/28/11 - modification to roughness length parameterizations based 
        on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
    
    
    Warning: 
        COARE4 bulk formula has not been officialy released by Fairall et al thus we do not provide this tool on behalf of them and we do not guarantee any results even if the code has been validated on many datasets.
    """
    list_input_needed = ['u','zu','t','zt','rh','zq','P','ts','Rs','Rl','lat','zi']
    inputs_complet = complete_missing_args(list_input_needed,inputs)
    
    u = inputs_complet['u']
    zu = inputs_complet['zu']
    t = inputs_complet['t']
    zt = inputs_complet['zt']
    rh = inputs_complet['rh']
    zq = inputs_complet['zq']
    P = inputs_complet['P']
    ts = inputs_complet['ts']
    Rs = inputs_complet['Rs']
    Rl = inputs_complet['Rl']
    lat = inputs_complet['lat']
    zi = inputs_complet['zi']
    
    
    jcool = 1
    
    
    # input variable u is assumed relative wind speed (magnitude of difference
    # between wind and surface current vectors). to follow orginal Fairall code, set
    # surface current speed us=0. if us data are available, construct u prior to
    # using this code.
    us = 0.*u
    logging.info('coare4 | RH: %s',ma.median(rh))
    # convert rh to specific humidity
    Qs = sea_humidity_method_qsat26sea(ts,P)/1000.0    # surface water specific humidity
    Q = air_humidity_method_qsat26air(t,P,rh)/1000.0  # specific humidity of air 
    #the rest of the code need to use mg/kg spec humidity
    logging.info('coare4 | Q:%s Qs:%s',ma.median(Q),ma.median(Qs))
    # set rain to zero
    rain = 0.*u # rain rate (mm/hr) - keep as option
    
    #***********  set constants **********************************************
    Beta = 1.2
    von  = 0.4
    fdg  = 1.00 # Turbulent Prandtl number
    tdk  = 273.16
    grav = grv(lat,shape_wanted=t.shape)
    
    #***********  air constants **********************************************
    Rgas = 287.1
    Le   = (2.501-.00237*ts)*1e6
    cpa  = 1004.67
    cpv  = cpa*(1+0.84*Q)
    rhoa = P*100./(Rgas*(t+tdk)*(1+0.61*Q))
    visa = 1.326e-5*(1.+6.542e-3*t+8.301e-6*power(t,2)-4.84e-9*power(t,3))
    
    #***********  cool skin constants  ***************************************
    Al   = 2.1e-5*power((ts+3.2),0.79)
    be   = 0.026
    cpw  = 4000.
    rhow = 1022.
    visw = 1e-6
    tcw  = 0.6
    logging.debug('coare4| shapes grav:%s',grav.shape)
    bigc = 16.*grav*cpw*power((rhow*visw),3)/(power(tcw,2)*power(rhoa,2))
#     if isinstance(bigc,ndarray)==False:
#         bigc = ma.array([bigc])
    wetc = 0.622*Le*Qs/(Rgas*power((ts+tdk),2))
    
    #***********  net radiation fluxes ***************************************
    Rns = 0.945*Rs # albedo correction
    # IRup = eps*sigma*T^4 + (1-eps)*IR
    # Rnl = IRup - IR
    # Rll = eps*sigma*T^4 - eps*IR  as below
    
    Rnl = 0.97*(5.67e-8*power((ts-0.3*jcool+tdk),4)-Rl) # initial value
    
    # IRup = Rnl + IR
    
    #****************  begin bulk loop ********************************************
    
    #***********  first guess ************************************************
    du = u-us
    logging.debug('coare4 | ts:%s t:%s zt:%s Qs:%s Q:%s',ts,t,zt,Qs,Q)
    dt = ts-t-.0098*zt
    dq = Qs-Q
    ta = t+tdk
    ug = 0.5
    dter  = 0.3
    ut    = sqrt(power(du,2)+power(ug,2))
    u10   = ut*log(10/1e-4)/log(zu/1e-4)
    usr   = 0.035*u10
    zo10  = 0.011*power(usr,2)/grav + 0.11*visa/usr
    Cd10  = power((von/log(10./zo10)),2)
    Ch10  = 0.00115
    Ct10  = Ch10/sqrt(Cd10)
    zot10 = 10./exp(von/Ct10)
    Cd    = power((von/log(zu/zo10)),2)
    Ct    = von/log(zt/zot10)
    CC    = von*Ct/Cd
    Ribcu = -zu/zi/.004/power(Beta,3)
    logging.debug('coare4 | grav:%s zu:%s ta:%s dt:%s dter:%s jcool:%s dq:%s ut:%s',grav,zu,ta,dt,dter,jcool,dq,ut)
    Ribu  = -grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/power(ut,2)
    zetu = CC*Ribu*(1.+27./9.*Ribu/CC)
    k50 = zetu>50 # stable with very thin M-O length relative to zu
    logging.debug('coare4 | k50 indice: %s',k50)
    logging.debug('coare4 | Ribu %s',Ribu)
    k = Ribu<0
    logging.debug('coare4 | k indice: %s %s',k,type(k))
#     if hasattr(zetu,'__size__'):
#         nb_elem = zetu.size
#     else:
#         nb_elem = 1
#     if nb_elem>1:
    zetu[k] = CC[k]*Ribu[k]/(1+Ribu[k]/Ribcu)
#     else:
#         if Ribu<0:
#             zetu = CC*Ribu/(1+Ribu/Ribcu)
    L10 = zu/zetu
    gf = ut/du
    usr = ut*von/(log(zu/zo10)-velocity_structure_method_psiu_26(zu/L10))
    tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot10)-temperature_structure_method_psit_26(zt/L10))
    qsr = -(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-temperature_structure_method_psit_26(zq/L10))
#     tkt = 0.001*ones((N,1))
    tkt = 0.001*ones(t.shape) #agrouaze matricial 
    
    #**********************************************************
    #  The following gives the new formulation for the
    #  Charnock variable
    #**********************************************************
    
#     charn = 0.011*ones((N,1))
    charn = 0.011*ones(t.shape) #agrouaze matricial 
    umax = 22
    a1 = 0.0016
    a2 = -0.0035
    charn = a1*u10+a2
    k = u10>umax
#     if hasattr(charn,'__size__'):
#         nb_elem = charn.size
#     else:
#         nb_elem = 1
#     if nb_elem>1:
    charn[k] = a1*umax+a2
#     else:
#         if u10>umax:
#             charn = a1*umax+a2
    nits = 6 # number of iterations
    
    #The following parameter determines which version of the moisture 
    #roughness length.
    #  0: Step this to 0 to use the form of the moisture rough derived by 
    #  Fairall et al. (2003).
    #
    #  1: Step this to 1 to use the form of the moisture rough determined by the
    #  CLIMODE, GASEX and CBLAST data.
    #
    #  Note that the thermal roughness length gives a Stanton number that is
    #  very similar to COARE 3.0.
    #
    climodeversion = 1
    
    #**************  bulk loop **************************************************
    
    for i in range(nits):
        zet = von*grav*zu/ta*(tsr +.61*ta*qsr)/power(usr,2)
        zo = charn*power(usr,2)/grav+0.11*visa/usr # surface roughness
        rr = zo*usr/visa
        logging.debug('coare4 | zu:%s zet:%s',zu,zet)
        L = zu/zet
        zot = minimum(1.0e-4/power(rr,0.55),2.4e-4/power(rr,1.2)) # temp roughness
        if climodeversion:
            zoq = minimum(2.0e-5/power(rr,0.22),1.1e-4/power(rr,.9))  # moisture roughness
        else:
            zoq = minimum(1.15e-4,5.5e-5/power(rr,0.60))         # moisture roughness
        cdhf = von/(log(zu/zo)-velocity_structure_method_psiu_26(zu/L))
        cqhf = von/(log(zq/zoq)-temperature_structure_method_psit_26(zq/L))
        cthf = von/(log(zt/zot)-temperature_structure_method_psit_26(zt/L))
        usr = ut*cdhf
#         if isinstance(usr,ndarray)==False:
#             usr = ma.array([usr])
        qsr = -(dq-wetc*dter*jcool)*cqhf
        logging.debug('dq:%s wetc:% dter:%s jcool:%s cqhf:%s',dq,wetc,dter,jcool,cqhf)
        tsr = -(dt-dter*jcool)*cthf
        tvsr = tsr+0.61*ta*qsr
        tssr = tsr+0.51*ta*qsr
        Bf = -grav/ta*usr*tvsr
#         ug = 0.2*ones((N,1))
        ug = 0.2*ones(t.shape) #agrouaze matricial 
        k = where(Bf>0)
#         if hasattr(ug,'__size__'):
#             nb_elem = ug.size
#         else:
#             nb_elem = 1
#         if nb_elem>1:
        ug[k] = maximum(.2,Beta*power((Bf[k]*zi),0.333))
#         else:
#             if Bf>0:
#                 ug = maximum(.2,Beta*power((Bf*zi),0.333))
        ut = sqrt(power(du,2)+power(ug,2))
        gf = ut/du
        hsb = -rhoa*cpa*usr*tsr
        hlb = -rhoa*Le*usr*qsr
        qout = Rnl+hsb+hlb
        dels = Rns*(0.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4)))
        qcol = qout-dels
        alq = Al*qcol+be*hlb*cpw/Le
#         if isinstance(alq,ndarray)==False:
#             alq = ma.array([alq])
#         xlamx = 6.0*ones((N,1))
        xlamx = 6.0*ones(t.shape) #agrouaze matricial 
        if isinstance(xlamx,ndarray)==False:
            xlamx = ma.array([xlamx])
        tkt = minimum(0.01, xlamx*visw/(sqrt(rhoa/rhow)*usr))#attention rhoa/rhow negatif!!!-> fillvalue masked
#         tkt = tkt
#         k = where(alq>0)[0]
        k = alq>0
        logging.debug('coare4 | k (alq>0) indice = %s',k.shape)
        tmp = 6./power(1+power(bigc[k]*alq[k]/power(usr[k],4),0.75),0.333)
#         logging.debug('coare4 | tmp = %s',tmp.shape)
#         if tmp.size>0:
        xlamx[k] = tmp
        logging.debug('coare4 | tkt:%s xlamx:%s rhoa%s usr%s',tkt[k].shape,xlamx[k].shape,rhoa[k].shape,usr[k].shape)
        tkt[k] = xlamx[k]*visw/(sqrt(rhoa[k]/rhow)*usr[k])
#         tkt[k] = 6./power(1+power(bigc[k]*alq[k]/power(usr[k],4),0.75),0.333)*visw/(sqrt(rhoa[k]/rhow)*usr[k]) #long version
        dter = qcol*tkt/tcw
        dqer = wetc*dter
        Rnl = 0.97*(5.67e-8*power((ts-dter*jcool+tdk),4)-Rl) # update dter
        if i==0: # save first iteration solution for case of zetu>50;
            usr50 = usr[k50]
            tsr50 = tsr[k50]
            qsr50 = qsr[k50]
            L50 = L[k50]
            zet50 = zet[k50]
            dter50 = dter[k50]
            dqer50 = dqer[k50]
            tkt50 = tkt[k50]
        u10 = ut + usr/von*(log(10./zu)-velocity_structure_method_psiu_26(10./L)+velocity_structure_method_psiu_26(zu/L))
        charn = a1*u10+a2
        k = u10>umax
        charn[k] = a1*umax+a2
    
    # insert first iteration solution for case with zetu>50
    if k50 != []:
        usr[k50] = usr50
        tsr[k50] = tsr50
        qsr[k50] = qsr50
        L[k50] = L50
        zet[k50] = zet50
        dter[k50] = dter50
        dqer[k50] = dqer50
        tkt[k50] = tkt50
    
    #****************  compute fluxes  ********************************************
    tau = rhoa*usr*usr/gf      # wind stress
    hsb = -rhoa*cpa*usr*tsr     # sensible heat flux
    hlb = -rhoa*Le*usr*qsr      # latent heat flux
    hbb = rhoa*cpa*usr*tvsr    # buoyancy flux
    hsbb = rhoa*cpa*usr*tssr   # sonic heat flux
    logging.debug('rhoa:%s Le:%s usr:%s qsr:%s tsr:%s cpa:%s',rhoa,Le,usr,qsr,tsr,cpa)
    #*****  compute transfer coeffs relative to ut @ meas. ht  ********************
    Cd = tau/rhoa/ut/maximum(0.1,du)
    Ch = -usr*tsr/ut/(dt-dter*jcool)
    Ce = -usr*qsr/(dq-dqer*jcool)/ut
    
    #***  compute 10-m neutral coeff relative to ut (output if needed) ************
    Cdn_10 = 1000.*power(von,2)/power(log(10./zo),2)
    Chn_10 = 1000.*power(von,2)*fdg/log(10./zo)/log(10./zot)
    Cen_10 = 1000.*power(von,2)*fdg/log(10./zo)/log(10./zoq)
    
    #***  compute 10-m neutral coeff relative to ut (output if needed) ************
    #  Find the stability functions
    #*********************************
    zrf_u = 10.             #User defined reference heights
    zrf_t = 10.
    zrf_q = 10.
    psi = velocity_structure_method_psiu_26(zu/L)
    psi10 = velocity_structure_method_psiu_26(10./L)
    psirf = velocity_structure_method_psiu_26(zrf_u/L)
    psiT = temperature_structure_method_psit_26(zt/L)
    psi10T = temperature_structure_method_psit_26(10./L)
    logging.debug('coare4 | zrf_t:%s zrf_q:%s zt:%s',zrf_t,zrf_q,zt)
    psirfT = temperature_structure_method_psit_26(zrf_t/L)
    psirfQ = temperature_structure_method_psit_26(zrf_q/L)
    gf = ut/du
    
    #*********************************************************
    #  Determine the wind speeds relative to ocean surface
    #  Note that usr is the friction velocity that includes 
    #  gustiness usr = sqrt(Cd) S, which is equation (18) in
    #  Fairall et al. (1996)
    #*********************************************************
    S = ut
    U = du
    S10 = S + usr/von*(log(10./zu)-psi10+psi)
    U10 = S10/gf
    # or U10 = U + usr./von./gf.*(log(10/zu)-psi10+psi);
    Urf = U + usr/von/gf*(log(zrf_u/zu)-psirf+psi)
    UN = U + psi*usr/von/gf
    U10N = U10 + psi10*usr/von/gf
    UrfN = Urf + psirf*usr/von/gf
    
    UN2 = usr/von/gf*log(zu/zo)
    U10N2 = usr/von/gf*log(10./zo)
    UrfN2  = usr/von/gf*log(zrf_u/zo)
    
    lapse = grav/cpa
    SST = ts-dter*jcool
    
    T = t
    T10 = T + tsr/von*(log(10./zt)-psi10T+psiT) + lapse*(zt-10.)
    Trf = T + tsr/von*(log(zrf_t/zt)-psirfT+psiT) + lapse*(zt-zrf_t)
    TN = T + psiT*tsr/von
    T10N = T10 + psi10T*tsr/von
    TrfN = Trf + psirfT*tsr/von
    
    TN2 = SST + tsr/von*log(zt/zot)-lapse*zt
    T10N2 = SST + tsr/von*log(10./zot)-lapse*10.
    TrfN2 = SST + tsr/von*log(zrf_t/zot)-lapse*zrf_t
    
    dqer = wetc*dter*jcool
    SSQ = Qs-dqer
    SSQ = SSQ*1000
    Q = Q*1000
    qsr = qsr*1000
    Q10 = Q + qsr/von*(log(10./zq)-psi10T+psiT)
    Qrf = Q + qsr/von*(log(zrf_q/zq)-psirfQ+psiT)
    QN = Q + psiT*qsr/von/sqrt(gf)
    Q10N = Q10 + psi10T*qsr/von
    QrfN = Qrf + psirfQ*qsr/von
    
    QN2 = SSQ + qsr/von*log(zq/zoq)
    Q10N2 = SSQ + qsr/von*log(10./zoq)
    QrfN2 = SSQ + qsr/von*log(zrf_q/zoq)
    RHrf = relative_humidity_method_fairall(Trf,P,Qrf/1000)
    logging.info('coare4 | usr:%s tsr:%s',ma.median(usr),ma.median(tsr))
    #****************  output  ****************************************************
    
#     A = [usr, tau, hsb, hlb, hbb, hsbb, tsr, qsr, zot, zoq, Cd, Ch, Ce,  L, zet, dter, dqer, tkt, Urf, Trf, Qrf, RHrf, UrfN, Rnl, Le, rhoa, UN, U10, U10N, Cdn_10, Chn_10, Cen_10]
    #   1   2   3   4   5   6    7   8   9  10  11 12 13 14  15  16   17   18  19  20  21  22   23  24  25  26  27  28  29     30     31    32
    res = {}
    res['usr'] = usr
    res['tau'] = tau
    res['hsb'] = hsb
    res['hlb'] = hlb
    res['hbb'] = hbb
    res['hsbb'] = hsbb
    res['tsr'] = tsr
    res['qsr'] = qsr
    res['zot'] = zot
    res['zoq'] = zoq
    res['Cd'] = Cd
    res['Ch'] = Ch
    res['Ce'] = Ce
    res['L'] = L
    res['zet'] = zet
    res['dter'] = dter
    res['dqer'] = dqer
    res['tkt'] = tkt
    res['Urf'] = Urf
    res['Trf'] = Trf
    res['Qrf'] = Qrf
    res['RHrf'] = RHrf
    res['UrfN'] = UrfN
    res['Rnl'] = Rnl
    res['Le'] = Le
    res['rhoa'] = rhoa
    res['UN'] = UN
    res['U10'] = U10
    res['U10N'] = U10N
    res['Chn_10'] = Chn_10
    res['Cen_10'] = Cen_10
    res['Cdn_10'] = Cdn_10
    return res