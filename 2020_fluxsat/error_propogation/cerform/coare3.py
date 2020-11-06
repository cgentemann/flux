# encoding: utf-8
"""
.. module::cerform.flux.coare3

COARE 3.0 bulk algorithm

:copyright: Copyright 2015 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Antoine Grouazel 

@creation: 23/06/2015
@purpose: compute FLux parameters
@context: OHF project
@note: results are ok regarding what is provided here: ftp://ftp1.esrl.noaa.gov/users/cfairall/wcrp_wgsf/computer_programs/cor3_0/coare30a_readme_1.pdf
@todo: validation with references OHF datasets
"""
from numpy import log,exp,sqrt,pi,minimum,power,maximum,ma,ones,amax,amin,ndarray,array
import logging
import pdb
from temperature_structure import temperature_structure_method_psit_30
from velocity_structure import velocity_structure_method_psiu_30
from gravity_constant import grv
from default_values import complete_missing_args

__all__ = ['coare3']

def coare3(inputs):
    """
    :disclaimer: 
        vectorial computing not yet validated
    :history: 
        transcoding from cor30a.m Fairall et al
        version with shortened iteration modified Rt and Rq
        uses wave information wave period in s and wave ht in m
        no wave, standard coare 2.6 charnock:  jwave=0 
        Oost et al.  zo=50/2/pi L (u*/c)^4.5 if jwave=1
        taylor and yelland  zo=1200 h*(L/h)^4.5 jwave=2
        x=[5.5 0 28.7 27.2 24.2 18.5 141 419 0 600 1010 15 15 15 0 1 1 5 1 ]
    
    Args:
        inputs (dic): inputs parameter containing u,us,ts,t,Qs,Q,Rs,Rl,rain,zi,P,zu,zt,zq,lat,jcool,jwave,twave,hwave fields
            u (float): wind speed (m/s)  at height zu (m)
            us (float): surface current speed in the wind direction (m/s)
            ts (float): bulk water temperature (C) if jcool=1, interface water T if jcool=0  
            t (float): bulk air temperature (C), height zt
            Qs (float): bulk water spec hum (g/kg) if jcool=1, ...
            Q (float): bulk air spec hum (g/kg), height zq
            Rs (float): downward solar flux (W/m^2)
            Rl (float): downard IR flux (W/m^2)
            rain (float): rain rate (mm/hr)
            zi (float): Planet Boundary Layer depth (m)
            P (float): Atmos surface pressure (mb)
            zu (float): wind speed measurement height (m)
            zt (float): air T measurement height (m)
            zq (float): air q measurement height (m)
            lat (float): latitude (deg, N=+)
            jcool (float): implement cool calculation skin switch, 0=no, 1=yes
            jwave (float): implement wave dependent roughness model
            twave (float): wave period (s)
            hwave (float): wave height (m)
    
    Hint:
        !! MIND THE CASE of the inputs keys !!
    
    Returns:
        A dict containing the following keys:
        
        {hsb (float) :    sensible heat flux (w/m^2),
        hlb (float) :    latent heat flux (w/m^2),
        RF (float) :     rain heat flux(w/m^2),
        wbar (float) :   webb mean w (m/s),
        tau (float) :    stress (nt/m^2),
        zo (float) :      velocity roughness length (m),
        zot (float) :      temperature roughness length (m),
        zoq (float) :      moisture roughness length (m),
        L (float) :        Monin_Obukhov stability length,
        usr (float) :      turbulent friction velocity (m/s), including gustiness,
        tsr (float) :      temperature scaling parameter (K),
        qsr (float) :      humidity scaling parameter (g/g),
        dter (float) :     cool skin temperature depression (K),
        dqer (float) :     cool skin humidity depression (g/g),
        tkt (float) :      cool skin thickness (m),
        Cd (float) :       velocity drag coefficient at zu, referenced to u,
        Ch (float) :       heat transfer coefficient at zt,
        Ce (float) :       moisture transfer coefficient at zq,
        Cdn_10 (float) :   10-m velocity drag coeeficient, including gustiness,
        Chn_10 (float) :   10-m heat transfer coeeficient, including gustiness,
        Cen_10 (float) :   10-m humidity transfer coeeficient, including gustiness,
        ug (float) :       geostrophic wind [m.s-1],}
        
    Warning: 
        vectorized version
    """
    list_input_needed = ['u','us','ts','t','Qs','Q','Rs','Rl','rain','zi','P','zu','zt','zq','lat','jcool','jwave','twave','hwave']
    inputs_complet = complete_missing_args(list_input_needed,inputs)
    u=inputs_complet['u']
    us=inputs_complet['us']
    ts=inputs_complet['ts']
    t=inputs_complet['t']
    Qs=inputs_complet['Qs']
    Q=inputs_complet['Q']
    Rs=inputs_complet['Rs']
    Rl=inputs_complet['Rl']
    rain=inputs_complet['rain']
    zi=inputs_complet['zi']
    P=inputs_complet['P']
    zu=inputs_complet['zu']
    zt=inputs_complet['zt']
    zq=inputs_complet['zq']
    lat=inputs_complet['lat']
    jcool=inputs_complet['jcool']
    jwave=inputs_complet['jwave']
    twave=inputs_complet['twave']
    hwave=inputs_complet['hwave']
    
    
    
    Qs = Qs/1000. #in original code: meant to convert g/kg -> mg/kg 
    Q = Q/1000.0 #the rest of the code need to use mg/kg spec humidity 
    logging.info('coare3 | Q:%s Qs:%s',ma.median(Q),ma.median(Qs))
    logging.debug('windspeed:%s currentspeed:%s sst:%s airt:%s water_spec_hum:%s \
    air_spec_hum:%s solar:%s IR:%s rain:%s P:%s levelWS:%s hs:%s',u,us,ts,t,Qs,Q,Rs,Rl,rain,P,zu,hwave)
    #############   set constants ############
    Beta = 1.2
    von = 0.4
    fdg = 1.00
    tdk = 273.16
#    print(lat.size())
    grav = grv(lat) #,shape_wanted=u.shape)#9.82
    #############  air constants ############
    Rgas = 287.1
    Le = (2.501-0.00237*ts)*1e6
    cpa = 1004.67
#     cpv = cpa*(1+0.84*Q) #unused
    rhoa = P*100./(Rgas*(t+tdk)*(1.+0.61*Q))
    logging.debug('rhoa:%s %s ',amax(rhoa),amin(rhoa))
    visa = 1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t)
    #############*  cool skin constants  ############
    Al = 2.1e-5*(ts+3.2)**0.79
    be = 0.026
    cpw = 4000.
    rhow = 1022.
    visw = 1e-6
    tcw = 0.6
#     logging.info('coare3 | shapes grav:%s %s',grav.shape,grav)
    bigc = 16*grav*cpw*power((rhow*visw),3)/(tcw*tcw*rhoa*rhoa)
    wetc = 0.622*Le*Qs/(Rgas*power((ts+tdk),2))

    #############   wave parameters  ############
    lwave = grav/2./pi*power(twave,2)
    cwave = grav/2./pi*twave
    
    #############  compute aux stuff ############
    Rns = Rs*.945
    Rnl = 0.97*(5.67e-8*power((ts-0.3*jcool+tdk),4)-Rl)
    
    #############   Begin bulk loop ############*
    
    #############  first guess ############
    du = u-us
    dt = ts-t-0.0098*zt
    dq = Qs-Q
    ta = t+tdk
    ug = 0.5
    dter = 0.3
    dqer = wetc*dter
    ut = sqrt(du*du+ug*ug)
    u10 = ut*log(10/1e-4)/log(zu/1e-4)
    usr = 0.035*u10
    zo10 = 0.011*usr*usr/grav+0.11*visa/usr
    Cd10 = power((von/log(10/zo10)),2)
    Ch10 = 0.00115
    Ct10 = Ch10/sqrt(Cd10)
    zot10 = 10.0/exp(von/Ct10)
    Cd = power((von/log(zu/zo10)),2)
    Ct = von/log(zt/zot10)
    CC = von*Ct/Cd
    Ribcu = -zu/zi/0.004/power(Beta,3)
    Ribu = -grav*zu/ta*((dt-dter*jcool)+0.61*ta*dq)/power(ut,2)
    nits = 3

    zetu = CC*Ribu/(1+Ribu/Ribcu)
    logging.debug('coare3 | zetu: %s %s',zetu,isinstance(zetu,(int,float,complex)))
    zetu[Ribu>=0] = CC[Ribu>=0]*Ribu[Ribu>=0]*(1+27.0/9.0*Ribu[Ribu>=0]/CC[Ribu>=0])
    
    
    L10 = zu/zetu
    logging.debug('coar3 | zetu:%s',zetu)
    if (zetu>50).any():
            nits = 1
    logging.debug('cor30a | nber of loop:%s',nits)
    
    usr = ut*von/(log(zu/zo10)-velocity_structure_method_psiu_30(zu/L10))
    tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot10)-temperature_structure_method_psit_30(zt/L10))
    qsr = -(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-temperature_structure_method_psit_30(zq/L10))
    
    tkt = 0.001
    charn = 0.011*ones(ut.shape)

    charn[ut>10] = 0.011+(ut[ut>10]-10)/(18-10)*(0.018-0.011)
    charn[ut>18] = 0.018
    
    #############  bulk loop ############
    for i in range(nits):
        zet = von*grav*zu/ta*(tsr*(1+0.61*Q)+0.61*ta*qsr)/(usr*usr)/(1+0.61*Q)
        if jwave==0:
            zo = charn*usr*usr/grav+0.11*visa/usr
        if jwave==1:
            zo = 50./2./pi*lwave*(usr/cwave)**4.5+0.11*visa/usr #Oost et al
        if jwave==2:
            zo = 1200.*hwave*(hwave/lwave)**4.5+0.11*visa/usr #Taylor and Yelland
        rr = zo*usr/visa
        L = zu/zet
        zoq = minimum(1.15e-4,5.5e-5/(rr**0.6))
        zot = zoq
        
        usr = ut*von/(log(zu/zo)-velocity_structure_method_psiu_30(zu/L))
        
        tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot)-temperature_structure_method_psit_30(zt/L))
        
        qsr = -(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-temperature_structure_method_psit_30(zq/L))
        
        logging.debug('dq:%s wetc:%s dter:%s jcool:%s cqhf:%s',dq,wetc,dter,jcool,von*fdg/(log(zq/zoq)-temperature_structure_method_psit_30(zq/L)))
        Bf = -grav/ta*usr*(tsr+.61*ta*qsr)
        
        ug = Beta*(Bf*zi)**0.333
        ug[Bf<=0] = 0.2
        ut = sqrt(du*du+ug*ug)
        Rnl = 0.97*(5.67e-8*power((ts-dter*jcool+tdk),4)-Rl)
        hsb = -rhoa*cpa*usr*tsr
        hlb = -rhoa*Le*usr*qsr
        qout = Rnl+hsb+hlb
        dels = Rns*(.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4)))    # Eq.16 Shortwave
        qcol = qout-dels
        alq = Al*qcol+be*hlb*cpw/Le                                 # Eq. 7 Buoy flux water
        
        xlamx = 6./(1+(bigc*alq/power(usr,4))**0.75)**0.333    # Eq 13 Saunders
        tkt = xlamx*visw/(sqrt(rhoa/rhow)*usr)                       #Eq.11 Sub. thk
        xlamx[alq<=0] = array([6.0])
        logging.debug('coare3 | xlamx:%s visw:%s rhoa:%s rhow:%s usr:%s',xlamx.shape,visw,rhoa[alq<=0].shape,rhow,usr[alq<=0].shape)
        tkt[alq<=0] = minimum(0.01,xlamx[alq<=0]*visw/(sqrt(rhoa[alq<=0]/rhow)*usr[alq<=0]))                      #Eq.11 Sub. thk
        dter = qcol*tkt/tcw#  Eq.12 Cool skin
        dqer = wetc*dter
    tau = rhoa*usr*usr*du/ut                #stress
    hsb = -rhoa*cpa*usr*tsr
    hlb = -rhoa*Le*usr*qsr
    logging.debug('rhoa:%s Le:%s usr:%s qsr:%s tsr:%s cpa:%s',rhoa,Le,usr,qsr,tsr,cpa)
    #############   rain heat flux ############
    
    dwat = 2.11e-5*((t+tdk)/tdk)**1.94 #! water vapour diffusivity
    dtmp = (1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa)      #!heat diffusivity
    alfac= 1./(1.+(wetc*Le*dwat)/(cpa*dtmp))           #! wet bulb factor
    RF = rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600.0
    #############   Webb et al. correection  ############
    wbar = 1.61*hlb/Le/(1+1.61*Q)/rhoa+hsb/rhoa/cpa/ta#formulation in hlb already includes webb
    #wbar=1.61*hlb/Le/rhoa+(1+1.61*Q)*hsb/rhoa/cpa/ta
    hl_webb = rhoa*wbar*Q*Le
    #############   compute transfer coeffs relative to ut @meas. ht ############
    Cd = tau/rhoa/ut/maximum(.1,du)
    Ch = -usr*tsr/ut/(dt-dter*jcool)
    Ce = -usr*qsr/(dq-dqer*jcool)/ut
    #############  10-m neutral coeff realtive to ut ############
    Cdn_10 = von*von/log(10./zo)/log(10./zo)
    Chn_10 = von*von*fdg/log(10./zo)/log(10./zot)
    Cen_10 = von*von*fdg/log(10./zo)/log(10./zoq)
    
    logging.info('hsb:%s hlb:%s',ma.median(hsb),ma.median(hlb))
    logging.debug('usr: %s tsr:%s',ma.median(usr),ma.median(tsr))
#     y = [hsb, hlb, tau, zo, zot, zoq, L, usr, tsr, qsr, dter, dqer, tkt, RF, wbar, Cd, Ch, Ce, Cdn_10, Chn_10, Cen_10, ug ]
    res = {}
    res['hsb'] = hsb
    res['hlb'] = hlb
    res['tau'] = tau
    res['zo'] = zo
    res['zot'] = zot
    res['zoq'] = zoq
    res['L'] = L
    res['usr'] = usr
    res['tsr'] = tsr
    res['qsr'] = qsr
    res['dter'] = dter
    res['dqer'] = dqer
    res['tkt'] = tkt
    res['RF'] = RF
    res['wbar'] = wbar
    res['Cd'] = Cd
    res['Ch'] = Ch
    res['Ce'] = Ce
    res['Cdn_10'] = Cdn_10
    res['Chn_10'] = Chn_10
    res['Cen_10'] = Cen_10
    res['ug'] = ug
    return res




    











