"""
:author: Antoine Grouazel
@purpose: give a single input/output interface to call coare3 or coare4 methods
creation: 31/08/2015
"""
from numpy import ma,ndarray
import logging
import os
# from coare4_transcoding import coare40vn
# from coare3_transcoding import cor30a
# from relative_humidity import get_relative_humidity
def default_input_values(shape,param=None):
    """
    purpose: give a matrix set with a default value for a parameter you dont have for flux computation
    Args:
        shape (tuple): give the dimension of input field latxlon (e.g 720x1440)
        param (str): name of the variable you want
    return:
        x (list or dict): containing default values to be used in coare3/4
    """
    vals = {}
    pwd = os.path.dirname(__file__)
    fid = open(os.path.join(pwd,'list_default_values.txt'),'r')
    lines = fid.readlines()
    fid.close()
    for ll in lines:
        ll = ll.strip('\n')
        key = ll.split('=')[0]
        value = float(ll.split('=')[1])
        if key in ['zu','zt','zq','jcool','jwave','zi']:
            vals[key] = value
        elif key=='lat':
            vals[key] = value*ma.ones((shape[0],))
        else:
            vals[key] = value*ma.ones(shape)
#     vals['zu'] = 10.0
#     vals['zt'] = 10.0
#     vals['zq'] = 10.0
#     vals['u'] = vals['u']*ma.ones(shape)
#     vals['us'] = vals['us']*ma.ones(shape)
#     vals['ts'] = 20.0*ma.ones(shape)
#     vals['t'] = 25.0*ma.ones(shape)
#     vals['Qs'] = 18.0*ma.ones(shape)#g/kg
#     vals['Q'] = 15.0*ma.ones(shape)#g/kg
#     vals['Rs'] = 150.0*ma.ones(shape)# ben
#     vals['Rl'] = 370.0*ma.ones(shape)# ben
#     vals['zi'] = 600.0 #600.0 since Fairall
#     vals['rain'] = 0.*ma.ones(shape)
# #     vals['jcool'] = 0
# #     vals['jwave']  = 0
#     vals['twave'] = 6.0*ma.ones(shape)
#     vals['hwave'] = 0.0*ma.ones(shape)
#     vals['lat'] = 0.*ma.ones((shape[0],))
#     vals['P'] = 1008.*ma.ones(shape)
#     vals['rh'] = 80.0*ma.ones(shape)
    if param is None:
        x = vals
#         x = u,us,ts,t,Qs,Q,Rs,Rl,rain,zi,P,zu,zt,zq,lat,jcool,jwave,twave,hwave
    else:
        if param in vals.keys():
            x = vals[param]
        else:
            raise Exception('%s is not defined.'%param)
    return x

def complete_missing_args(list_input_needed,inputs_given,shape_wanted=None):
    """
    put default value when the given args are not fully-filled
    Args:
        list_input_needed (list):
        inputs_given (dic):
    Returns:
        inputs_complet (dic):
    """
    inputs_complet = {}
    shape_wanted = (1,)
    if shape_wanted is None:
        #find the shape of the fields
        max_length = 0
        for ii in inputs_given.keys():
            if isinstance(inputs_given[ii],ndarray)==False and ii not in ['jwave','jcool']:
                raise Exception('%s argument must be numpy masked array (even unique value) here: %s'%(ii,inputs_given[ii]))
            if isinstance(inputs_given[ii],ndarray):
                leni = len(inputs_given[ii].shape)
                if leni>max_length:
                        shape_wanted = inputs_given[ii].shape
                        max_length = leni
        logging.info('shape of input and output fields: %s',shape_wanted)
        
    #fill missing args
    for input_needed in list_input_needed:
        if input_needed in inputs_given.keys():
            if isinstance(inputs_given[input_needed],ndarray)==False:
                inputs_given[input_needed] = ma.array([inputs_given[input_needed]])
            inputs_complet[input_needed] = inputs_given[input_needed]
            logging.info('you have provided %s = %s',input_needed,inputs_complet[input_needed].shape)
            logging.debug('provided %s = %s',input_needed,inputs_complet[input_needed])
        else:
            inputs_complet[input_needed] = default_input_values(shape_wanted,input_needed)
            logging.info(' %s is default value',input_needed)
            logging.debug(' %s is default value: %s',input_needed,inputs_complet[input_needed])
            
    return inputs_complet

# def input_output_homogeneisor(version,shape_wanted,u=None,zu=None,t=None,zt=None,rh=None,zq=None,P=None,ts=None,Rs=None,Rl=None,lat=None,zi=None,us=None,Qs=None,Q=None,rain=None,jcool=None,jwave=None,twave=None,hwave=None):
#     """
#     purpose: wrapper to call coare4 or coare3 with same inputs/outputs
#     Args:
#         version (str): coare version 3 or 4
#         shape_wanted (tuple): shape of the inputs matrix (usualy (lon,lat))
#     Return:
#         res (dic): dictionary of float arrays containing the output parameter of COARE3 or COARE4
#     """
# #     def_input = default_input_values(shape_wanted)
# #     u,us,ts,t,Qs,Q,Rs,Rl,rain,zi,P,zu,zt,zq,lat,jcool,jwave,twave,hwave = def_input
# #     u_0,zu_0,t_0,zt_0,rh_0,zq_0,P_0,ts_0,Rs_0,Rl_0,lat_0,zi_0,us_0,Qs_0,Q_0,rain_0,jcool,jwave,twave,hwave = def_input
#     inputs = [u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,us,Qs,Q,rain,jcool,jwave,twave,hwave]
#     list_input = ['u','zu','t','zt','rh','zq','P','ts','Rs','Rl','lat','zi','us','Qs','Q','rain','jcool','jwave','twave','hwave']
#     inputs_dic = {}
#     #put default values for parameter not given
#     for ii in range(len(inputs)):
#         if inputs[ii] is None:
#             inputs_dic[list_input[ii]] = default_input_values(shape_wanted,list_input[ii])
#         else:
#             logging.info('you have provided %s',list_input[ii])
#             inputs_dic[list_input[ii]] = inputs[ii]
#     #special case RH
#     if Q is not None and t is not None:
#         inputs_dic['rh'] = get_relative_humidity(air_temperature=t,air_specific_humdity=Q/1000.0)*100.
#     res = {}
#     list_output = ['usr', 'tau', 'hsb', 'hlb', 'hbb', 'hsbb', 'tsr', 'qsr', 'zot', 'zoq', 'Cd', 'Ch', 'Ce', 'L', 'zet', 'dter', 'dqer', 'tkt', 'Urf', 'Trf', 'Qrf', 'RHrf', 'UrfN', 'Rnl', 'Le', 'rhoa', 'UN', 'U10', 'U10N', 'Cdn_10', 'Chn_10', 'Cen_10','zo','RF', 'wbar', 'ug' ]
#     for pp in list_output:
#         res[pp] = None#set to None all possible output
#     if version == '4':
#         A = coare40vn(u=inputs_dic['u'],zu=inputs_dic['zu'],t=inputs_dic['t'],zt=inputs_dic['zt'],rh=inputs_dic['rh'],zq=inputs_dic['zq'],P=inputs_dic['P'],ts=inputs_dic['ts'],Rs=inputs_dic['Rs'],Rl=inputs_dic['Rl'],lat=inputs_dic['lat'],zi=inputs_dic['zi'])
#         [usr, tau, hsb, hlb, hbb, hsbb, tsr, qsr, zot, zoq, Cd, Ch, Ce,  L, zet, dter, dqer, tkt, Urf, Trf, Qrf, RHrf, UrfN, Rnl, Le, rhoa, UN, U10, U10N, Cdn_10, Chn_10, Cen_10] = A
#         res['usr'] = usr
#         res['tau'] = tau
#         res['hsb'] = hsb
#         res['hlb'] = hlb
#         res['hbb'] = hbb
#         res['hsbb'] = hsbb
#         res['tsr'] = tsr
#         res['qsr'] = qsr
#         res['zot'] = zot
#         res['zoq'] = zoq
#         res['Cd'] = Cd
#         res['Ch'] = Ch
#         res['Ce'] = Ce
#         res['L'] = L
#         res['zet'] = zet
#         res['dter'] = dter
#         res['dqer'] = dqer
#         res['tkt'] = tkt
#         res['Urf'] = Urf
#         res['Trf'] = Trf
#         res['Qrf'] = Qrf
#         res['RHrf'] = RHrf
#         res['UrfN'] = UrfN
#         res['Rnl'] = Rnl
#         res['Le'] = Le
#         res['rhoa'] = rhoa
#         res['UN'] = UN
#         res['U10'] = U10
#         res['U10N'] = U10N
#         res['Chn_10'] = Chn_10
#         res['Cen_10'] = Cen_10
#         res['Cdn_10'] = Cdn_10
#     elif version == '3':
# #         x = u,us,ts,t,Qs,Q,Rs,Rl,rain,zi,P,zu,zt,zq,lat,jcool,jwave,twave,hwave
#         y = cor30a(u=inputs_dic['u'],us=inputs_dic['us'],ts=inputs_dic['ts'],t=inputs_dic['t'],Qs=inputs_dic['Qs'],Q=inputs_dic['Q'],Rs=inputs_dic['Rs'],Rl=inputs_dic['Rl'],rain=inputs_dic['rain'],zi=inputs_dic['zi'],P=inputs_dic['P'],zu=inputs_dic['zu'],zt=inputs_dic['zt'],zq=inputs_dic['zq'],lat=inputs_dic['lat'],jcool=inputs_dic['jcool'],jwave=inputs_dic['jwave'],twave=inputs_dic['twave'],hwave=inputs_dic['hwave'])
#         [hsb, hlb, tau, zo, zot, zoq, L, usr, tsr, qsr, dter, dqer, tkt, RF, wbar, Cd, Ch, Ce, Cdn_10, Chn_10, Cen_10, ug ] = y
#         res['usr'] = usr
#         res['tau'] = tau
#         res['hlb'] = hlb
#         res['hsb'] = hsb
#         res['tsr'] = tsr
#         res['qsr'] = qsr
#         res['zot'] = zot
#         res['zoq'] = zoq
#         res['Cd'] = Cd
#         res['Ch'] = Ch
#         res['Ce'] = Ce
#         res['L'] = L
#         res['dter'] = dter
#         res['dqer'] = dqer
#         res['tkt'] = tkt
#         res['Chn_10'] = Chn_10
#         res['Cen_10'] = Cen_10
#         res['Cdn_10'] = Cdn_10
#         res['zo'] = zo
#         res['RF'] = RF
#         res['wbar'] = wbar
#         res['ug'] = ug
# #     v3 = [zo,RF, wbar, ug ]
# #     out = [usr, tau, hsb, hlb, hbb, hsbb, tsr, qsr, zot, zoq, Cd, Ch, Ce,  L, zet, dter, dqer, tkt, Urf, Trf, Qrf, RHrf, UrfN, Rnl, Le, rhoa, UN, U10, U10N, Cdn_10, Chn_10, Cen_10]+v3
#     return res