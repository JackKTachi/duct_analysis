import pyspedas
import pytplot 
from pytplot.tplot import tplot

# xarrayにはmepeのdatatype='pa'のFEDUデータをdata_quantsでxarrayに取り出してから入れる

def file_name():
    return "'erg_mepe_l3_pa_FEDU_87.5keV', 'erg_mepe_l3_pa_FEDU_72.6keV', 'erg_mepe_l3_pa_FEDU_60.4keV', 'erg_mepe_l3_pa_FEDU_50.3keV', 'erg_mepe_l3_pa_FEDU_42.0keV','erg_mepe_l3_pa_FEDU_35.0keV','erg_mepe_l3_pa_FEDU_29.3keV','erg_mepe_l3_pa_FEDU_24.5keV','erg_mepe_l3_pa_FEDU_20.5keV','erg_mepe_l3_pa_FEDU_17.1keV','erg_mepe_l3_pa_FEDU_14.3keV','erg_mepe_l3_pa_FEDU_12.0keV','erg_mepe_l3_pa_FEDU_10.0keV','erg_mepe_l3_pa_FEDU_8.4keV','erg_mepe_l3_pa_FEDU_7.0keV'"

def each_eV_plot(xarray):
    new_coords = {'v2_dim': [5.,  15.,  25.,  35.,  45.,  55.,  65.,  75.,  85.,  95., 105.,
       115., 125., 135., 145., 155., 165., 175.]}
    xarray = xarray.assign_coords(new_coords)
    data1 = xarray
    pytplot.store_data('erg_mepe_l3_pa_FEDU_87.5keV', data = {'x':data1['time'], 'y':data1[:,1,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_87.5keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'87.5keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_72.6keV', data = {'x':data1['time'], 'y':data1[:,2,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_72.6keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'72.6keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_60.4keV', data = {'x':data1['time'], 'y':data1[:,3,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_60.4keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'60.4keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_50.3keV', data = {'x':data1['time'], 'y':data1[:,4,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_50.3keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'50.3keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_42.0keV', data = {'x':data1['time'], 'y':data1[:,5,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_42.0keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'42.0keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_35.0keV', data = {'x':data1['time'], 'y':data1[:,6,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_35.0keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'35.0keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_29.3keV', data = {'x':data1['time'], 'y':data1[:,7,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_29.3keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'29.3keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_24.5keV', data = {'x':data1['time'], 'y':data1[:,8,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_24.5keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'24.5keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_20.5keV', data = {'x':data1['time'], 'y':data1[:,9,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_20.5keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'20.5keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_17.1keV', data = {'x':data1['time'], 'y':data1[:,10,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_17.1keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'17.1keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_14.3keV', data = {'x':data1['time'], 'y':data1[:,11,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_14.3keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'14.3keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_12.0keV', data = {'x':data1['time'], 'y':data1[:,12,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_12.0keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'12.0keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_10.0keV', data = {'x':data1['time'], 'y':data1[:,13,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_10.0keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'10.0keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_8.4keV', data = {'x':data1['time'], 'y':data1[:,14,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_8.4keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'8.4keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    pytplot.store_data('erg_mepe_l3_pa_FEDU_7.0keV', data = {'x':data1['time'], 'y':data1[:,15,:], 'v':data1['v2_dim']})
    pytplot.options('erg_mepe_l3_pa_FEDU_7.0keV', opt_dict={'spec':1,'zlog':1, 'ytitle':'7.0keV', 'ysubtitle':'PA[deg]', 'ztitle':'','zsubtitle':'[/s/cm^2/sr/keV]'})
    return 'erg_mepe_l3_pa_FEDU_87.5keV', 'erg_mepe_l3_pa_FEDU_72.6keV', 'erg_mepe_l3_pa_FEDU_60.4keV', 'erg_mepe_l3_pa_FEDU_50.3keV', 'erg_mepe_l3_pa_FEDU_42.0keV','erg_mepe_l3_pa_FEDU_35.0keV','erg_mepe_l3_pa_FEDU_29.3keV','erg_mepe_l3_pa_FEDU_24.5keV','erg_mepe_l3_pa_FEDU_20.5keV','erg_mepe_l3_pa_FEDU_17.1keV','erg_mepe_l3_pa_FEDU_14.3keV','erg_mepe_l3_pa_FEDU_12.0keV','erg_mepe_l3_pa_FEDU_10.0keV','erg_mepe_l3_pa_FEDU_8.4keV','erg_mepe_l3_pa_FEDU_7.0keV'

