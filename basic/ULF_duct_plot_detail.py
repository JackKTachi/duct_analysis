import pyspedas
import pytplot 
from pytplot.tplot import tplot
import numpy as np
import each_eV_plot_mepe as eep
import compressional_B as cb

def ULF_duct_load(download_trange): 


    lst = download_trange[0].split('-')
    pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/hfa/l3/1min/'+lst[0]+'/'+lst[1]+'/erg_pwe_hfa_l3_1min_'+lst[0]+lst[1]+lst[2]+'_v03_07.cdf')    
        
    from pyspedas.erg import pwe_ofa
    pwe_ofa(trange=download_trange)

    pytplot.cdf_to_tplot('./erg_data/satellite/erg/mepe/l3/pa/'+lst[0]+'/'+lst[1]+'/erg_mepe_l3_pa_'+lst[0]+lst[1]+lst[2]+'_v01_01.cdf')
    data1=pytplot.data_quants['FEDU']
    'erg_mepe_l3_pa_FEDU_87.5keV', 'erg_mepe_l3_pa_FEDU_72.6keV', 'erg_mepe_l3_pa_FEDU_60.4keV', 'erg_mepe_l3_pa_FEDU_50.3keV', 'erg_mepe_l3_pa_FEDU_42.0keV','erg_mepe_l3_pa_FEDU_35.0keV','erg_mepe_l3_pa_FEDU_29.3keV','erg_mepe_l3_pa_FEDU_24.5keV','erg_mepe_l3_pa_FEDU_20.5keV','erg_mepe_l3_pa_FEDU_17.1keV','erg_mepe_l3_pa_FEDU_14.3keV','erg_mepe_l3_pa_FEDU_12.0keV','erg_mepe_l3_pa_FEDU_10.0keV','erg_mepe_l3_pa_FEDU_8.4keV','erg_mepe_l3_pa_FEDU_7.0keV' == eep.each_eV_plot(data1)
    
    pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/ofa/l3/property/'+lst[0]+'/'+lst[1]+'/erg_pwe_ofa_l3_property_dsi_'+lst[0]+lst[1]+lst[2]+'_v01_03.cdf')
    #pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/ofa/l3/property/2017/03/erg_pwe_ofa_l3_property_dsi_20170329_v01_03.cdf')
    # 2017Mar30 is the file of v01_04, others are v01_03

    'B_perp_para', 'B1_total', 'delta_z', 'delta_xy', 'ULF_wna', 'erg_mgf_l2_magt_8sec', 'cos_th_intpl' == cb.ULFwna(download_trange)
    #'delta_z', 'erg_mgf_l2_mag_8sec_MAF_x&y', 'erg_mgf_l2_magt_8sec' == cb.ULFwna(download_trange)

    from pyspedas.erg import mepe
    mepe(trange=download_trange)

    from pyspedas.erg import orb
    orb(trange=download_trange)
    
    return 'ne_mgf',\
        'kvec_polar_132',\
        'erg_pwe_ofa_l2_spec_B_spectra_132','erg_pwe_ofa_l2_spec_E_spectra_132',\
        'erg_mepe_l2_omniflux_FEDO',\
        'erg_mepe_l3_pa_FEDU_87.5keV', 'erg_mepe_l3_pa_FEDU_72.6keV', 'erg_mepe_l3_pa_FEDU_60.4keV','erg_mepe_l3_pa_FEDU_50.3keV', 'erg_mepe_l3_pa_FEDU_42.0keV','erg_mepe_l3_pa_FEDU_35.0keV','erg_mepe_l3_pa_FEDU_29.3keV','erg_mepe_l3_pa_FEDU_24.5keV','erg_mepe_l3_pa_FEDU_20.5keV','erg_mepe_l3_pa_FEDU_17.1keV','erg_mepe_l3_pa_FEDU_14.3keV','erg_mepe_l3_pa_FEDU_12.0keV',\
            'B_perp_para', 'B1_total', 'delta_z', 'delta_xy', 'ULF_wna', 'erg_mgf_l2_magt_8sec',\
            'erg_orb_l2_pos_rmlatmlt','erg_orb_l2_pos_eq'

def ULF_evensearch_load(download_trange): 


    lst = download_trange[0].split('-')
        
    from pyspedas.erg import pwe_ofa
    pwe_ofa(trange=download_trange)

    pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/ofa/l3/property/'+lst[0]+'/'+lst[1]+'/erg_pwe_ofa_l3_property_dsi_'+lst[0]+lst[1]+lst[2]+'_v01_03.cdf')
    #pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/ofa/l3/property/2017/03/erg_pwe_ofa_l3_property_dsi_20170329_v01_03.cdf')
    # 2017Mar30 is the file of v01_04, others are v01_03

    'delta_z', 'erg_mgf_l2_mag_8sec_MAF_x&y', 'erg_mgf_l2_magt_8sec' == cb.compressional_B(download_trange)

    from pyspedas.erg import orb
    orb(trange=download_trange)
    
    return 'kvec_polar_132',\
        'erg_pwe_ofa_l2_spec_B_spectra_132','erg_pwe_ofa_l2_spec_E_spectra_132',\
            'delta_z', 'erg_mgf_l2_magt_8sec','erg_mgf_l2_mag_8sec_MAF_x&y',\
            'erg_orb_l2_pos_rmlatmlt','erg_orb_l2_pos_eq'

def ULF_eventsearch_plot(wna,ofa_B,deltaB,magxy,magt, orb_rmlatmlt, plot_trange):
    
    pytplot.options(wna, opt_dict={'ytitle':'wna','ysubtitle':'[kHz]','ylog':1, 'spec':1})

    pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    labels = pytplot.split_vec( orb_rmlatmlt )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
    pytplot.xlim( plot_trange[0], plot_trange[1] )

    pytplot.tplot( [wna,ofa_B,deltaB,magt], var_label=labels, xsize=10, ysize=10)

def ULF_duct_plot(density,wna,ofa_B,ofa_E,integral,deltaB,magxy,magt,mepe_50,mepe_35,mepe_24,mepe_12, orb_rmlatmlt, plot_trange):
    
    pytplot.options(wna, opt_dict={'ytitle':'wna','ysubtitle':'[kHz]','ylog':1, 'spec':1})

    pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options(ofa_E, opt_dict={'ytitle':'ofa-E','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options(integral, opt_dict={'ytitle':'integral','ylog':1})
    labels = pytplot.split_vec( orb_rmlatmlt )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
    pytplot.xlim( plot_trange[0], plot_trange[1] )

    pytplot.tplot( [density,wna,ofa_B,ofa_E,integral,deltaB,magxy,magt,mepe_50,mepe_35,mepe_24,mepe_12], var_label=labels, xsize=10, ysize=15)


def plot_B_int_z_xy(ofa_B, integral, z, xy, orb_rmlatmlt, plot_trange):

    pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options(integral, opt_dict={'ytitle':'integral','ylog':1})
    labels = pytplot.split_vec( orb_rmlatmlt )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
    pytplot.xlim( plot_trange[0], plot_trange[1] )

    pytplot.tplot( [ofa_B,integral,z,xy], var_label=labels, xsize=10, ysize=15)        


def plot_B_int_density(ofa_B,density, integral, z, xy, orb_rmlatmlt, plot_trange):

    pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options(integral, opt_dict={'ytitle':'integral','ylog':1})
    labels = pytplot.split_vec( orb_rmlatmlt )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
    pytplot.xlim( plot_trange[0], plot_trange[1] )

    pytplot.tplot( [ofa_B,density,integral,z,xy], var_label=labels, xsize=10, ysize=15)  


def kvec_Bspec_intgrl_rein_delB_magt_Bxy( wna, ofa_B, integral, deltaB, magt, magxy, orb_rmlatmlt, plot_trange, download_trange,f, theta):

    import plasma_params as pp

    from pyspedas.erg import mgf
    mgf(trange=download_trange)

    #lst = download_trange[0].split('-')
    #pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/hfa/l3/1min/'+lst[0]+'/'+lst[1]+'/erg_pwe_hfa_l3_1min_'+lst[0]+lst[1]+lst[2]+'_v03_07.cdf')
    #pyspedas.tinterpol('ne_mgf', 'erg_mgf_l2_magt_8sec', newname='ne_mgf_interpolated')
    B = pytplot.data_quants['erg_mgf_l2_magt_8sec']*1e-9
    #N = pytplot.data_quants['ne_mgf_interpolated']*1e6
    fc = (pp.Q*B*f/pp.ME)**(1/2)/2/np.pi
    #fp = (80.6*N)**(1/2)/2/np.pi
    #re_in = fp/(fc*f*(np.cos(theta)-f/fc))**0.5  

    N0 = 2e6

    fc = (pp.Q*B/pp.ME)/2/np.pi
    fp = ((pp.Q**2*N0)/(pp.EPS*pp.ME))**(1/2)/2/np.pi
    mu = fp*(fc*f*(np.cos(theta)-(f/fc)))**(-1/2)

    pytplot.store_data('re_in_', data={'x': pytplot.data_quants['erg_mgf_l2_magt_8sec']['time'], 'y': mu})
    pytplot.options('re_in_', opt_dict={'ytitle':'$\mu (N0)$','ylog':0})    
    #pytplot.store_data('re_in', data={'x': pytplot.data_quants['ne_mgf_interpolated']['time'], 'y': re_in})
    #pytplot.options('re_in', opt_dict={'ytitle':'$\mu (Ntime)$','ylog':0})



    pytplot.options(wna, opt_dict={'ytitle':'wna','ysubtitle':'[kHz]','ylog':1, 'spec':1})
    pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options(integral, opt_dict={'ytitle':'integral','ylog':1})
    #pytplot.options(re_in, opt_dict={'ytitle':'Refractive Index','ylog':0})
    pytplot.options(deltaB, opt_dict={'ytitle':'$\delta B$','ylog':0})
    pytplot.options(magt, opt_dict={'ytitle':'$B_{total}$','ylog':0})
    pytplot.options(magxy, opt_dict={'ytitle':'$B_{xy}$','ylog':0})
    labels = pytplot.split_vec( orb_rmlatmlt )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
    pytplot.options('ne_mgf', opt_dict={'ytitle':'density','ylog':0})
    pytplot.xlim( plot_trange[0], plot_trange[1] )

    pytplot.tplot( ['ne_mgf',wna,ofa_B,integral,'re_in_',deltaB, magt,magxy], var_label=labels, xsize=10, ysize=15)


def RIset(wna, ofa_B, integral, deltaB, magt, magxy, orb_rmlatmlt,path, plot_trange, download_trange, f):

    import plasma_params as pp
    import get_density as gd
    import calc_RI as cr
    import pyspedas

    from pyspedas.erg import mgf
    mgf(trange=download_trange)


    'ne' == gd.get_density(path)

    pyspedas.tinterpol('ne', 'erg_mgf_l2_magt_8sec', newname='ne_intpl')

    N = pytplot.data_quants['ne_intpl']*1e6
    B = pytplot.data_quants['erg_mgf_l2_magt_8sec']*1e-9

    mu = cr.calc_RI(N, B, f)
    pytplot.store_data('mu', data={'x':B['time'], 'y':mu})


    pytplot.options(wna, opt_dict={'ytitle':'wna','ysubtitle':'[kHz]','ylog':1, 'spec':1})
    pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options(integral, opt_dict={'ytitle':'integral','ylog':1})
    pytplot.options(deltaB, opt_dict={'ytitle':'$\delta B$','ylog':0})
    pytplot.options(magt, opt_dict={'ytitle':'$B_{total}$','ylog':0})
    pytplot.options(magxy, opt_dict={'ytitle':'$B_{xy}$','ylog':0})
    labels = pytplot.split_vec( orb_rmlatmlt )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
    pytplot.options('ne_mgf', opt_dict={'ytitle':'density','ylog':0})
    pytplot.xlim( plot_trange[0], plot_trange[1] )

    pytplot.tplot( [wna,ofa_B,integral,'mu', deltaB, 'ne', magt,magxy], var_label=labels, xsize=10, ysize=15)





""" 
    if density is not None:
        B = pytplot.data_quants['erg_mgf_l2_magt_8sec']*1e-9
        N = density
        fc = (pp.Q*B*f/pp.ME)**(1/2)
        fp = fp = (80.6*N)**(1/2)
        re_in = fp/(fc*f*(np.cos(theta)-f/fc))**0.5

        pytplot.store_data('re_in', data={'x': pytplot.data_quants['ne_mgf_interpolated']['time'], 'y': re_in})
        pytplot.options('re_in', opt_dict={'ytitle':'$\mu (N0)$','ylog':0})
        pytplot.options(wna, opt_dict={'ytitle':'wna','ysubtitle':'[kHz]','ylog':1, 'spec':1})
        pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
        pytplot.options(integral, opt_dict={'ytitle':'integral','ylog':1})
        #pytplot.options(re_in, opt_dict={'ytitle':'Refractive Index','ylog':0})
        pytplot.options(deltaB, opt_dict={'ytitle':'$\delta B$','ylog':0})
        pytplot.options(magt, opt_dict={'ytitle':'$B_{total}$','ylog':0})
        pytplot.options(magxy, opt_dict={'ytitle':'$B_{xy}$','ylog':0})
        labels = pytplot.split_vec( orb_rmlatmlt )
        pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
        pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
        pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
        pytplot.xlim( plot_trange[0], plot_trange[1] )

        pytplot.tplot( [wna,ofa_B,integral,'re_in',deltaB,magt,magxy], var_label=labels, xsize=10, ysize=15)

    else:
         """
