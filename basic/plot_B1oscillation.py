import calc_RI as cr
import datetime
import ULF_duct_plot_detail as udp
import calc_integral_chorus as cic
import ULF_duct_plot_ver2 as udp2
import calc_stokes as cs
import calc_refractive_index as cri
import pytplot
import numpy as np
import plasma_params as pp
import compressional_B as cb

def plotB1(tr, f):
    download_trange = []

    for item in tr:
        dt = datetime.datetime.strptime(item, '%Y-%m-%d %H:%M:%S')
        date_only = dt.strftime('%Y-%m-%d')
        download_trange.append(date_only)

    lst = download_trange[0].split('-')
    #pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/hfa/l3/1min/'+lst[0]+'/'+lst[1]+'/erg_pwe_hfa_l3_1min_'+lst[0]+lst[1]+lst[2]+'_v03_07.cdf')    
        
    from pyspedas.erg import pwe_ofa
    pwe_ofa(trange=download_trange)

    #pytplot.cdf_to_tplot('./erg_data/satellite/erg/mepe/l3/pa/'+lst[0]+'/'+lst[1]+'/erg_mepe_l3_pa_'+lst[0]+lst[1]+lst[2]+'_v01_01.cdf')
    #data1=pytplot.data_quants['FEDU']
    #'erg_mepe_l3_pa_FEDU_87.5keV', 'erg_mepe_l3_pa_FEDU_72.6keV', 'erg_mepe_l3_pa_FEDU_60.4keV', 'erg_mepe_l3_pa_FEDU_50.3keV', 'erg_mepe_l3_pa_FEDU_42.0keV','erg_mepe_l3_pa_FEDU_35.0keV','erg_mepe_l3_pa_FEDU_29.3keV','erg_mepe_l3_pa_FEDU_24.5keV','erg_mepe_l3_pa_FEDU_20.5keV','erg_mepe_l3_pa_FEDU_17.1keV','erg_mepe_l3_pa_FEDU_14.3keV','erg_mepe_l3_pa_FEDU_12.0keV','erg_mepe_l3_pa_FEDU_10.0keV','erg_mepe_l3_pa_FEDU_8.4keV','erg_mepe_l3_pa_FEDU_7.0keV' == eep.each_eV_plot(data1)
    
    pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/ofa/l3/property/'+lst[0]+'/'+lst[1]+'/erg_pwe_ofa_l3_property_dsi_'+lst[0]+lst[1]+lst[2]+'_v01_03.cdf')
    #pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/ofa/l3/property/2017/03/erg_pwe_ofa_l3_property_dsi_20170329_v01_03.cdf')
    # 2017Mar30 is the file of v01_04, others are v01_03

    'B_perp_para', 'B1_total', 'delta_z', 'delta_xy', 'ULF_wna', 'erg_mgf_l2_magt_8sec', 'cos_th_intpl' == cb.ULFwna(download_trange)
    #'delta_z', 'erg_mgf_l2_mag_8sec_MAF_x&y', 'erg_mgf_l2_magt_8sec' == cb.ULFwna(download_trange)

    #from pyspedas.erg import mepe
    #mepe(trange=download_trange)

    from pyspedas.erg import orb
    orb(trange=download_trange)

    B = pytplot.data_quants['erg_mgf_l2_magt_8sec']*1e-9
    N0 = 2e6
    mu = cr.calc_RI(N0, B, f)
    pytplot.store_data('mu', data={'x':B['time'], 'y':mu})

    min = 0.1
    max = 20
    'erg_pwe_ofa_l2_spec_B_chorus_integrate' == cic.integral_chorus_spec(download_trange,tr, min, max)

    pytplot.options('kvec_polar_132', opt_dict={'ytitle':'wna','ysubtitle':'[kHz]','ylog':1, 'spec':1})
    pytplot.options('erg_pwe_ofa_l2_spec_B_spectra_132', opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options('erg_pwe_ofa_l2_spec_B_chorus_integrate', opt_dict={'ytitle':'chorus intensity', 'ysubtitle':'[$pT^2$]','ylog':1})
    pytplot.options('delta_z', opt_dict={'ytitle':'$\delta B$','ylog':0})
    pytplot.options('erg_mgf_l2_magt_8sec', opt_dict={'ytitle':'$B_{total}$','ylog':0})
    pytplot.options('mu', opt_dict={'ytitle':'$\mu$','ylog':1})
    labels = pytplot.split_vec( 'erg_orb_l2_pos_rmlatmlt' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )

    pytplot.xlim( tr[0], tr[1] )

    pytplot.tplot(['kvec_polar_132', 'erg_pwe_ofa_l2_spec_B_spectra_132', 'erg_pwe_ofa_l2_spec_B_chorus_integrate', 'mu', 'delta_z', 'erg_mgf_l2_magt_8sec', 'delta_xy', 'B1_total'],var_label=labels, xsize=15, ysize=20)














