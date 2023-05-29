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

    'delta_z', 'erg_mgf_l2_mag_8sec_MAF_x&y', 'erg_mgf_l2_magt_8sec' == cb.compressional_B(download_trange)

    from pyspedas.erg import mepe
    mepe(trange=download_trange)

    from pyspedas.erg import orb
    orb(trange=download_trange)
    
    return 'ne_mgf',\
        'kvec_polar_132',\
        'erg_pwe_ofa_l2_spec_B_spectra_132','erg_pwe_ofa_l2_spec_E_spectra_132',\
        'erg_mepe_l2_omniflux_FEDO',\
        'erg_mepe_l3_pa_FEDU_87.5keV', 'erg_mepe_l3_pa_FEDU_72.6keV', 'erg_mepe_l3_pa_FEDU_60.4keV','erg_mepe_l3_pa_FEDU_50.3keV', 'erg_mepe_l3_pa_FEDU_42.0keV','erg_mepe_l3_pa_FEDU_35.0keV','erg_mepe_l3_pa_FEDU_29.3keV','erg_mepe_l3_pa_FEDU_24.5keV','erg_mepe_l3_pa_FEDU_20.5keV','erg_mepe_l3_pa_FEDU_17.1keV','erg_mepe_l3_pa_FEDU_14.3keV','erg_mepe_l3_pa_FEDU_12.0keV',\
            'delta_z', 'erg_mgf_l2_magt_8sec','erg_mgf_l2_mag_8sec_MAF_x&y',\
            'erg_orb_l2_pos_rmlatmlt','erg_orb_l2_pos_eq'
            
def ULF_duct_plot(ofa_B,ofa_E,mepe_50,delta_z,magt,magxy,orb_rmlatmlt, plot_trange):
    
    pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options(ofa_E, opt_dict={'ytitle':'ofa-E','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})

    labels = pytplot.split_vec( orb_rmlatmlt )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
    pytplot.xlim( plot_trange[0], plot_trange[1] )

    pytplot.tplot( [ofa_B,ofa_E,mepe_50,delta_z,magt,magxy], var_label=labels, xsize=10, ysize=10)


def ULF_duct_plot_12to50(density, wna,ofa_B,ofa_E,omniflux,mepe87, mepe72, mepe60,mepe_50, mepe_42, mepe_35, mepe_24, mepe_12,delta_z,magt,magxy,orb_rmlatmlt, plot_trange, download_trange=None):
    
    pytplot.options(wna, opt_dict={'ytitle':'wna','ysubtitle':'[kHz]','ylog':1, 'spec':1})

    pytplot.options(ofa_B, opt_dict={'ytitle':'ofa-B','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})
    pytplot.options(ofa_E, opt_dict={'ytitle':'ofa-E','ysubtitle':'[kHz]','ylog':1, 'zlog':1, 'spec':1})

    pytplot.options(omniflux, opt_dict={'ylog': 0, 'ytitle': 'omniflux', 'ysubtitle': '[keV]'})

    labels = pytplot.split_vec( orb_rmlatmlt )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_x', opt_dict={'ytitle':'L [Re]','ytitle_location':'left'} )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_y', 'ytitle', 'MLat [deg]' )
    pytplot.options( 'erg_orb_l2_pos_rmlatmlt_z', 'ytitle', 'MLT [h]' )
    pytplot.xlim( plot_trange[0], plot_trange[1] )

    if download_trange is not None:
        from pyspedas.erg import mgf
        mgf(trange=download_trange)

        import plasma_params as pp
        import numpy as np
        data = pytplot.data_quants['erg_mgf_l2_magt_8sec']
        fce1 = np.zeros(data['time'].size)
        fce2 = np.zeros(data['time'].size)
        for i in range(data['time'].size):
            fce2[i] = pp.Q*data[i]/pp.ME/1e9/2/np.pi/1e3*0.2
            fce1[i] = pp.Q*data[i]/pp.ME/1e9/2/np.pi/1e3*0.1

        pytplot.store_data('fce2', data={'x': data['time'], 'y': fce2})
        pytplot.store_data('fce1', data={'x': data['time'], 'y': fce1})
        pytplot.options('fce2', opt_dict={'ytitle':'ofa-B','ytitle': 'fce2 [kHz]', 'yrange':[1e-1, 3e-0],'ylog': 1, 'line_style': '-', 'Color': 'darkblue', 'thick': 1, 'alpha':1})
        pytplot.options('fce1', opt_dict={'ytitle': 'fce1 [kHz]', 'ylog': 1, 'yrange':[1e-1, 3e-0],'line_style': '--', 'Color': 'black', 'thick': 1, 'alpha':1})

        pytplot.store_data('ofa_B_fce', data=[ofa_B,'fce1','fce2'])
        #pytplot.store_data('ofa_E_fce', data=[ofa_E,'fce1','fce2'])
        #pytplot.store_data('kvec', data= [wna, 'fce1','fce2'])
        pytplot.options('ofa_B_fce', opt_dict={'yrange':[1e-1, 3e-0]})

        from pyspedas.erg import pwe_hfa
        pwe_hfa(trange=download_trange)
        pytplot.options('erg_pwe_hfa_l2_low_spectra_esum', opt_dict={'ytitle':'HFA-esum',})

        pytplot.tplot( [wna,'ofa_B_fce',ofa_E,mepe87, mepe72, mepe60,mepe_50, mepe_42, mepe_35, mepe_24, mepe_12,delta_z,magt,magxy,'erg_pwe_hfa_l2_low_spectra_esum'], var_label=labels, xsize=10, ysize=25)        

    else:
        pytplot.tplot( [density, wna,ofa_B,ofa_E,omniflux,mepe87, mepe72, mepe60,mepe_50, mepe_42, mepe_35, mepe_24, mepe_12,delta_z,magt,magxy], var_label=labels, xsize=10, ysize=30)
