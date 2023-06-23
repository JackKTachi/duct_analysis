import pytplot
import pyspedas
import plasma_params as pp
import numpy as np

def calc_refractive_index(download_trange, f, theta, density=None):
    from pyspedas.erg import mgf
    mgf(trange=download_trange)

    if density is not None:
        B = pytplot.data_quants['erg_mgf_l2_magt_8sec']
        N = density
        fc = (pp.Q*B*f/pp.ME)**(1/2)
        fp = fp = (80.6*N)**(1/2)
        re_in = fp*np.cos(theta)/(fc*f*(np.cos(theta)-f/fc))**0.5

    else:
        
        lst = download_trange[0].split('-')
        pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/hfa/l3/1min/'+lst[0]+'/'+lst[1]+'/erg_pwe_hfa_l3_1min_'+lst[0]+lst[1]+lst[2]+'_v03_07.cdf')
        pyspedas.tinterpol('ne_mgf', 'erg_mgf_l2_magt_8sec', newname='ne_mgf_interpolated')
        B = pytplot.data_quants['erg_mgf_l2_magt_8sec']
        N = pytplot.data_quants['ne_mgf_interpolated']
        fc = (pp.Q*B*f/pp.ME)**(1/2)
        fp = fp = (80.6*N)**(1/2)
        re_in = fp*np.cos(theta)/(fc*f*(np.cos(theta)-f/fc))**0.5    


    pytplot.store_data('re_in', data={'x': pytplot.data_quants['ne_mgf_interpolated']['time'], 'y': re_in})
    pytplot.options('re_in', opt_dict={'ytitle':'Refractive Index','ylog':0})

    return 're_in'



