import pyspedas
import pytplot 
import numpy as np
from pytplot.tplot import tplot

def compressional_B(tr):

    from pyspedas.erg import mgf
    mgf(trange = tr)

    from pyspedas.erg import orb
    orb(trange = tr)

    import numpy as np
    import pytplot

    pytplot.split_vec('erg_mgf_l2_mag_8sec_gsm')

    data_x = pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_x']
    data_y = pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_y']
    data_z = pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_z']
    average_spin = 25

    data_200s_total = np.zeros(data_x['time'].size)
    data_200s_x = np.zeros(data_x['time'].size)
    data_200s_y = np.zeros(data_x['time'].size)
    data_200s_z = np.zeros(data_x['time'].size)

    for i in range(data_x['time'].size):
        idx = np.arange(i-average_spin//2, i+average_spin//2+1).astype(int)
        idx = np.clip(idx, 0, data_x['time'].size-1)

        data_200s_x[i] = np.sum(data_x[idx]) / average_spin
        data_200s_y[i] = np.sum(data_y[idx]) / average_spin
        data_200s_z[i] = np.sum(data_z[idx]) / average_spin
        data_200s_total[i] = np.sqrt(data_200s_x[i]**2 + data_200s_y[i]**2 + data_200s_z[i]**2)

    pytplot.split_vec('erg_orb_l2_pos_rmlatmlt')

    cos_th = data_200s_z/ data_200s_total
    sin_th = np.sqrt(data_200s_x**2 + data_200s_y**2) / data_200s_total

    cos_ph = data_200s_x/ np.sqrt(data_200s_x**2 + data_200s_y**2)
    sin_ph = data_200s_y / np.sqrt(data_200s_x**2 + data_200s_y**2)

    beta = pytplot.data_quants['erg_orb_l2_pos_rmlatmlt_y'] * np.pi / 180.0
    alph = np.pi * pytplot.data_quants['erg_orb_l2_pos_rmlatmlt_z'] / 12.0

    cos_alph = np.cos(alph)
    sin_alph = np.sin(alph)
    cos_beta = np.cos(beta)
    sin_beta = np.sin(beta)

    pytplot.store_data('cos_th', data={'x': data_x['time'], 'y': cos_th})
    pytplot.store_data('sin_th', data={'x': data_x['time'], 'y': sin_th})
    pytplot.store_data('cos_ph', data={'x': data_x['time'], 'y': cos_ph})
    pytplot.store_data('sin_ph', data={'x': data_x['time'], 'y': sin_ph})
    pytplot.store_data('cos_alph', data={'x': pytplot.data_quants['erg_orb_l2_pos_rmlatmlt_x']['time'], 'y': cos_alph})
    pytplot.store_data('sin_alph', data={'x': pytplot.data_quants['erg_orb_l2_pos_rmlatmlt_x']['time'], 'y': sin_alph})
    pytplot.store_data('cos_beta', data={'x': pytplot.data_quants['erg_orb_l2_pos_rmlatmlt_x']['time'], 'y': cos_beta})
    pytplot.store_data('sin_beta', data={'x': pytplot.data_quants['erg_orb_l2_pos_rmlatmlt_x']['time'], 'y': sin_beta})

    pyspedas.tinterpol('cos_th', 'cos_alph', newname='cos_th_intpl')
    pyspedas.tinterpol('sin_th', 'cos_alph', newname='sin_th_intpl')
    pyspedas.tinterpol('cos_ph', 'cos_alph', newname='cos_ph_intpl')
    pyspedas.tinterpol('sin_ph', 'cos_alph', newname='sin_ph_intpl')

    gamma = np.arctan((pytplot.data_quants['cos_beta']*(pytplot.data_quants['cos_alph']*pytplot.data_quants['sin_ph_intpl']-pytplot.data_quants['sin_alph']*pytplot.data_quants['cos_ph_intpl']))/(pytplot.data_quants['cos_th_intpl']*pytplot.data_quants['cos_beta']*(pytplot.data_quants['cos_alph']*pytplot.data_quants['cos_ph_intpl']+pytplot.data_quants['sin_alph']*pytplot.data_quants['sin_ph_intpl'])+pytplot.data_quants['sin_th_intpl']*pytplot.data_quants['sin_beta']))

    pyspedas.tinterpol('erg_mgf_l2_mag_8sec_gsm_x', 'cos_alph', newname='erg_mgf_l2_mag_8sec_gsm_x_intpl')
    pyspedas.tinterpol('erg_mgf_l2_mag_8sec_gsm_y', 'cos_alph', newname='erg_mgf_l2_mag_8sec_gsm_y_intpl')
    pyspedas.tinterpol('erg_mgf_l2_mag_8sec_gsm_z', 'cos_alph', newname='erg_mgf_l2_mag_8sec_gsm_z_intpl')

    B_MFA_z = pytplot.data_quants['sin_th_intpl'] * pytplot.data_quants['cos_ph_intpl'] * pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_x_intpl'] + pytplot.data_quants['sin_th_intpl'] * pytplot.data_quants['sin_ph_intpl'] * pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_y_intpl'] + pytplot.data_quants['cos_th_intpl'] * pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_z_intpl']

    delta_B = np.zeros(B_MFA_z.size)

    for i in range(B_MFA_z.size):
        idx1 = np.arange(i-average_spin//2, i+average_spin//2+1).astype(int)
        idx1 = np.clip(idx1, 0, B_MFA_z.size-1)
        delta_B[i] =  B_MFA_z[i] - np.sum(B_MFA_z[idx1]) / average_spin

    B_MFA_x = pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_x_intpl'] * (np.cos(gamma) * pytplot.data_quants['cos_th_intpl'] * pytplot.data_quants['cos_ph_intpl'] + np.sin(gamma) * pytplot.data_quants['sin_ph_intpl']) + pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_y_intpl'] * (np.cos(gamma) * pytplot.data_quants['cos_th_intpl'] * pytplot.data_quants['sin_ph_intpl'] - np.sin(gamma) * pytplot.data_quants['cos_ph_intpl']) - pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_z_intpl'] * np.cos(gamma) * pytplot.data_quants['sin_th_intpl']
    B_MFA_y = pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_x_intpl'] * (np.sin(gamma) * pytplot.data_quants['cos_th_intpl'] * pytplot.data_quants['cos_ph_intpl'] - np.cos(gamma) * pytplot.data_quants['sin_ph_intpl']) + pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_y_intpl'] * (np.sin(gamma) * pytplot.data_quants['cos_th_intpl'] * pytplot.data_quants['sin_ph_intpl'] + np.cos(gamma) * pytplot.data_quants['cos_ph_intpl']) - pytplot.data_quants['erg_mgf_l2_mag_8sec_gsm_z_intpl'] * np.sin(gamma) * pytplot.data_quants['sin_th_intpl']

    B_MFA_xy = np.zeros((B_MFA_x['time'].size, 2))
    for i in range(B_MFA_x['time'].size):
        B_MFA_xy[i,0] = B_MFA_x[i]
        B_MFA_xy[i,1] = B_MFA_y[i]

    pytplot.store_data('delta_z', data={'x': pytplot.data_quants['sin_th_intpl']['time'], 'y': delta_B})
    pytplot.store_data('erg_mgf_l2_mag_8sec_MAF_x&y', data={'x': pytplot.data_quants['sin_th_intpl']['time'], 'y': B_MFA_xy})

    pytplot.options('erg_mgf_l2_mag_8sec_MAF_x&y', opt_dict={'legend_names': ['$B_x$', '$B_y$'], 'ytitle': '$B_x & B_y$', 'ysubtitle': '[nT]'})
    pytplot.options('delta_z', opt_dict={'legend_names': [r'$\delta B_z$'], 'ytitle': r'$\delta B_z$', 'ysubtitle': '[nT]'})
    pytplot.options('erg_mgf_l2_magt_8sec', opt_dict={'ytitle': r'$Btotal$', 'ysubtitle': '[nT]'})

    return 'delta_z', 'erg_mgf_l2_mag_8sec_MAF_x&y', 'erg_mgf_l2_magt_8sec'
