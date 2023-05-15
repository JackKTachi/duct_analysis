import pyspedas
import pytplot 
import numpy as np
from pytplot.tplot import tplot

def compressional_B(tr):
    from pyspedas.erg import mgf
    mgf(trange = tr)
    pytplot.split_vec('erg_mgf_l2_mag_8sec_dsi')
    # 前後25ポイントのデータを平均
    data_x = pytplot.data_quants['erg_mgf_l2_mag_8sec_dsi_x']
    data_y = pytplot.data_quants['erg_mgf_l2_mag_8sec_dsi_y']
    data_z = pytplot.data_quants['erg_mgf_l2_mag_8sec_dsi_z']
    average_spin = 25

    data_200s_x = np.zeros(data_x['time'].size)
    data_200s_y = np.zeros(data_x['time'].size)
    data_200s_z = np.zeros(data_x['time'].size)

    for i in range(data_x['time'].size):
        idx = np.arange(i-average_spin//2, i+average_spin//2+1).astype(int)
        idx = np.clip(idx, 0, data_x['time'].size-1)

        data_200s_x[i] = np.sum(data_x[idx]) / average_spin
        data_200s_y[i] = np.sum(data_y[idx]) / average_spin
        data_200s_z[i] = np.sum(data_z[idx]) / average_spin

    rotmat=np.zeros((3, 3, data_200s_x.size))

    for i in range(data_200s_x.size):
        bvec = [data_200s_x[i], data_200s_y[i], data_200s_z[i]]
        zz = [0, 0, 1]

        yhat = np.cross(zz, bvec) # in dsi coordinate, z-axis is roughly anti-sunward
        xhat = np.cross(yhat, bvec) # right-handed orthogonal coordinate system
        zhat = bvec

        # 単位ベクトルに変換
        yhat = yhat / np.linalg.norm(yhat)
        xhat = xhat / np.linalg.norm(xhat)
        zhat = zhat / np.linalg.norm(zhat)

        # 回転行列を作成
        rotmat[:,:,i] = np.array([xhat, yhat, zhat])

        data_rot_200s = np.zeros((data_200s_x.size, 3))

    for i in range(data_200s_x.size):
        data_200s = [data_200s_x[i], data_200s_y[i], data_200s_z[i]]
        data_rot_200s[i,:] = np.dot(rotmat[:,:,i], data_200s)

    data_8sec_x = pytplot.data_quants['erg_mgf_l2_mag_8sec_dsi_x']
    data_8sec_y = pytplot.data_quants['erg_mgf_l2_mag_8sec_dsi_y']
    data_8sec_z = pytplot.data_quants['erg_mgf_l2_mag_8sec_dsi_z']

    # rotate above dates
    data_rot_8sec = np.zeros((data_8sec_x['time'].size, 3))

    for i in range(data_8sec_x['time'].size):
        data = [data_8sec_x[i], data_8sec_y[i], data_8sec_z[i]]
        data_rot_8sec[i,:] = np.dot(rotmat[:,:,i], data)
    
    delta_z = np.zeros(data_rot_8sec[:,0].size)

    for i in range(data_rot_8sec[:,0].size):
        delta_z[i] = data_rot_200s[i,2] - data_rot_8sec[i,2]

    data_rot_8sec_xy = np.zeros((data_rot_8sec[:,0].size, 2))
    for i in range(data_rot_8sec[:,0].size):
        data_rot_8sec_xy[i,0] = data_rot_8sec[i,0]
        data_rot_8sec_xy[i,1] = data_rot_8sec[i,1]

    pytplot.store_data('delta_z', data={'x': data_8sec_x['time'], 'y': delta_z})
    pytplot.store_data('erg_mgf_l2_mag_8sec_MAF_x&y', data={'x': data_8sec_x['time'], 'y': data_rot_8sec_xy})
    pytplot.options('erg_mgf_l2_mag_8sec_MAF_x&y', opt_dict={'legend_names': ['$B_x$', '$B_y$'], 'ytitle': '$B_x & B_y$', 'ysubtitle': '[nT]'})
    pytplot.options('delta_z', opt_dict={'legend_names': [r'$\delta B_z$'], 'ytitle': r'$\delta B_z$', 'ysubtitle': '[nT]'})
    pytplot.options('erg_mgf_l2_magt_8sec', opt_dict={'ytitle': r'$Btotal$', 'ysubtitle': '[nT]'})

    return 'delta_z', 'erg_mgf_l2_mag_8sec_MAF_x&y', 'erg_mgf_l2_magt_8sec'


    






