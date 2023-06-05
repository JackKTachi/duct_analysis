import numpy as np
import pyspedas
import pytplot
import math

def calc_stokes(download_trange, ptr):
    from pyspedas.erg import mgf
    mgf(trange=download_trange, datatype='8sec')

    from pyspedas.erg import mgf
    mgf(trange=ptr, datatype='64hz', coord='sgi')

    mag = pytplot.get_data('erg_mgf_l2_magt_8sec')

    # Calculate fce
    fce = mag[1] * 1.0 / (10 ** 9) * 1.6 * (10 ** -19) / (9.1093 * (10 ** -31)) / 2. / math.pi / 1000.  # in kHz

    # Calculate fci
    fci = mag[1] * 1.0 / (10 ** 9) * 1.6 * (10 ** -19) / (1.6726 * (10 ** -27)) / 2. / math.pi / 1000.  # in kHz

    # Calculate flh
    import numpy as np
    flh = np.sqrt(fce * fci)

    # Store data
    pytplot.store_data('fce', data={'x': mag[0], 'y': fce})
    pytplot.store_data('fce2', data={'x': mag[0], 'y': fce * 0.5})
    pytplot.store_data('flh', data={'x': mag[0], 'y': flh})

    pytplot.options('fce', opt_dict={'ytitle': 'fce [kHz]', 'Color': 'red', 'thick': 1})
    pytplot.options('fce2', opt_dict={'ytitle': 'fce [kHz]', 'Color': 'red', 'thick': 1, 'line_style': '--'})
    pytplot.options('flh', opt_dict={'ytitle': 'flh [kHz]', 'Color': 'black', 'thick': 1, 'line_style': '--'})

    from pyspedas.erg import pwe_ofa
    pwe_ofa(trange=download_trange, datatype='spec')

    lst = download_trange[0].split('-')
    pytplot.cdf_to_tplot('./erg_data/satellite/erg/pwe/ofa/l2/matrix/sgi/'+lst[0]+'/'+lst[1]+'/erg_pwe_ofa_l2_matrix_sgi_'+lst[0]+lst[1]+lst[2]+'_v01_03.cdf')

    from pyspedas.erg import pwe_hfa
    pwe_hfa(trange=download_trange, mode='low')
    from pyspedas.erg import pwe_hfa
    pwe_hfa(trange=download_trange, mode='high')
    sE00 = pytplot.data_quants['Ex_Ex_132']
    sE01 = pytplot.data_quants['Ex_Ey_132']
    sE10 = pytplot.data_quants['Ey_Ex_132']
    sE11 = pytplot.data_quants['Ey_Ey_132']
    rrE = np.zeros((2, 2, len(sE00['time']), len(sE00['v1_dim']), 2))

    rrE[0, 0, :, :, :] = sE00.data
    rrE[1, 0, :, :, :] = sE01.data
    rrE[0, 1, :, :, :] = sE10.data
    rrE[1, 1, :, :, :] = sE11.data

    pytplot.split_vec('erg_mgf_l2_mag_64hz_sgi')
    pytplot.store_data('interpotime', data={'x': sE00['time'], 'y': sE00[:,0,0]})

    pyspedas.tinterpol(names='erg_mgf_l2_mag_64hz_sgi_x', interp_to='interpotime', newname='erg_mgf_l2_mag_64hz_sgi_x_interp')
    pyspedas.tinterpol(names='erg_mgf_l2_mag_64hz_sgi_y', interp_to='interpotime', newname='erg_mgf_l2_mag_64hz_sgi_y_interp')
    pyspedas.tinterpol(names='erg_mgf_l2_mag_64hz_sgi_z', interp_to='interpotime', newname='erg_mgf_l2_mag_64hz_sgi_z_interp')

    data_x = pytplot.data_quants['erg_mgf_l2_mag_64hz_sgi_x_interp']
    data_y = pytplot.data_quants['erg_mgf_l2_mag_64hz_sgi_y_interp']
    data_z = pytplot.data_quants['erg_mgf_l2_mag_64hz_sgi_z_interp']

    import numpy as np

    rotmatE = np.zeros((2, 2, len(data_x['time'])))
    rotmat_tE = np.zeros((2, 2, len(data_x['time'])))

    for i in range(len(data_x['time'])):
        bvec = [data_x[i], data_y[i], data_z[i]]
        zz = [0., 0., 1.]
        
        xhat = bvec - np.dot(bvec, zz)
        yhat = np.cross(zz, xhat)
        
        xhat = xhat / np.linalg.norm(xhat)
        yhat = yhat / np.linalg.norm(yhat)
        
        rotmatE[:, :, i] = np.vstack((xhat[:2], yhat[:2]))
        rotmat_tE[:, :, i] = np.transpose(np.vstack((xhat[:2], yhat[:2])))

    for i in range (len(data_x['time'])):
        for j in range (132):
            for k in range (2):
                rrE[:,:,i,j,k] = np.dot(np.dot(rotmatE[:,:,i], rrE[:,:,i,j,k]), rotmat_tE[:,:,i])

    data_xx = pytplot.data_quants['Ex_Ex_132']
    pytplot.store_data('Ex_Ex_re', data={'x': data_xx['time'], 'y':data_xx[:,:,0], 'v': data_xx['v1_dim']})
    pytplot.store_data('Ex_Ex_im', data={'x': data_xx['time'], 'y':data_xx[:,:,1], 'v': data_xx['v1_dim']})

    data_yy = pytplot.data_quants['Ey_Ey_132']
    pytplot.store_data('Ey_Ey_re', data={'x': data_yy['time'], 'y':data_yy[:,:,0], 'v': data_yy['v1_dim']})
    pytplot.store_data('Ey_Ey_im', data={'x': data_yy['time'], 'y':data_yy[:,:,1], 'v': data_yy['v1_dim']})

    pytplot.options('Ex_Ex_re', opt_dict={'ytitle': 'Ex_Ex_re', 'ysubtitle': '[kHz]', 'ztitle':'[pT/Hz]', 'yrange':[0.064, 20], 'zrange':[1e-3,1e2]})
    pytplot.options('Ex_Ex_im', opt_dict={'ytitle': 'Ex_Ex_im', 'ysubtitle': '[kHz]', 'ztitle':'[pT/Hz]', 'yrange':[0.064, 20], 'zrange':[1e-3,1e2]})
    pytplot.options('Ey_Ey_re', opt_dict={'ytitle': 'Ey_Ey_re', 'ysubtitle': '[kHz]', 'ztitle':'[pT/Hz]', 'yrange':[0.064, 20], 'zrange':[1e-3,1e2]})
    pytplot.options('Ey_Ey_im', opt_dict={'ytitle': 'Ey_Ey_im', 'ysubtitle': '[kHz]', 'ztitle':'[pT/Hz]', 'yrange':[0.064, 20], 'zrange':[1e-3,1e2]})

    import numpy as np

    st_I = np.zeros((len(sE00['time']), len(sE00['v1_dim'])))
    st_Q = np.zeros((len(sE00['time']), len(sE00['v1_dim'])))
    st_U = np.zeros((len(sE00['time']), len(sE00['v1_dim'])))
    st_V = np.zeros((len(sE00['time']), len(sE00['v1_dim'])))

    for i in range(len(sE00['time'])):
        for j in range(len(sE00['v1_dim'])):
            st_I[i, j] = rrE[0, 0, i, j, 0] + rrE[1, 1, i, j, 0]
            st_Q[i, j] = rrE[0, 0, i, j, 0] - rrE[1, 1, i, j, 0]
            st_U[i, j] = 2 * rrE[0, 1, i, j, 0]
            st_V[i, j] = 2 * rrE[0, 1, i, j, 1]
    st_chi = np.arctan2(st_U, st_Q) / 2

    pytplot.store_data('st_I', data={'x': sE00['time'], 'y': st_I, 'v': sE00['v1_dim']})
    pytplot.store_data('st_Q/I', data={'x': sE00['time'], 'y': st_Q/st_I, 'v': sE00['v1_dim']})
    pytplot.store_data('st_U/I', data={'x': sE00['time'], 'y': st_U/st_I, 'v': sE00['v1_dim']})
    pytplot.store_data('st_V/I', data={'x': sE00['time'], 'y': st_V/st_I, 'v': sE00['v1_dim']})
    pytplot.store_data('st_chi', data={'x': sE00['time'], 'y': st_chi, 'v': sE00['v1_dim']})

    import add_spec_bins0 as asb

    data_I = asb.add_spec_bins(pytplot.data_quants['st_I'])
    data_Q = asb.add_spec_bins(pytplot.data_quants['st_Q/I'])
    data_U = asb.add_spec_bins(pytplot.data_quants['st_U/I'])
    data_V = asb.add_spec_bins(pytplot.data_quants['st_V/I'])
    data_chi = asb.add_spec_bins(pytplot.data_quants['st_chi'])

    pytplot.store_data('st_I_spec', data={'x': data_I['time'], 'y': data_I, 'v': data_I['v_dim']})
    pytplot.store_data('st_Q/I_spec', data={'x': data_I['time'], 'y': data_Q, 'v': data_I['v_dim']})
    pytplot.store_data('st_U/I_spec', data={'x': data_I['time'], 'y': data_U, 'v': data_I['v_dim']})
    pytplot.store_data('st_V/I_spec', data={'x': data_I['time'], 'y': data_V, 'v': data_I['v_dim']})
    pytplot.store_data('st_chi_spec', data={'x': data_I['time'], 'y': data_chi, 'v': data_I['v_dim']})

    pytplot.options('st_I_spec', opt_dict={'ytitle': 'st_I', 'ysubtitle': '[kHz]', 'ztitle':'', 'ylog':1,'yrange':[0.064, 20], 'spec':1, 'zlog':1,'zrange':[1e-11, 1e0]})
    pytplot.options('st_Q/I_spec', opt_dict={'ytitle': 'st_Q/I', 'ysubtitle': '[kHz]', 'ztitle':'', 'ylog':1,'yrange':[0.064, 20], 'spec':1,'zrange':[-1, 1], 'Colormap': 'RdBu'})
    pytplot.options('st_U/I_spec', opt_dict={'ytitle': 'st_U/I', 'ysubtitle': '[kHz]', 'ztitle':'', 'ylog':1,'yrange':[0.064, 20], 'spec':1,'zrange':[-1,1], 'Colormap': 'RdBu'})
    pytplot.options('st_V/I_spec', opt_dict={'ytitle': 'st_V/I', 'ysubtitle': '[kHz]', 'ztitle':'', 'ylog':1,'yrange':[0.064, 20], 'spec':1,'zrange':[-1,1], 'Colormap': 'RdBu'})
    pytplot.options('st_chi_spec', opt_dict={'ytitle': 'st_chi', 'ysubtitle': '[kHz]', 'ztitle':'', 'ylog':1,'yrange':[0.064, 20], 'spec':1,'zrange':[-1,1], 'Colormap': 'RdBu'})

    pytplot.store_data( 'erg_pwe_ofa_l2_spec_E_spectra_merged', data = ['erg_pwe_ofa_l2_spec_E_spectra_132','fce', 'fce2', 'flh','erg_pwe_ofa_l2_spec_E_spectra_132'])
    pytplot.store_data( 'erg_pwe_ofa_l2_spec_B_spectra_merged', data = ['erg_pwe_ofa_l2_spec_B_spectra_132','fce', 'fce2', 'flh', 'erg_pwe_ofa_l2_spec_B_spectra_132'])
    pytplot.store_data( 'erg_pwe_hfa_l2_spectra_e_mix', data = ['erg_pwe_hfa_l2_high_spectra_e_mix', 'fce', 'fce2', 'flh','erg_pwe_hfa_l2_high_spectra_e_mix'])

    pytplot.options('erg_pwe_ofa_l2_spec_E_spectra_merged', opt_dict={'ytitle': 'E-spectra', 'ysubtitle': '[kHz]', 'ylog':1, 'yrange':[0.064, 20], 'spec':1})
    pytplot.options('erg_pwe_ofa_l2_spec_B_spectra_merged', opt_dict={'ytitle': 'B-spectra', 'ysubtitle': '[kHz]', 'ylog':1, 'yrange':[0.064, 20], 'spec':1})
    pytplot.options('erg_pwe_hfa_l2_spectra_e_mix', opt_dict={'ytitle': 'e-spectra', 'ysubtitle': '[kHz]', 'ylog':1, 'yrange':[20, 100], 'spec':1})

    pytplot.store_data('st_Q/I_spec_log0', data={'x': data_I['time'], 'y': data_Q, 'v': data_I['v_dim']})
    pytplot.options('st_Q/I_spec_log0', opt_dict={'ytitle': 'st_Q/I', 'ysubtitle': '[kHz]', 'ztitle':'', 'ylog':0,'yrange':[0, 20], 'spec':1,'zrange':[-1, 1], 'Colormap': 'RdBu'})

    pytplot.ylim( 'erg_pwe_hfa_l2_high_spectra_e_mix', 20, 100)
    """ pytplot.ylim( 'fce', 0.064, 20.0)
    pytplot.ylim( 'fce2', 0.064, 20.0)
    pytplot.ylim( 'flh', 0.064, 20.0) """
    pytplot.zlim('st_I_spec', 1e-10, 1e0)
    #pytplot.ylim( 'erg_pwe_hfa_l2_low_spectra_e_mix', 0.064, 20.0)
    pytplot.ylim('erg_pwe_hfa_l2_spectra_e_mix', 20, 100)
    pytplot.ylim('erg_pwe_ofa_l2_spec_E_spectra_merged', 0.064, 20.0)
    pytplot.ylim('erg_pwe_ofa_l2_spec_B_spectra_merged', 0.064, 20.0)
    pytplot.xlim(ptr[0], ptr[1])
    #pytplot.tplot( ['erg_pwe_hfa_l2_spectra_e_mix','erg_pwe_ofa_l2_spec_E_spectra_merged', 'erg_pwe_ofa_l2_spec_B_spectra_merged', 'st_I_spec','st_V/I_spec', 'st_Q/I_spec','st_U/I_spec','st_chi_spec'], xsize=10, ysize=25)

    return 'erg_pwe_hfa_l2_spectra_e_mix','erg_pwe_ofa_l2_spec_E_spectra_merged', 'erg_pwe_ofa_l2_spec_B_spectra_merged', 'st_I_spec','st_V/I_spec', 'st_Q/I_spec','st_Q/I_spec_log0','st_U/I_spec','st_chi_spec', 'erg_mgf_l2_mag_64hz_sgi_x_interp', 'erg_mgf_l2_mag_64hz_sgi_y_interp', 'erg_mgf_l2_mag_64hz_sgi_z_interp'














