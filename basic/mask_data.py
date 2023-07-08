import pytplot
import pyspedas
import numpy as np
import xarray as xr

def mask_data(ofa_B, wna_input, pT2perHz):

    B = pytplot.data_quants[ofa_B]

    # 値がpT2perHzを下回る要素をマスク
    mask_data = xr.where(B > pT2perHz, B, np.nan)

    pytplot.store_data('B_mask', data={'x': mask_data['time'], 'y': mask_data, 'v': B['v_dim']})
    B_mask = pytplot.data_quants['B_mask']        
    """ 
        B = pytplot.data_quants[ofa_B]
        mask_data = B.where(B > pT2perHz, drop=True)
        pytplot.store_data('B_mask', data={'x': mask_data['time'], 'y': mask_data, 'v': B['v_dim']})
        B_mask = pytplot.data_quants['B_mask']
    """
    new_coords = {'v_dim': [0.064,  0.128,  0.192,  0.256,  0.32 ,  0.384,  0.448,  0.512,
            0.576,  0.64 ,  0.704,  0.768,  0.832,  0.896,  0.96 ,  1.024,
            1.088,  1.152,  1.216,  1.28 ,  1.344,  1.408,  1.472,  1.536,
            1.6  ,  1.664,  1.728,  1.792,  1.856,  1.92 ,  1.984,  2.048,
            2.112,  2.176,  2.24 ,  2.304,  2.368,  2.432,  2.496,  2.56 ,
            2.624,  2.688,  2.752,  2.816,  2.88 ,  2.944,  3.008,  3.072,
            3.136,  3.2  ,  3.264,  3.328,  3.392,  3.456,  3.52 ,  3.584,
            3.648,  3.712,  3.776,  3.84 ,  3.904,  3.968,  4.032,  4.096,
            4.16 ,  4.224,  4.288,  4.352,  4.416,  4.48 ,  4.544,  4.608,
            4.672,  4.736,  4.8  ,  4.864,  4.928,  4.992,  5.056,  5.12 ,
            5.184,  5.248,  5.312,  5.376,  5.44 ,  5.504,  5.568,  5.632,
            5.696,  5.76 ,  5.824,  5.888,  5.952,  6.016,  6.08 ,  6.144,
            6.208,  6.272,  6.336,  6.4  ,  6.464,  6.528,  6.592,  6.656,
            6.72 ,  6.784,  6.848,  6.912,  7.104,  7.424,  7.776,  8.128,
            8.48 ,  8.864,  9.248,  9.664, 10.112, 10.56 , 11.04 , 11.552,
        12.064, 12.576, 13.12 , 13.728, 14.368, 15.008, 15.648, 16.352,
        17.088, 17.824, 18.624, 19.456]}
    data = B_mask.assign_coords(new_coords)
    pytplot.store_data('B_mask_add', data={'x': data['time'], 'y': data, 'v': data['v_dim']})

    pytplot.options('B_mask_add', opt_dict={'ylog': True, 'zlog': True, 'zrange': [1e-4, 100], 'spec': True,  'yrange': [0.064, 20], 'ztitle': '[nT^2/Hz]', 'ytitle': 'OFA_B', 'ysubtitle': '[kHz]'})
    import pyspedas
    pyspedas.tinterpol(wna_input, 'B_mask_add', newname='wna_intrp')
    wna = pytplot.data_quants['wna_intrp']

    indexes = xr.DataArray(B>pT2perHz, dims=['time', 'v_dim'])
    wna_mask = wna.where(indexes, drop=None)
    pytplot.store_data('wna_mask', data={'x': wna_mask['time'], 'y': wna_mask, 'v': wna_mask['v_dim']})
    wna_masked = pytplot.data_quants['wna_mask']
    data_wna = wna_masked.assign_coords(new_coords)
    pytplot.store_data('wna_mask_add', data={'x': data_wna['time'], 'y': data_wna, 'v': data_wna['v_dim']})
    pytplot.options('wna_mask_add', opt_dict={'ylog': True, 'zrange': [0,20], 'spec': True,  'yrange': [0.064, 20], 'ztitle': '[degree]', 'ytitle': 'WNA', 'ysubtitle': '[kHz]'})

    return 'B_mask_add', 'wna_mask_add'



