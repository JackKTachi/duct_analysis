import numpy as np
import pyspedas
import pytplot

def integrate_chorus(download_trange,r1, min, max):

    from pyspedas.erg import pwe_ofa
    pwe_ofa(trange=download_trange)

    data1 = pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132']['spec_bins']

    index_min = -1
    for i, value in enumerate(data1):
        if value > min:
            index_min = i
            break

    index_max = -1
    for i, value in enumerate(data1):
        if value < max:
            index_max = i
        else:
            break


    import datetime
    # 元の時間配列
    time_array = pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132']['time']

    # 指定された範囲
    r1_start = datetime.datetime.strptime(r1[0], '%Y-%m-%d %H:%M:%S')
    r1_end = datetime.datetime.strptime(r1[1], '%Y-%m-%d %H:%M:%S')
    # 始まりのインデックスを取得
    start_index = -1
    for i, timestamp in enumerate(time_array):
        if timestamp >= np.datetime64(r1_start):
            start_index = i
            break

    # 終わりのインデックスを取得
    end_index = -1
    for i in range(len(time_array) - 1, -1, -1):
        timestamp = time_array[i]
        if timestamp <= np.datetime64(r1_end):
            end_index = i
            break

    data = np.zeros(len(pytplot.get_data('erg_pwe_ofa_l2_spec_B_spectra_132')[0]))

    for i in range(start_index, end_index+1):
        data[i] = np.sum(pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132'][i][index_min:index_max+1])

    pytplot.store_data('erg_pwe_ofa_l2_spec_B_chorus_integrate', data={'x': pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132']['time'], 'y': data})
    return 'erg_pwe_ofa_l2_spec_B_chorus_integrate'

def integral_chorus_spec(download_trange,r1, min, max):

    from pyspedas.erg import pwe_ofa
    pwe_ofa(trange=download_trange)

    import pytplot
    
    data1 = pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132']['spec_bins']

    index_min = -1
    for i, value in enumerate(data1):
        if value > min:
            index_min = i
            break

    index_max = -1
    for i, value in enumerate(data1):
        if value < max:
            index_max = i
        else:
            break


    import datetime
    # 元の時間配列
    time_array = pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132']['time']

    # 指定された範囲
    r1_start = datetime.datetime.strptime(r1[0], '%Y-%m-%d %H:%M:%S')
    r1_end = datetime.datetime.strptime(r1[1], '%Y-%m-%d %H:%M:%S')
    # 始まりのインデックスを取得
    start_index = -1
    import numpy as np
    for i, timestamp in enumerate(time_array):
        if timestamp >= np.datetime64(r1_start):
            start_index = i
            break

    # 終わりのインデックスを取得
    end_index = -1
    for i in range(len(time_array) - 1, -1, -1):
        timestamp = time_array[i]
        if timestamp <= np.datetime64(r1_end):
            end_index = i
            break

    data = np.zeros(len(pytplot.get_data('erg_pwe_ofa_l2_spec_B_spectra_132')[0]))

    for i in range(start_index, end_index+1):
        for j in range(index_min, index_max+1):
            data[i] += pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132'][i][j]*pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132']['spec_bins'][j]*1e3
    import numpy as np
    import pytplot
    average_spin1 = 8
    data_int1 = np.zeros(data.size)
    for i in range(data.size):
        idx1 = np.arange(i-average_spin1//2, i+average_spin1//2+1).astype(int)
        idx1 = np.clip(idx1, 0, data.size-1)
        data_int1[i] = np.sum(data[idx1]) / average_spin1

    pytplot.store_data('erg_pwe_ofa_l2_spec_B_chorus_integrate', data={'x': pytplot.data_quants['erg_pwe_ofa_l2_spec_B_spectra_132']['time'], 'y': data_int1})



