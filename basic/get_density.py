import numpy as np
import pytplot


def get_density(file_path):

    # ファイルを開く
    with open(file_path, "r") as file:
        # ヘッダー行をスキップする
        next(file)

        # DateとNeの列を格納するリスト
        dates = []
        times = []
        ne_values = []

        # ファイルの各行を読み込む
        for line in file:
            # 行をスペースで分割
            columns = line.split()

            # Date列とNe列の値を取得してリストに追加
            date = columns[0]
            time = columns[1]
            ne = columns[-1]
            dates.append(date)
            times.append(time)
            ne_values.append(ne)

    dates.pop(0)
    times.pop(0)
    ne_values.pop(0)

    hrs = [float(value) for value in times]
    ncc = [float(value) for value in ne_values]

    ts = []
    n = len(dates)
    for i in range(n):
        yr = dates[i][:4]
        mn = dates[i][4:6]
        dy = dates[i][6:8]
        hr = int(hrs[i])
        mm = int(hrs[i] * 60) - hr * 60
        sc = int(hrs[i] * 3600) - hr * 3600 - mm * 60
        ts.append(f"{yr}-{mn}-{dy}/{hr:02}:{mm:02}:{sc:02}")

    from datetime import datetime
    tds = []  # 変換後の数値の時刻データを格納するリスト
    dts = []
    # 文字列をdatetimeオブジェクトに変換
    for i in range(n):
        dt = datetime.strptime(ts[i], "%Y-%m-%d/%H:%M:%S")
        dts.append(dt)

    import pytplot
    pytplot.store_data('ne', data={'x':dts, 'y':ncc})

    return 'ne'


