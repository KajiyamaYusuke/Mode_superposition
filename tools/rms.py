import numpy as np

def calculate_spl(filename):
    print(f"--- Processing: {filename} ---")
    try:
        # 1. データを読み込む
        # comments='#' でヘッダ行をスキップします
        data = np.loadtxt(filename, comments='#')
        dt = 5e-5
        
        # 2列目(インデックス1)が圧力データと仮定
        # もし1列しかないデータの場合は pressure = data とする
        if data.ndim == 1:
            pressure = data
        else:
            pressure = data[:, 1]
        
        t_start = 0.15
        t_end   = 0.4

        start_idx = int(t_start / dt)
        end_idx   = int(t_end / dt)

        # エラー防止（データが短い場合）
        if end_idx > len(pressure):
            print(f"警告: 指定された終了時間({t_end}s)がデータ長を超えています。末尾まで使用します。")
            end_idx = len(pressure)

        valid_pressure = pressure[start_idx:end_idx]
        
        # 3. DC成分（直流バイアス）の除去
        # RMSは交流成分の強さを見るものなので、中心を0に合わせます
        valid_pressure = valid_pressure - np.mean(valid_pressure)
        
        # 4. RMS (Root Mean Square: 実効値) の計算
        # 式: sqrt( (x1^2 + x2^2 + ... + xn^2) / n )
        rms_pressure = np.sqrt(np.mean(valid_pressure**2))
        
        # 5. SPL (Sound Pressure Level) への変換
        # 基準音圧 P0 = 20 μPa (20e-6 Pa)
        p0 = 20e-6
        
        if rms_pressure > 0:
            spl = 20 * np.log10(rms_pressure / p0)
        else:
            spl = -np.inf # 完全な無音の場合
            
        print(f"RMS Pressure: {rms_pressure:.4f} Pa")
        print(f"SPL         : {spl:.2f} dB")
        
        return spl

    except Exception as e:
        print(f"Error: {e}")
        return None

# ================================
# 実行部分
# ================================
# 計算したいファイル名を指定して実行
spl = calculate_spl("../output/pressure_vt.dat")
