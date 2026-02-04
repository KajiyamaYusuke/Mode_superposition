import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# =========================
# 設定
# =========================
filename = "../output/pressure_vt.dat"  # ファイルパス
sim_dt = 1.0e-5
output_interval = 5
dt = sim_dt * output_interval

# =========================
# データ処理
# =========================
try:
    data = np.loadtxt(filename, comments='#')
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
    valid_pressure = valid_pressure - np.mean(valid_pressure)
    window = np.hanning(len(valid_pressure))
    
    N = len(valid_pressure)
    freq = np.fft.rfftfreq(N, d=dt)
    fft_val = np.fft.rfft(valid_pressure * window)
    
    # dB変換
    p0 = 20e-6
    db_amplitude = 20 * np.log10(np.abs(fft_val) / N * 2 + 1e-12 / p0)
    
    # =========================
    # ★重要：ロバストなF0検出
    # =========================
    # 1. 低周波ノイズ(50Hz以下)を除外した範囲を作る
    min_freq = 50.0
    mask = freq > min_freq
    masked_freq = freq[mask]
    masked_db   = db_amplitude[mask]
    
    # 2. ピーク（山）をすべて探す
    #    条件: 最大値から 20dB 以内にあるピークのみを候補とする
    max_db = np.max(masked_db)
    threshold = max_db - 20.0 
    
    peaks, _ = find_peaks(masked_db, height=threshold, distance=5)
    
    f0 = 0.0
    if len(peaks) > 0:
        # 3. 見つかった候補の中で「一番周波数が低いもの」を選ぶ
        # peaks はインデックスのリストなので、[0]番目が最低周波数
        first_peak_idx = peaks[0]
        f0 = masked_freq[first_peak_idx]
        f0_db = masked_db[first_peak_idx]
        
        # 確認用: 最大ピークがF0でない場合の警告
        max_idx = np.argmax(masked_db)
        if first_peak_idx != max_idx:
            print(f"【注意】最大ピーク({masked_freq[max_idx]:.1f}Hz)はF0ではありません。")
            print(f"倍音の方が強いですが、正しくF0({f0:.1f}Hz)を検出しました。")
        else:
            print(f"F0検出: {f0:.1f} Hz (最大ピークと一致)")
            
    else:
        print("有効なピークが見つかりませんでした。")

    # =========================
    # プロット
    # =========================
    plt.figure(figsize=(10, 5))
    plt.plot(freq, db_amplitude, label='Spectrum')
    plt.xlim(0, 3000)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("SPL [dB]")
    plt.grid(True)
    plt.title(f"Fundamental Frequency: {f0:.1f} Hz")
    
    # F0の位置を赤丸で表示
    if f0 > 0:
        plt.plot(f0, f0_db, 'ro', markersize=10, label='Detected F0')
        # その他の候補（倍音など）を青バツで表示
        plt.plot(masked_freq[peaks], masked_db[peaks], 'bx', label='Harmonics')
        
    plt.legend()
    plt.show()

except Exception as e:
    print(f"Error: {e}")