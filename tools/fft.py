import numpy as np
import matplotlib.pyplot as plt

# =========================
# 設定
# =========================
filename = "../output/pressure_vt.dat"   # ファイル名を合わせました
sim_dt = 1.0e-5             # シミュレーションdt
output_interval = 5       # 出力間隔
dt = sim_dt * output_interval  # サンプリング周期

# =========================
# データ読み込みと解析
# =========================
# データ読み込み
data = np.loadtxt(filename, comments='#')

# 1列目はステップ数(時間)、2列目以降が圧力データと仮定
steps = data[:, 0]          # 横軸（ステップ）
pressure = data[:, 1] # 縦軸計算用（2列目以降すべて）

# 各行の最小値を計算（これを波形とする）
#pressure = np.min(pressure_cols, axis=1)

# ★重要：データの後半（振動している部分）だけを使う
start_idx = int(len(pressure) * 0.5) 
valid_pressure = pressure[start_idx:]

# DC成分（平均のズレ）を除去
valid_pressure = valid_pressure - np.mean(valid_pressure)

# ハニング窓をかける
window = np.hanning(len(valid_pressure))
valid_pressure_windowed = valid_pressure * window

# FFT計算
N = len(valid_pressure)
freq = np.fft.rfftfreq(N, d=dt)
fft_val = np.fft.rfft(valid_pressure_windowed)
amplitude = np.abs(fft_val) / N * 2
p0 = 20e-6 
db_amplitude = 20 * np.log10(amplitude / p0)

# 1. 探索範囲 (20Hz以上)
mask = freq > 20
masked_freq = freq[mask]
masked_amp = amplitude[mask]

# 2. 「山（ピーク）」をすべて見つける
# 条件: 左隣より大きく、かつ右隣より大きい点
if len(masked_amp) > 2:
    is_peak = (masked_amp[1:-1] > masked_amp[:-2]) & (masked_amp[1:-1] > masked_amp[2:])
    # 配列サイズを合わせるため前後にFalseを追加
    is_peak = np.r_[False, is_peak, False]
    
    peak_indices = np.where(is_peak)[0]
    peak_freqs = masked_freq[peak_indices]
    peak_amps = masked_amp[peak_indices]
    
    # 3. 閾値フィルタ (最大ピークの10%以上あるものだけ残す)
    max_amp = np.max(masked_amp)
    threshold = 0.1 * max_amp 
    
    significant_peaks_idx = np.where(peak_amps > threshold)[0]
    
    if len(significant_peaks_idx) > 0:
        # 4. その中で「一番低い周波数」を選ぶ
        first_peak_idx = significant_peaks_idx[0]
        f0 = peak_freqs[first_peak_idx]
        f0_amp = peak_amps[first_peak_idx]
    else:
        # 万が一閾値を超えるものがない場合は最大値を使う
        idx = np.argmax(masked_amp)
        f0 = masked_freq[idx]
        f0_amp = masked_amp[idx]
else:
    f0 = 0
    f0_amp = 0

print(f"Detected Fundamental Frequency (F0): {f0:.2f} Hz")

fig, ax = plt.subplots(figsize=(10, 6))

# FFT結果のプロット
ax.plot(freq, db_amplitude)

# タイトルやラベルの設定
ax.set_title(f"Frequency Domain (F0 = {f0:.2f} Hz)")
ax.set_xlabel("Frequency [Hz]")
ax.set_ylabel("SPL [dB]")

# 見たい範囲（0Hz ～ 2000Hz）
ax.set_xlim(0, 6000) 

# グリッド表示
ax.grid(which='both', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()