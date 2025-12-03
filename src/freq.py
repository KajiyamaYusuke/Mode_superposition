import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt
from scipy.fft import fft, fftfreq

# 1. データ読み込み
harea = np.loadtxt("../output/area.dat")
labels = harea[:, 0]
values = harea[:, 1:]

x = labels * 1e-5  # 秒
row_min = np.min(values, axis=1)

# ---- ★ ここで解析区間を限定 ★ ----
t_start = 0.15
t_end   = 0.20
idx = np.where((x >= t_start) & (x <= t_end))[0]

x_seg = x[idx]
row_seg = row_min[idx]

# ---- サンプリング周波数 ----
fs = 1 / (x_seg[1] - x_seg[0])
nyq = fs / 2

# 2. ノイズ除去（低域通過フィルタ例）
cutoff = 1000  # Hz
b, a = butter(4, cutoff / nyq, btype='low')
row_seg_smooth = filtfilt(b, a, row_seg)

# 3. FFT（区間データで）
N = len(row_seg_smooth)
yf = fft(row_seg_smooth)
xf = fftfreq(N, d=1/fs)

# 正の周波数だけ
pos = xf > 0
xf = xf[pos]
yf = np.abs(yf[pos])

peak_freq = xf[np.argmax(yf)]
print(f"[0.15-0.20 s] FFT主要周波数: {peak_freq:.2f} Hz")

# 4. ピーク検出
peaks, _ = find_peaks(row_seg_smooth, prominence=0.001)
peak_times = x_seg[peaks]

if len(peak_times) >= 2:
    periods = np.diff(peak_times)
    frequency = 1 / np.mean(periods)
    print(f"[0.15-0.20 s] ピーク間平均周波数: {frequency:.2f} Hz")
else:
    print("⚠️ 区間内でピークが十分検出できませんでした。")

# ---- 可視化 ----
plt.figure(figsize=(12,5))
plt.plot(x_seg, row_seg, alpha=0.5, label="raw (seg)")
plt.plot(x_seg, row_seg_smooth, label="smoothed (seg)")
plt.plot(peak_times, row_seg_smooth[peaks], "ro", label="peaks")
plt.xlabel("time [s]")
plt.ylabel("area")
plt.title("0.15–0.20 s glottal area")
plt.grid(True)
plt.legend()
plt.show()

# FFTプロット
plt.figure(figsize=(10,4))
limit = 1500
mask = xf <= limit

plt.plot(xf[mask], yf[mask])
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude")
plt.title("FFT (0.15–0.20 s)")
plt.grid(True)
plt.xlim(0,limit)
plt.show()
