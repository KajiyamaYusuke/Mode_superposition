import numpy as np
import matplotlib.pyplot as plt

# =========================
# 設定
# =========================
filename = "../output/pressure_vt.dat"   # datファイル名
dt = 1.0e-5                 # 1ステップあたりの時間 [s] ← 必ず調整！

# =========================
# データ読み込み
# =========================
data = np.loadtxt(filename)

step = data[:, 0]
pressure = data[:, 1]

# DC成分除去（重要）
pressure = pressure - np.mean(pressure)

# =========================
# FFT
# =========================
N = len(pressure)
freq = np.fft.rfftfreq(N, d=dt)
fft_val = np.fft.rfft(pressure)

amplitude = np.abs(fft_val) / N * 2   # 片側スペクトル

# =========================
# プロット
# =========================
plt.figure()
plt.plot(freq, amplitude)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude")
plt.title("Pressure FFT")
plt.xlim(0, 2000)   # 必要に応じて
plt.grid(True)
plt.show()
