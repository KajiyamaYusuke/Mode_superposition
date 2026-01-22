import numpy as np
import matplotlib.pyplot as plt
# 1. datファイルを読み込む（空白区切りやタブ区切りを想定）
harea = np.loadtxt("../output/airflow_vt.dat")

labels = harea[:, 0]        # ← 1列目（横軸のラベル）
values = harea[:, 1]  

# 2. 各行の最小値を求める
#row_min = np.min(values, axis=1)  # 各行の最小値を1次元配列で返す

x = labels * 1e-5
start_time = 0.15           # 切り出したい開始時刻 (秒) ※波形を見て調整してください
duration   = 0.02

plt.figure(figsize=(15,4))
plt.plot(x, values, linestyle='-', color='g')  # ← マーカーなし線だけ
plt.xlabel('time')
plt.ylabel('airflow')
plt.title('T6')
plt.grid(True)
plt.xlim(start_time, start_time + duration)
plt.ylim(-0.001, 0.001)
plt.tight_layout()
plt.legend()
plt.show()
