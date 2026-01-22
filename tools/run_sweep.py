import os
import subprocess
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt

# ==========================================
# 1. 設定
# ==========================================
EXE_PATH = "../build/simulation"  # 実行ファイルのパス
PARAM_FILE = "../input/param.txt"
OUTPUT_DIR = "../output"
AREA_FILE = "area.dat"

# 解析パラメータ
DT = 1.0e-5
START_TIME = 0.1
CUTOFF_FREQ = 500

# --- 基準パラメータ ---
# kc=[kc1, kc2, kc3], zeta固定, ps基準値
base_params = {"kc": [10, 1000, 3000], "zeta": 0.1, "ps": 1600}

experiments = {}

# 1. Baseline（基準）
experiments["Baseline"] = base_params.copy()

# 2. 駆動圧 (Ps) の変動
ps_range = [1000, 1200, 1400, 1800, 2000, 2200]
for p in ps_range:
    new_p = base_params.copy()
    new_p["ps"] = p
    experiments[f"Ps_{p}"] = new_p

# 3. 接触剛性 (kc) の個別変動
kc_labels = ["kc1", "kc2", "kc3"]
kc_scales = [0.1, 0.5, 2.0, 10.0]

for i, label in enumerate(kc_labels):
    for scale in kc_scales:
        new_p = base_params.copy()
        current_kc = list(base_params["kc"])
        current_kc[i] = int(current_kc[i] * scale)
        new_p["kc"] = current_kc
        experiments[f"{label}_x{scale}"] = new_p

# ==========================================
# 2. 関数定義
# ==========================================
def update_param_file(filepath, params):
    """param.txtを書き換える"""
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    new_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        new_lines.append(line)
        if "# damping coefficient zeta" in line:
            if i + 1 < len(lines):
                new_lines.append(f"{params['zeta']}\n")
                i += 2; continue
        elif "# contact coefficients kc1" in line:
            if i + 1 < len(lines):
                k = params['kc']
                new_lines.append(f"{k[0]}    {k[1]}    {k[2]}\n")
                i += 2; continue
        elif "# subglottal pressure Ps" in line:
            if i + 1 < len(lines):
                new_lines.append(f"{params['ps']}\n")
                i += 2; continue
        i += 1
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)

def analyze_waveform(filepath, dt, start_time):
    """波形解析を行う関数"""
    if not os.path.exists(filepath):
        return None

    try:
        # データの読み込み
        data = np.loadtxt(filepath, skiprows=1)
        # 2列目以降の最小値を声門面積とする
        waveform = np.min(data[:, 1:], axis=1)
        time = data[:, 0] * dt
        
        # 解析区間の抽出（ここから下のインデントに注意）
        mask = time >= start_time
        if np.sum(mask) < 100: mask = time >= 0
        t_seg = time[mask]
        y_seg = waveform[mask]
        
        # フィルタ処理
        fs = 1.0 / dt
        nyq = fs / 2.0
        b, a = butter(4, CUTOFF_FREQ / nyq, btype='low')
        y_smooth = filtfilt(b, a, y_seg)
        
        # ピーク検出
        peaks, _ = find_peaks(y_smooth, distance=int(fs/1000))
        
        if len(peaks) < 2:
            freq = 0.0
            jitter = 0.0
        else:
            periods = np.diff(t_seg[peaks])
            avg_per = np.mean(periods)
            freq = 1.0 / avg_per if avg_per > 0 else 0
            jitter = (np.mean(np.abs(np.diff(periods))) / avg_per) * 100 if len(periods) > 1 else 0

        # Gitter (Noise)
        y_s = pd.Series(y_seg)
        mean_area = np.mean(y_seg)
        if mean_area > 1e-12:
            local_avg = y_s.rolling(5, center=True).mean().fillna(y_s)
            diff_abs = np.abs(y_s - local_avg)
            gitter = (diff_abs.mean() / mean_area) * 100
        else:
            gitter = 0.0

        return {"Frequency [Hz]": freq, "Jitter [%]": jitter, "Gitter [%]": gitter}

    except Exception as e:
        print(f"Analysis Error: {e}")
        return None

# ==========================================
# 3. メイン処理
# ==========================================
results = []

# パラメータファイルのバックアップ
if os.path.exists(PARAM_FILE):
    shutil.copy(PARAM_FILE, PARAM_FILE + ".bak")

try:
    print(f"--- Starting Independent Sensitivity Analysis ({len(experiments)} cases) ---")
    
    for case, params in experiments.items():
        print(f"Running: {case} ... ", end="", flush=True)
        
        # パラメータ更新
        update_param_file(PARAM_FILE, params)
        
        # シミュレーション実行
        ret = subprocess.run([EXE_PATH], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if ret.returncode != 0:
            print("Failed.")
            continue
            
        # 解析
        dat = os.path.join(OUTPUT_DIR, AREA_FILE)
        res = analyze_waveform(dat, DT, START_TIME)
        
        if res:
            res["Case"] = case
            # グループ分け
            if "Ps" in case: res["Group"] = "Ps"
            elif "kc1" in case: res["Group"] = "kc1"
            elif "kc2" in case: res["Group"] = "kc2"
            elif "kc3" in case: res["Group"] = "kc3"
            else: res["Group"] = "Baseline"
            
            results.append(res)
            print(f"Freq: {res['Frequency [Hz]']:.1f} Hz")
        else:
            print("No output.")

finally:
    # パラメータファイルの復元
    if os.path.exists(PARAM_FILE + ".bak"):
        shutil.move(PARAM_FILE + ".bak", PARAM_FILE)

# ==========================================
# 4. 可視化
# ==========================================
if results:
    df = pd.DataFrame(results)
    # CSV保存
    csv_name = "sensitivity_independent.csv"
    df.to_csv(csv_name, index=False)
    print(f"\nSaved results to '{csv_name}'")
    
    # グラフ描画
    groups = ["Ps", "kc1", "kc2", "kc3"]
    fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharey=True)
    
    base_row = df[df["Case"]=="Baseline"]
    base_freq = base_row["Frequency [Hz]"].values[0] if not base_row.empty else 0
    
    for i, grp in enumerate(groups):
        sub_df = df[(df["Group"] == grp) | (df["Group"] == "Baseline")].copy()
        sub_df = sub_df.sort_values("Case")
        
        colors = ['gray' if c == "Baseline" else 'skyblue' for c in sub_df["Case"]]
        axes[i].bar(sub_df["Case"], sub_df["Frequency [Hz]"], color=colors)
        axes[i].axhline(base_freq, color='red', linestyle='--', alpha=0.5)
        
        axes[i].set_title(f"Effect of {grp}")
        axes[i].tick_params(axis='x', rotation=90)
        axes[i].grid(axis='y', linestyle='--', alpha=0.7)
    
    axes[0].set_ylabel("Frequency [Hz]")
    plt.tight_layout()
    plt.show()
else:
    print("No results to display.")