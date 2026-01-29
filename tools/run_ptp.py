import os
import subprocess
import shutil
import numpy as np

# ==========================================
# 1. 設定
# ==========================================
EXE_PATH = "../build/simulation"   # 実行ファイルのパス
PARAM_FILE = "../input/param.txt"
OUTPUT_DIR = "../output"
TARGET_FILE = "pressure_vt.dat"    # 解析対象ファイル (pressure_vt)

# 探索設定
PS_MIN = 200.0     # 探索範囲の下限 [Pa]
PS_MAX = 1000.0    # 探索範囲の上限 [Pa]
TOLERANCE = 10.0   # 探索を終了する精度（この幅以下になるまで絞り込む）

# 発声判定パラメータ
START_TIME = 0.1  # 判定開始時間 [s] (過渡応答を避ける)
AMP_THRESH = 50.0  # 発声とみなす振幅の閾値 [Pa] (これを超えたら発声)

# ==========================================
# 2. 関数定義
# ==========================================

def update_pressure(filepath, ps_value):
    """param.txt の Ps (声門下圧) を書き換える"""
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    new_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        new_lines.append(line)
        if "# subglottal pressure Ps" in line:
            if i + 1 < len(lines):
                # 小数点以下も書き込む
                new_lines.append(f"{ps_value:.2f}\n")
                i += 2
                continue
        i += 1
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)

def check_oscillation(filepath, dt=1.0e-5):
    """
    pressure_vt.dat を読み込み、振幅が基準を超えているか判定する
    """
    if not os.path.exists(filepath):
        return False, 0.0

    try:
        # データの読み込み
        # 1列目: ステップ数, 2列目: 圧力 (simulation.cppの出力形式に基づく)
        data = np.loadtxt(filepath, skiprows=1)
        
        if data.ndim == 1 or data.shape[0] < 100:
            return False, 0.0

        time_seq = data[:, 0] * dt
        pressure = data[:, 1]
        
        # 後半の安定部分のみ抽出
        mask = time_seq >= START_TIME
        if np.sum(mask) == 0: 
            # データが短い場合は後半50%を使う
            mask = time_seq >= (time_seq[-1] * 0.5)
            
        y_seg = pressure[mask]
        
        # 振幅計算: (最大値 - 最小値) / 2
        amplitude = (np.max(y_seg) - np.min(y_seg)) / 2.0
        
        # 判定
        is_oscillating = amplitude > AMP_THRESH
        
        return is_oscillating, amplitude

    except Exception as e:
        print(f"  [Error] Analysis failed: {e}")
        return False, 0.0

# ==========================================
# 3. メイン処理（二分探索）
# ==========================================
# バックアップ作成
if os.path.exists(PARAM_FILE):
    shutil.copy(PARAM_FILE, PARAM_FILE + ".bak")

try:
    low = PS_MIN    # 発声しない圧力（の下限）
    high = PS_MAX   # 発声する圧力（の上限）
    
    best_ptp = None # 見つかったPTP候補
    
    print(f"--- Starting PTP Search (Target: {TARGET_FILE}) ---")
    print(f"Range: {low} - {high} Pa, Threshold Amp: {AMP_THRESH} Pa")
    
    iteration = 1
    
    # 探索ループ
    while (high - low) > TOLERANCE:
        # 中間点を試す
        mid = (low + high) / 2.0
        
        print(f"\n[Iter {iteration}] Testing Ps = {mid:.2f} Pa ... ", end="", flush=True)
        
        # 1. パラメータ更新
        update_pressure(PARAM_FILE, mid)
        
        # 2. シミュレーション実行 (ログは非表示)
        ret = subprocess.run([EXE_PATH], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        if ret.returncode != 0:
            print("Simulation Failed.")
            break
            
        # 3. 判定
        dat_path = os.path.join(OUTPUT_DIR, TARGET_FILE)
        is_osc, amp_val = check_oscillation(dat_path)
        
        print(f"-> Amp = {amp_val:.2f} Pa [{'OSC' if is_osc else 'SILENCE'}]")
        
        if is_osc:
            # 発声した -> もっと低い圧力でもいけるかもしれない
            high = mid
            best_ptp = mid
        else:
            # 発声しなかった -> 圧力が足りない
            low = mid
            
        iteration += 1

    print("\n==========================================")
    if best_ptp is not None:
        print(f" Result: PTP is approx. {best_ptp:.2f} Pa")
        print(f" (Search ended with range {low:.2f} - {high:.2f} Pa)")
    else:
        print(f" Failed: No oscillation found even at {high} Pa.")
    print("==========================================")

finally:
    # パラメータを元に戻す
    if os.path.exists(PARAM_FILE + ".bak"):
        shutil.move(PARAM_FILE + ".bak", PARAM_FILE)