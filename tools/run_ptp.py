import os
import subprocess
import shutil
import numpy as np
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks, get_window

# ==========================================
# 1. 設定
# ==========================================
EXE_PATH = "../build/simulation"   # 実行ファイルのパス
PARAM_FILE = "../input/param.txt"
OUTPUT_DIR = "../output"
TARGET_FILE = "pressure_vt.dat"    # 解析対象

# 探索設定
PS_MIN = 200.0     # 探索範囲の下限 [Pa]
PS_MAX = 1000.0    # 探索範囲の上限 [Pa]
TOLERANCE = 10.0   # 探索精度 [Pa]

# 判定条件パラメータ
START_TIME = 0.1  # 解析開始時間 [s] (過渡応答無視)
MIN_AMP_PA = 10.0   # 【条件1】最小圧力振幅 [Pa]
MIN_FREQ = 50.0     # ノイズ/DC成分除外用の最低周波数 [Hz]

# ==========================================
# 2. 関数定義
# ==========================================

def get_f0_robust(xf, yf, min_freq=50, amp_threshold_ratio=0.2):
    """
    倍音の方が強くてもF0を正しく検出する関数
    
    xf: 周波数軸配列
    yf: 振幅スペクトル配列
    min_freq: これ以下の周波数は無視（ノイズ対策）
    amp_threshold_ratio: 最大振幅の何割以上を「主要なピーク」とみなすか (例: 0.2 = 20%)
    """
    # 1. 最大振幅で正規化
    max_amp = np.max(yf)
    if max_amp == 0: return 0.0
    yf_norm = yf / max_amp
    
    # 2. ピークを検出
    # height: 最大値の20%以上の高さがあるピークだけ抽出
    peaks, _ = find_peaks(yf_norm, height=amp_threshold_ratio)
    
    # 3. ピークの周波数をリスト化
    peak_freqs = xf[peaks]
    
    # 4. 指定した最低周波数 (min_freq) 以下のノイズを除外
    valid_freqs = peak_freqs[peak_freqs > min_freq]
    
    if len(valid_freqs) > 0:
        # 5. 有効なピークの中で「最も周波数が低いもの」をF0とする
        f0 = np.min(valid_freqs)
        return f0
    else:
        # ピークが見つからない場合は、仕方なく最大値の周波数を返すか0を返す
        return xf[np.argmax(yf)]

def update_pressure(filepath, ps_value):
    """param.txt の Ps を書き換える"""
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    new_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        new_lines.append(line)
        if "# subglottal pressure Ps" in line:
            if i + 1 < len(lines):
                new_lines.append(f"{ps_value:.2f}\n")
                i += 2
                continue
        i += 1
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)

def analyze_oscillation(filepath, dt=1.0e-5):
    """
    発声判定を行う関数
    戻り値: (判定結果bool, 振幅, F0周波数)
    """
    if not os.path.exists(filepath):
        return False, 0.0, 0.0

    try:
        data = np.loadtxt(filepath, skiprows=1)
        if data.ndim == 1 or data.shape[0] < 100:
            return False, 0.0, 0.0

        time_seq = data[:, 0] * dt
        pressure = data[:, 1]
        
        # 1. 解析区間の抽出（後半の安定部）
        mask = time_seq >= START_TIME
        if np.sum(mask) == 0: mask = time_seq >= (time_seq[-1] * 0.5)
        
        y_seg = pressure[mask]
        
        # --- 【条件1】時間領域の振幅チェック ---
        amplitude = (np.max(y_seg) - np.min(y_seg)) / 2.0
        if amplitude < MIN_AMP_PA:
            # 振幅が小さすぎる時点で不合格
            return False, amplitude, 0.0
            
        # --- 【条件2】FFTによるピーク(F0)チェック ---
        N = len(y_seg)
        
        # 直流成分除去 & 窓関数適用
        y_detrend = y_seg - np.mean(y_seg)
        window = get_window('hann', N)
        y_windowed = y_detrend * window
        
        # FFT計算
        yf = np.abs(fft(y_windowed))[:N//2]
        xf = fftfreq(N, d=dt)[:N//2]
        
        # 振幅スペクトルを正規化（最大値を1に）
        if np.max(yf) == 0: return False, amplitude, 0.0
        yf_norm = yf / np.max(yf)
        
        # ピーク検出 (scipy.signal.find_peaks)
        # prominence: 周辺よりどれだけ突出しているか (0.1 = 最大値の10%以上の突出が必要)
        # width: ノイズの鋭いスパイクを除外するために幅を指定してもよい
        peaks, properties = find_peaks(yf_norm, prominence=0.1)
        
        valid_peaks = []
        for p in peaks:
            freq = xf[p]
            if freq >= MIN_FREQ: # DCや低周波ドリフトを除外
                valid_peaks.append((freq, yf_norm[p]))
        
        # 有効なピークがなければ発声なし（ただのノイズ）
        if not valid_peaks:
            return False, amplitude, 0.0
            
        f0_freq = get_f0_robust(xf, yf)
        
        # ここまで来たら「振幅十分」かつ「明確なピークあり」
        return True, amplitude, f0_freq

    except Exception as e:
        print(f"  [Analysis Error] {e}")
        return False, 0.0, 0.0

# ==========================================
# 3. メイン処理（二分探索）
# ==========================================
if os.path.exists(PARAM_FILE):
    shutil.copy(PARAM_FILE, PARAM_FILE + ".bak")

try:
    low = PS_MIN
    high = PS_MAX
    best_ptp = None
    best_freq = 0.0
    
    print(f"--- PTP Search with FFT Analysis ---")
    print(f"Conditions: Amp >= {MIN_AMP_PA} Pa AND Clear Peak >= {MIN_FREQ} Hz")
    print(f"Range: {low} - {high} Pa")
    
    iteration = 1
    
    while (high - low) > TOLERANCE:
        mid = (low + high) / 2.0
        print(f"\n[Iter {iteration}] Testing Ps = {mid:.2f} Pa ... ", end="", flush=True)
        
        # 1. パラメータ更新
        update_pressure(PARAM_FILE, mid)
        
        # 2. シミュレーション実行 (エラー時は表示されるよう設定)
        ret = subprocess.run([EXE_PATH], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        if ret.returncode != 0:
            print("Simulation Failed.")
            break
            
        # 3. 解析
        dat_path = os.path.join(OUTPUT_DIR, TARGET_FILE)
        is_osc, amp, freq = analyze_oscillation(dat_path)
        
        status = "OSC" if is_osc else "NO"
        print(f"-> {status} (Amp: {amp:.1f} Pa, F0: {freq:.1f} Hz)")
        
        if is_osc:
            # 発声成功 -> もっと低い圧力を試す
            high = mid
            best_ptp = mid
            best_freq = freq
        else:
            # 発声失敗 -> 圧力を上げる
            low = mid
            
        iteration += 1

    print("\n==========================================")
    if best_ptp is not None:
        print(f" Result: PTP is approx. {best_ptp:.2f} Pa")
        print(f" (F0 at threshold: {best_freq:.1f} Hz)")
    else:
        print(f" Failed: No oscillation found within range.")
    print("==========================================")

finally:
    if os.path.exists(PARAM_FILE + ".bak"):
        shutil.move(PARAM_FILE + ".bak", PARAM_FILE)