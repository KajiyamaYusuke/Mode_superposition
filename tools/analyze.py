import os

import subprocess

import numpy as np

import optuna

import shutil

from scipy.fft import fft, fftfreq

from scipy.signal import find_peaks, get_window

import matplotlib.pyplot as plt



# OptunaのMatplotlibバックエンドを使用して図を出力

from optuna.visualization.matplotlib import plot_param_importances

from optuna.visualization.matplotlib import plot_slice



# ==========================================

# 1. 設定・目標値 (実験データのターゲット)

# ==========================================

EXE_PATH = os.path.abspath("../build/simulation")



# コピー元となる input フォルダの場所

BASE_INPUT_DIR = "../input"



# 並列用の作業フォルダ群を置くベースディレクトリ（★これが抜けてエラーになりました）

WORKSPACE_BASE = "../workspaces"

PARAM_FILE = "../input/param.txt"

AREA_FILE = "../output/airflow_vt.dat"

DT = 1.0e-5



TARGET_F0 = 132.9      # [Hz]

TARGET_JITTER = 2.6     # [%]

TARGET_CQ = 0.182        # 閉鎖率 (Closed Quotient)



# 損失関数の重み

W_F0 = 1.0

W_JITTER = 0.5

W_CQ = 1.0



MIN_AMP_AREA = 1e-6



# ==========================================

# 2. 関数定義

# ==========================================

def update_param_file(filepath, p_dict):

    """パラメータを書き換える"""

    with open(filepath, 'r', encoding='utf-8') as f:

        lines = f.readlines()

   

    new_lines = []

    i = 0

    while i < len(lines):

        line = lines[i]

        new_lines.append(line)

        if "# damping coefficient zeta" in line:

            if i + 1 < len(lines):

                new_lines.append(f"{p_dict['zeta']:.4f}\n")

                i += 2; continue

        elif "# contact coefficients kc1" in line:

            if i + 1 < len(lines):

                new_lines.append(f"{p_dict['kc1']:.1f}    {p_dict['kc2']:.1f}    {p_dict['kc3']:.1f}\n")

                i += 2; continue

        elif "# subglottal pressure Ps" in line:

            if i + 1 < len(lines):

                new_lines.append(f"{p_dict['ps']:.1f}\n")

                i += 2; continue

        i += 1

    with open(filepath, 'w', encoding='utf-8') as f:

        f.writelines(new_lines)



def calc_f0_autocorr(y, dt, fmin=50, fmax=500):

    y = y - np.mean(y)

    corr = np.correlate(y, y, mode='full')

    corr = corr[len(corr)//2:]



    lag_min = int(1/(fmax*dt))

    lag_max = int(1/(fmin*dt))



    peak = np.argmax(corr[lag_min:lag_max]) + lag_min

    f0 = 1.0 / (peak * dt)

    return f0



def calc_jitter_local_area(y, dt, t,

                           start_time=0.05,

                           window_samples=5,

                           min_mean_amp=1e-6):

    """

    局所移動平均との差から正規化ジッターを計算

    """



    # 指定時刻以降を抽出

    mask = t > start_time

    if np.sum(mask) < window_samples:

        return None



    y_sub = y[mask]



    mean_amp = np.mean(y_sub)

    if mean_amp < min_mean_amp:

        return None  # 発声していない



    # 移動平均（畳み込みで高速化）

    kernel = np.ones(window_samples) / window_samples

    local_avg = np.convolve(y_sub, kernel, mode="same")



    diff = np.abs(y_sub - local_avg)



    jitter = np.mean(diff) / mean_amp



    return jitter



def analyze_waveform(filepath, dt):



    if not os.path.exists(filepath):

        return None



    try:

        data = np.loadtxt(filepath, skiprows=1)

        if data.ndim == 1 or data.shape[0] < 100:

            return None



        time = data[:, 0] * dt

        y = data[:, 1]



        dt_out = time[1] - time[0]



        # 0.1秒以降を安定区間とする

        mask = time >= 0.1

        if np.sum(mask) == 0:

            mask = time >= (time[-1] * 0.5)



        t_seg = time[mask]

        y_seg = y[mask]



        # 発声判定

        flow_max = np.max(y_seg)
        flow_min = np.min(y_seg)
        
        # 1. 振幅のチェック（ノイズ弾き）
        amp = (flow_max - flow_min) / 2.0
        if amp < MIN_AMP_FLOW:
            return None
            
        # 2. ★新規追加：「開きっぱなし」チェック
        # 最小流量が最大流量の 50% 以上ある場合は、声帯が衝突せずにヒラヒラしているだけとみなす
        if flow_min > flow_max * 0.5:
            return None



        # ===== F0（自己相関推奨） =====

        f0 = calc_f0_autocorr(y_seg, dt_out)

        if not (50 <= f0 <= 500):

            return None



        # ===== 新ジッター =====

        jitter = calc_jitter_local_area(

            y_seg, dt_out, t_seg,

            start_time=0.05,

            window_samples = max(5, int((1/f0) / dt * 0.05))

        )



        if jitter is None:

            return None



        # ===== CQ =====

        close_threshold = 1e-8  # 絶対値基準

        closed_time = np.sum(y_seg < close_threshold) * dt_out

        total_time = len(y_seg) * dt_out

        cq = closed_time / total_time if total_time > 0 else 0.0



        return {"F0": f0, "Jitter": jitter, "CQ": cq}



    except Exception:

        return None



# ==========================================

# 3. Optuna の目的関数

# ==========================================

def objective(trial):

    trial_id = trial.number

    work_dir = os.path.join(WORKSPACE_BASE, f"trial_{trial_id}")

   

    # 必要な4つのディレクトリのパスを定義

    input_dir = os.path.join(work_dir, "input")

    output_dir = os.path.join(work_dir, "output")

    tools_dir = os.path.join(work_dir, "tools")

   

    # ディレクトリを一気に作成

    os.makedirs(input_dir, exist_ok=True)

    os.makedirs(output_dir, exist_ok=True)

    os.makedirs(tools_dir, exist_ok=True)



    # 探索するパラメータの定義

    p = {

        "kc1": trial.suggest_float("kc1", 1.0, 100.0, log=True),

        "kc2": trial.suggest_float("kc2", 100.0, 10000.0, log=True),

        "kc3": trial.suggest_float("kc3", 100.0, 10000.0, log=True),

        "zeta": trial.suggest_float("zeta", 0.005, 0.5),

        "ps": trial.suggest_float("ps", 500.0, 2000.0)

    }



    # 1. param.txt の実体コピーと書き換え

    param_src = os.path.join(BASE_INPUT_DIR, "param.txt")

    param_dst = os.path.join(input_dir, "param.txt")

    shutil.copy(param_src, param_dst)

    update_param_file(param_dst, p)

   

    # 2. 巨大な「M5」フォルダは、シンボリックリンクを作成

    m5_src = os.path.abspath(os.path.join(BASE_INPUT_DIR, "M5"))

    m5_dst = os.path.join(input_dir, "M5")

    if not os.path.exists(m5_dst):

        os.symlink(m5_src, m5_dst)



    # 3. C++の実行

    ret = subprocess.run(

        [EXE_PATH],

        cwd=tools_dir,           # ★重要: tools フォルダ内で実行する！

        stdout=subprocess.DEVNULL,

        stderr=subprocess.DEVNULL,     # ★重要: 出力を取得する設定

    )

   

    if ret.returncode != 0:

        print(f"\n[{trial_id}] C++ Error! Return code: {ret.returncode}")

        #stderr_msg = ret.stderr if ret.stderr else "No STDERR output"
        shutil.rmtree(work_dir, ignore_errors=True)
        print(f"[{trial_id}] STDERR: {stderr_msg}")

        return 10000.0, 10000.0, 10000.0

       

    area_file = os.path.join(output_dir, "airflow_vt.dat")

    res = analyze_waveform(area_file, DT)

   

    if res is None:

        # 発声しなかった場合も3つの大きなペナルティを返す
        shutil.rmtree(work_dir, ignore_errors=True)
        return 5000.0, 5000.0, 5000.0

       

    err_f0 = abs(np.log(res["F0"]/TARGET_F0))

    err_jitter = abs(res["Jitter"]*100 - TARGET_JITTER) / (TARGET_JITTER + 1e-6)

    err_cq = abs(res["CQ"] - TARGET_CQ) / (TARGET_CQ + 1e-6)

   

    # 記録

    trial.set_user_attr("F0", res["F0"])

    trial.set_user_attr("Jitter", res["Jitter"])

    trial.set_user_attr("CQ", res["CQ"])



    shutil.rmtree(work_dir, ignore_errors=True)

   

    return err_f0, err_jitter, err_cq



# ==========================================

# 4. メイン処理 (最適化 & 図の出力)

# ==========================================

if __name__ == "__main__":

    print("--- Starting Optimization (Parallel Ready) ---")

   

    os.makedirs(WORKSPACE_BASE, exist_ok=True)

   

    # ★変更点: SQLiteデータベースを指定してStudyを作成/読み込み

    study_name = "vocal_fold_optimization"

    storage_url = "sqlite:///optuna_study.db"

   

    study = optuna.create_study(

        study_name=study_name,

        storage=storage_url,

        load_if_exists=True,  # 既にDBがあれば続きから始める

        directions=["minimize","minimize","minimize"],

        sampler=optuna.samplers.TPESampler(n_startup_trials=20)

    )

   

    # このスクリプト1つあたり何回計算するか (例: 50回)

    # 4プロセス立ち上げれば、全体で約200回計算が進みます。

    study.optimize(objective, n_trials=50)

   

    print("\n=== Finished this process ===")



    # ==========================================

    # 5. 感度解析の図を出力 (PNG保存)

    # ==========================================

    print("\nSaving Sensitivity Analysis Plots...")

   

# 例：F0誤差 (目的関数の0番目) に対する感度を出力する

    fig1 = plot_param_importances(study, target=lambda t: t.values[0])

    plt.title("Parameter Importances for F0 Error")

    plt.tight_layout()

    plt.savefig("sensitivity_importance_F0.png", dpi=300)



    # ===== Jitter =====

    plt.figure()

    plot_param_importances(study, target=lambda t: t.values[1])

    plt.title("Parameter Importances for Jitter Error")

    plt.tight_layout()

    plt.savefig("sensitivity_importance_Jitter.png", dpi=300)

    plt.close()



    plt.figure()

    plot_slice(study, target=lambda t: t.values[1])

    plt.title("Parameter Slice Plot for Jitter Error")

    plt.tight_layout()

    plt.savefig("sensitivity_slice_Jitter.png", dpi=300)

    plt.close()



    # ===== Pareto =====

    plt.figure()

    optuna.visualization.matplotlib.plot_pareto_front(study)

    plt.tight_layout()

    plt.savefig("pareto_front.png", dpi=300)

    plt.close()

   

    # 例：CQ誤差 (目的関数の2番目) に対するスライス図を出力する

    fig2 = plot_slice(study, target=lambda t: t.values[2])

    plt.title("Parameter Slice Plot for CQ Error")

    plt.tight_layout()

    plt.savefig("sensitivity_slice_CQ.png", dpi=300)



    optuna.visualization.matplotlib.plot_pareto_front(study)

   

    print("Done. Check 'sensitivity_importance.png' and 'sensitivity_slice.png'.")