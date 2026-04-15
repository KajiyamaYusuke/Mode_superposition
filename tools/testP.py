import os
import subprocess
import numpy as np
import optuna
import matplotlib.pyplot as plt

# ==========================================
# 設定
# ==========================================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE_PATH = os.path.join(BASE_DIR, "build", "simulation")
PARAM_FILE = os.path.join(BASE_DIR, "input", "param.txt")
AREA_FILE = os.path.join(BASE_DIR, "output", "airflow_vt.dat")
DB_URL = "sqlite:///optuna_study.db"
STUDY_NAME = "vocal_fold_optimization"
SIM_DT = 1.0e-5

def update_param_file(filepath, p_dict):
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

# 1. DBから一番良さそうなパレート解を1つ取得
study = optuna.load_study(study_name=STUDY_NAME, storage=DB_URL)
best_trials = study.best_trials

if len(best_trials) == 0:
    print("No successful trials found in DB.")
    exit()

# 適当に最初のパレート解を選択
target_trial = best_trials[0]
p = target_trial.params

print("=== Re-running Simulation with Pareto-best parameters ===")
for k, v in p.items():
    print(f"  {k}: {v:.4f}")

# 2. パラメータを書き換えて実行
update_param_file(PARAM_FILE, p)
print("Running C++ simulation...")
subprocess.run([EXE_PATH], cwd=os.path.join(BASE_DIR, "tools"), stdout=subprocess.DEVNULL)

# 3. 結果の読み込みとプロット
if not os.path.exists(AREA_FILE):
    print("Error: Output file not found.")
    exit()

data = np.loadtxt(AREA_FILE)
time = data[:, 0] * SIM_DT
flow = data[:, 1]

plt.figure(figsize=(10, 4))
plt.plot(time, flow, label="Glottal Flow", color="blue")
plt.axvline(x=0.1, color='red', linestyle='--', label="Analysis Start (0.1s)")
plt.xlabel("Time [s]")
plt.ylabel("Airflow [m^3/s]")
plt.title(f"Waveform of Pareto Solution (Trial {target_trial.number})")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("waveform_check.png", dpi=300)
print("\nDone! Please open 'waveform_check.png' to visually inspect the phonation.")