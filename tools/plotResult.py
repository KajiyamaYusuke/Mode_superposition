import optuna
import matplotlib.pyplot as plt
from optuna.visualization.matplotlib import plot_param_importances, plot_slice

# 1. 元のStudyを読み込む
study = optuna.load_study(
    study_name="vocal_fold_optimization",
    storage="sqlite:///optuna_study.db"
)

# 2. 成功したTrialだけを抽出した「クリーンなStudy」をメモリ上に作成する
clean_study = optuna.create_study(directions=["minimize", "minimize", "minimize"])
for trial in study.trials:
    # values が存在し、かつペナルティ値(5000など)ではないTrialのみ追加
    if trial.values is not None and all(v < 1000.0 for v in trial.values):
        clean_study.add_trial(trial)

print(f"Filtered Study: {len(clean_study.trials)} successful trials out of {len(study.trials)}.")

# 3. クリーンなデータで感度解析とスライスプロットを出力
target_names = ["F0 Error", "Jitter Error", "CQ Error"]

for i, name in enumerate(target_names):
    # 感度（重要度）のプロット
    plt.figure()
    plot_param_importances(clean_study, target=lambda t: t.values[i], target_name=name)
    plt.title(f"Parameter Importances for {name} (Clean Data)")
    plt.tight_layout()
    # ファイル名からスペースを消して保存
    filename_imp = f"clean_importance_{name.replace(' ', '')}.png"
    plt.savefig(filename_imp, dpi=300)
    plt.close()

    # スライスプロット（傾向の散布図）
    plt.figure()
    plot_slice(clean_study, target=lambda t: t.values[i], target_name=name)
    plt.title(f"Parameter Slice Plot for {name} (Clean Data)")
    plt.tight_layout()
    filename_slice = f"clean_slice_{name.replace(' ', '')}.png"
    plt.savefig(filename_slice, dpi=300)
    plt.close()

print("Done! Check the 'clean_***.png' files.")