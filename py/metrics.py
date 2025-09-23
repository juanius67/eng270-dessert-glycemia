import numpy as np, os

HERE = os.path.dirname(__file__)
OUTDIR = os.path.abspath(os.path.join(HERE, "..", "tables"))

def compute_metrics(name, results, outdir=OUTDIR):
    t, G, I = results["t"], results["G"], results["I"]
    dt = t[1]-t[0]
    Gb = G[0]

    peak_G = G.max()
    t_peak_G = t[G.argmax()]
    auc_G = np.trapz(np.clip(G-Gb, 0, None), t)

    peak_I = I.max()
    t_peak_I = t[I.argmax()]
    auc_I = np.trapz(np.clip(I-I[0], 0, None), t)

    os.makedirs(outdir, exist_ok=True)
    line = f"{name},{peak_G:.2f},{t_peak_G:.1f},{auc_G:.1f},{peak_I:.2f},{t_peak_I:.1f},{auc_I:.1f}\n"
    header = "dessert,peak_G(t0 units),t_peak_G(min),AUC_G,peak_I,t_peak_I(min),AUC_I\n"
    path = os.path.join(outdir, "summary.csv")
    if not os.path.exists(path):
        with open(path, "w") as f: f.write(header)
    with open(path, "a") as f: f.write(line)
    print(f"[metrics] {name}: peak_G={peak_G:.1f} at {t_peak_G:.1f} min; AUC_G={auc_G:.1f}")
