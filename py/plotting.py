import os
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
OUTDIR = os.path.abspath(os.path.join(HERE, "..", "figures"))

def plot_results(name, results, outdir=OUTDIR):
    os.makedirs(outdir, exist_ok=True)
    t, G, I = results["t"], results["G"], results["I"]

    # Glucose
    plt.figure()
    plt.plot(t, G, label=name)
    plt.xlabel("Time (min)"); plt.ylabel("Glucose (mg/dL)")
    plt.title(f"Glucose response: {name}")
    plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{name}_glucose.png"))
    plt.close()

    # Insulin
    plt.figure()
    plt.plot(t, I, label=name)
    plt.xlabel("Time (min)"); plt.ylabel("Insulin (µU/mL)")
    plt.title(f"Insulin response: {name}")
    plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{name}_insulin.png"))
    plt.close()
