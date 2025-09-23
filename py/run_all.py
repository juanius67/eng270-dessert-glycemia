import os, yaml
from desserts import desserts
from interface import run_simulation
from plotting import plot_results
from metrics import compute_metrics

HERE = os.path.dirname(__file__)
CONFIG = os.path.join(HERE, "..", "configs", "params.yaml")

def main():
    with open(CONFIG) as f:
        base = yaml.safe_load(f)

    # reset summary file
    summ = os.path.join(HERE, "..", "tables", "summary.csv")
    if os.path.exists(summ): os.remove(summ)

    for name, (A,k) in desserts.items():
        print(f"Simulating {name}...")
        res = run_simulation(base, A, k)
        plot_results(name, res)
        compute_metrics(name, res)

    print("Done. See figures/ and tables/.")

if __name__ == "__main__":
    main()
