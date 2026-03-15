import sys
import csv
import matplotlib.pyplot as plt
import subprocess

def load_csv(path):
    k, E, F = [], [], []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            k.append(int(row["k"]))
            E.append(float(row["energy"]))
            F.append(float(row["fidelity"]))
    return k, E, F

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_matrix.py data.csv")
        sys.exit(1)
    print("Plot")
    csv_path = sys.argv[1]
    k, E, F = load_csv(csv_path)

    plt.figure()
    plt.plot(k, E, marker="o")
    plt.xlabel("k")
    plt.ylabel("Energy <X>")
    plt.grid()
    plt.title("DB-QITE (matrix) Energy")
    out1 = "matrix_energy.png"
    plt.savefig(out1, dpi=200, bbox_inches="tight")

    plt.figure()
    plt.plot(k, F, marker="o")
    plt.xlabel("k")
    plt.ylabel("Fidelity to |->")
    plt.grid()
    plt.title("DB-QITE (matrix) Fidelity")
    out2 = "matrix_fidelity.png"
    plt.savefig(out2, dpi=200, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()