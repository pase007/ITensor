import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt


# ------------- Extraction of Data in filename (Latter maybe add headers?) --------------------
def extract_s_step(filename):
    match = re.search(r"data([0-9]*\.?[0-9]+)", filename)
    if match:
        return float(match.group(1))
    else:
        raise ValueError(f"Could not extract s_step from {filename}")

def extract_data_full(filename):
    match = re.search(r"data([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)", filename)
    if match:
        return float(match.group(1))
    else:
        raise ValueError(f"Could not extract s_step from {filename}")

# ---------------------- Read csv-files ----------------------
def read_csv(path):
    ks, E, F = [], [], []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            ks.append(int(row["k"]))
            E.append(float(row["Energy"]))
            F.append(float(row["Fidelity"]))
    return ks, E, F

def read_csv_Opt(path):
    s, k, If = [], [], []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            s.append(int(row["s"]))
            k.append(float(row["k"]))
            If.append(float(row["Infidelity"]))
    return s, k, If


# ----------------------------- Plot Data in different ways --------------
def plot_single_curves(out_png, csv_path, mode="F"):
    # Numerical data
    ks, E, F = read_csv(csv_path)
    If = []
    for f in F:
        If.append(1-f)

    # Analytic comparison
    #s_step = extract_s_step(csv_path)
    s_step = 0.6
    tau = np.array(ks) * s_step
    F_exact = 0.5 * (1.0 + np.tanh(2.0 * tau))
    If_exact = 1.0 - F_exact

    # Fidelity plot
    plt.figure()
    plt.ylabel("Fidelity to |->")
    plt.title("Fidelity vs k")
    if mode == "IF":
        fid_y = If
        fid_theo = If_exact
        plt.ylabel("Infidelity to |->")
        plt.title("Infidelity vs k")
    elif mode == "F":
        fid_y = F
        fid_theo = F_exact
    else:
        print("Unknown mode in Plot - default: Fidelity")
        fid_y = F
        fid_theo = F_exact
    plt.plot(ks, fid_y, marker="o", label="DB-QITE (numerical)", markersize=2.5)
    #plt.plot(ks, fid_theo, linestyle="--", label="Analytic ITE", linewidth =1.5)
    plt.xlabel("k")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png.replace(".png", "_fidelity.png"), dpi=200)

    # Energy plot
    plt.figure()
    plt.plot(ks, E, marker="o", markersize = 2.5)
    plt.xlabel("k")
    plt.ylabel("Energy <X>")
    plt.title("Energy vs k")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png.replace(".png", "_energy.png"), dpi=200)

    print(f"Saved {out_png.replace('.png','_energy.png')} and {out_png.replace('.png','_fidelity.png')}")
    plt.show()


def plot_multiple_fids(out_png, csv_files, mode="IF" ):
    plt.figure()

    for csv_path in csv_files:
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        ks = data[:,0]
        F  = data[:,2]
        #If = data[:,3]
        If = 1-F

        s_step = extract_s_step(csv_path)
        tau = ks #* s_step
        k = np.linspace(0, int(max(tau)), 12)
        f_exact = 0.5 * (1.0 + np.tanh(2.0 * k*s_step))

        label = f"s={s_step}"
        if mode == "IF":
            fid_y = If
            f_theo = 1 - f_exact
        elif mode == "F":
            fid_y = F
            f_theo = f_exact
        else:
            fid_y = F
            f_theo = f_exact
        plt.plot(tau, fid_y, marker="o", label=label, markersize=2.5)
        plt.plot(k, f_theo, linestyle="--", linewidth=1, label="Analytic_"+f"{s_step}")

    # Analytic curve (use smallest step for smooth reference)
    tau_ref = np.linspace(0, max(tau), 500)
    F_exact = 0.5 * (1.0 + np.tanh(2.0 * tau_ref))
    If_exact = 1.0 - F_exact


    plt.ylabel("Infidelity to |->")
    plt.title("Infidelity vs Steps")
    if mode == "IF":
        fid_theo = If_exact
    elif mode == "F":
        fid_theo = F_exact
        plt.ylabel("Fidelity to |->")
        plt.title("Fidelity vs Time")
    else:
        print("Unknown mode in Plot - default: Fidelity")
        fid_theo = F_exact
        plt.ylabel("Fidelity to |->")
        plt.title("Fidelity vs Time")

    #plt.plot(tau_ref, fid_theo, linestyle="--", linewidth=1, label="Analytic")
    #plt.xlabel("normalized time τ = s k")
    plt.xlabel("Number of steps k")
    plt.yscale("log")
    plt.xscale("log")
    #plt.xlim(-0.1, 8.1)
    #plt.ylim(0.49, 1.01)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)


def plot_optimization(out_png, csv_files):
    plt.figure()

    for csv_path in csv_files:
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        s = data[:,0]
        k = data[:,1]

        epsilon = extract_data_full(csv_path)
        label = f"e={epsilon}"
        plt.plot(s, k, marker="o", label=label, markersize=2.0)

    plt.ylabel("Number of steps k")
    plt.title("Optimization for fixed target epsilon")
    plt.xlabel("Step size s")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

def plot_infidelity_trace(out_png, csv_files):
    plt.figure()

    for csv_path in csv_files:
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        e = data[:,0]
        k = data[:,1]

        s = extract_data_full(csv_path)
        label = f"s={s}"
        plt.plot(e, k, marker="o", label=label, markersize=2.0)

    plt.ylabel("Number of steps k")
    plt.title("#Steps to reach infidelity target for s")
    plt.xlabel("Infidelity target epsilon")
    plt.xscale("log")
    plt.yticks(np.arange(min(k), max(k), 2))
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

# ------------------- Available Plot functions and command dict -----------------------
commands = {
    "single": plot_single_curves,
    "infid_trace": plot_infidelity_trace,
    "multi": plot_multiple_fids,
    "opti": plot_optimization
}
def main():
    command = sys.argv[1]

    if command == "single":
        plot_single_curves(sys.argv[2], sys.argv[3])

    elif command == "multi":
        plot_multiple_fids(sys.argv[2], sys.argv[3:])

    elif command == "opti":
        plot_optimization(sys.argv[2], sys.argv[3:])

    elif command == "infid_trace":
        plot_infidelity_trace(sys.argv[2], sys.argv[3:])

    else:
        print("Unknown command")

if __name__ == "__main__":
    main()