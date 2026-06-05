import sys
import csv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat as un
import FitStuff
from FitStuff import daten_fitten, linear_funct


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

def read_csv(path):
    x, y = [], []

    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            x.append(float(row["g"]))
            y.append(float(row["gap"]))
    return x, y

def read_csv_spectrum(path):
    g = []
    rows_of_energies = []
    gap, gap2, gap3 = [], [], []

    with open(path, newline="") as f:
        reader = csv.reader(f)

        next(reader, None)  # skip header

        for row in reader:
            if not row:
                continue

            gi = float(row[0])
            energies = [float(x) for x in row[1:]]

            g.append(gi)
            rows_of_energies.append(energies)

            # compute gap per row
            if len(energies) >= 2:
                n = len(energies)
                #print(n)
                gap.append(energies[n-2] - energies[n-1])
                gap2.append(energies[n-5] - energies[n-1])
                #gap3.append(energies[n-15] - energies[n-1])
            else:
                gap.append(math.nan)
                gap2.append(math.nan)
                #gap3.append(math.nan)


    # transpose for spectrum plotting
    max_len = max(len(row) for row in rows_of_energies)
    padded = [row + [math.nan]*(max_len - len(row)) for row in rows_of_energies]
    Espec = list(map(list, zip(*padded)))

    return g, Espec, gap, gap2, gap3


def read_csv_spectrumPert(path):
    g = []
    rows_of_energies = []

    with open(path, newline="") as f:
        reader = csv.reader(f)
        next(reader, None)  # skip header

        for row in reader:
            if not row:
                continue
            gi = float(row[0])
            energies = [float(x) for x in row[1:]]
            g.append(gi)
            rows_of_energies.append(energies)

    # transpose for spectrum plotting
    max_len = max(len(row) for row in rows_of_energies)
    padded = [row + [math.nan]*(max_len - len(row)) for row in rows_of_energies]
    Espec = list(map(list, zip(*padded)))

    return g, Espec

def read_csv_gap_analysis(path):
    g, E0, E1, gap, var0, var1 = [], [], [], [], [], []

    with open(path, newline="") as f:
        r = csv.DictReader(f)
        has_var0 = "var0" in (r.fieldnames or [])
        has_var1 = "var1" in (r.fieldnames or [])

        for row in r:
            g.append(float(row["g"]))
            E0.append(float(row["E0"]))
            E1.append(float(row["E1"]))
            gap.append(float(row["gap"]))
            var0.append(float(row["var0"]) if has_var0 and row["var0"] else math.nan)
            var1.append(float(row["var1"]) if has_var1 and row["var1"] else math.nan)

    return g, E0, E1, gap, var0, var1

def read_csv_gapConvergence(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)

        next(reader, None)  # header
        g_vals = [float(x) for x in next(reader)]
        N_row = [int(x) for x in next(reader)]

        # --- ADDED variable-length gap convergence format ---
        # New row 2 format: N_start,Nmax_for_g0,Nmax_for_g1,...
        # Old row 2 format: N0,N1,N2,... shared by every g curve.
        if len(N_row) == len(g_vals) + 1:
            N_start = N_row[0]
            N_by_g = {
                g: list(range(N_start, N_max + 1))
                for g, N_max in zip(g_vals, N_row[1:])
            }
        else:
            N_by_g = {g: list(N_row) for g in g_vals}
        # --- END ADDED variable-length gap convergence format ---

        E0_by_g = {g: [] for g in g_vals}
        E1_by_g = {g: [] for g in g_vals}
        gaps_by_g = {g: [] for g in g_vals}

        for g in g_vals:
            for N in N_by_g[g]:
                row = next(reader)
                E0, E1, gap = map(float, row)
                E0_by_g[g].append(E0)
                E1_by_g[g].append(E1)
                gaps_by_g[g].append(gap)

    return g_vals, N_by_g, E0_by_g, E1_by_g, gaps_by_g

# ----------------------------- Plot Data in different ways --------------
def plot_spectrum(out_png, csv_path, s_path):
    g, Espec, gap, gap2, gap3 = read_csv_spectrum(csv_path)
    gx, sy = read_csv(s_path)

    # --- Spectrum plot ---
    plt.figure()
    plt.xlabel("coupling g")
    plt.ylabel("Energy spectrum")
    plt.title("Spectrum variation with coupling g")

    x = np.linspace(0.4, max(g), 40)
    y_gsq = x**2 - 1/(4*x**2)
    y2_gsq = 2*x**2
    plt.plot(x, y_gsq, linestyle="--", linewidth=1.0, label="g^2")
    plt.plot(x, y2_gsq, linestyle="--", linewidth=1.0, label="2g^2")


    for i, E in enumerate(Espec):
        plt.plot(g, E, marker="o", markersize=2.5, label=f"E{i}")

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png.replace(".png", "_spectrum.png"), dpi=200)

    # --- Gap plot ---
    plt.figure()
    plt.xlabel("coupling g")
    plt.ylabel("Δ = E1 - E0")
    plt.title("Energy gap vs coupling g")
    plt.plot(x, y_gsq, linestyle="--", linewidth=1.0, label="g^2")
    #plt.plot(x, y2_gsq, linestyle="--", linewidth=1.0, label="2g^2")


    plt.plot(g, gap, marker="o", markersize=2.0, label="gap")
    plt.plot(g, gap2, marker="o", markersize=2.0, label="gap2")
    #plt.plot(g, gap3, marker="o", markersize=2.0, label="gap3")

    plt.plot(gx, sy, marker="o", markersize=2.0, label="opt s")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png.replace(".png", "_gap.png"), dpi=200)

    print(f"Saved {out_png.replace('.png','_spectrum.png')} and {out_png.replace('.png','_gap.png')}")


def plot_spectrumPert(out_png, csv_path):
    g, Espec = read_csv_spectrumPert(csv_path)

    # --- Spectrum plot ---
    plt.figure()
    plt.xlabel("coupling g")
    plt.ylabel("Energy spectrum")
    plt.title("Spectrum variation with coupling g")

    x = np.linspace(min(g), max(g), 40)
    y_gsq = x**2
    y2_gsq = 2*x**2
    E1dash = x**2 + 1/(4*x**2)
    E2dash = x**2 - 1/(4*x**2) - 1/(4*x**6)
    E0dash = - 3/(32*x**6)
    plt.plot(x, y_gsq, linestyle="--", linewidth=1.0, label="E=g^2")
    #plt.plot(x, y2_gsq, linestyle="--", linewidth=1.0, label="2g^2")
    plt.plot(x, E1dash, linestyle="--", linewidth=1.0, label="E+")
    plt.plot(x, E2dash, linestyle="--", linewidth=1.0, label="E-")
    plt.plot(x, E0dash, linestyle="--", linewidth=1.0, label="E0")



    for i, E in enumerate(Espec):
        plt.plot(g, E, marker="o", markersize=2.5, label=f"E_{len(Espec)-i-1}")

    plt.grid(True)
    plt.legend()
    plt.ylim(-6, 6)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)



def plot_spectrum_gap(out_png, csv_path):
    g, E0, E1, gap, var0, var1 = read_csv_gap_analysis(csv_path)

    # --- Spectrum plot ---
    plt.figure()
    plt.xlabel("coupling g")
    plt.ylabel("Energy spectrum")
    plt.title("Spectrum variation with coupling g")

    plt.plot(g, E0, marker="o", markersize=2.5, label="E_0")
    plt.plot(g, E1, marker="o", markersize=2.5, label="E_1")
    plt.plot(g, gap, marker="o", markersize=2.5, label="Delta")

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

    has_variance = any(math.isfinite(v) for v in var0 + var1)
    if has_variance:
        variance_png = out_png.replace(".png", "_variance.png")
        plt.figure()
        plt.xlabel("coupling g")
        plt.ylabel("energy variance")
        plt.title("DMRG state variance")
        plt.semilogy(g, var0, marker="o", markersize=2.5, label="var(E_0)")
        plt.semilogy(g, var1, marker="o", markersize=2.5, label="var(E_1)")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(variance_png, dpi=200)
        print(f"Saved {variance_png}")


def plot_gap_convergence(out_png, csv_path):
    g_vals, N_by_g, E0_by_g, E1_by_g, gaps_by_g = read_csv_gapConvergence(csv_path)

    plt.figure()
    for g in g_vals:
        N_vals = N_by_g[g]
        plt.plot(N_vals, gaps_by_g[g], marker="o", label=f"g={g}", markersize=2.0, linewidth=1.0)

    plt.xlabel("N")
    plt.xticks(np.arange(3, 22, 3))

    plt.ylabel("gap")
    plt.title("Gap convergence with system size")
    plt.grid(True)
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)


def plot_gap_convergence_fit(out_png, csv_path):
    g_vals, N_by_g, E0_by_g, E1_by_g, gaps_by_g = read_csv_gapConvergence(csv_path)

    plt.figure()
    print("Finite-size mass-gap fits: gap(N) = m + A/sqrt(N) exp(-N/xi)")
    fit_g_vals = []
    fit_m_vals = []
    fit_dm_vals = []
    fit_xi_vals = []
    fit_dxi_vals = []

    for g in g_vals:
        N_vals = N_by_g[g]
        gaps = gaps_by_g[g]
        line, = plt.plot(N_vals, gaps, marker="o", linestyle="", label=f"g={g} data", markersize=2.5)

        try:
            pars, stdevs, N_fit, gap_fit = FitStuff.fit_mass_gap_convergence(N_vals, gaps, fl = 0.05)
        except Exception as exc:
            print(f"g={g}: fit failed ({exc})")
            continue

        m, A, xi = pars
        dm, dA, dxi = stdevs
        xi_m = xi * m
        d_xi_m = math.sqrt((xi * dm)**2 + (m * dxi)**2)
        fit_g_vals.append(g)
        fit_m_vals.append(m)
        fit_dm_vals.append(dm)
        fit_xi_vals.append(xi)
        fit_dxi_vals.append(dxi)
        plt.plot(N_fit, gap_fit, linestyle="--", color=line.get_color(), label=f"g={g} fit", linewidth=1.0)
        print(
            f"g={g}: "
            f"A={A:.8g} +/- {dA:.2g}, "
            f"m={m:.8g} +/- {dm:.2g}, "
            f"xi={xi:.8g} +/- {dxi:.2g}, "
            f"xi*m={xi_m:.8g} +/- {d_xi_m:.2g}"
        )

    plt.xlabel("N")
    plt.xticks(np.arange(3, 85, 5))
    plt.ylabel("gap")
    plt.title("Gap convergence with finite-size fits")
    plt.grid(True)
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

    if fit_g_vals:
        params_png = out_png.replace(".png", "_fit_params.png")
        g_arr = np.asarray(fit_g_vals, dtype=float)
        g_squared_arr = g_arr**2
        m_arr = np.asarray(fit_m_vals, dtype=float)
        dm_arr = np.asarray(fit_dm_vals, dtype=float)
        xi_arr = np.asarray(fit_xi_vals, dtype=float)
        dxi_arr = np.asarray(fit_dxi_vals, dtype=float)
        g_squared_xi_arr =  xi_arr / g_squared_arr
        g_squared_dxi_arr = dxi_arr / g_squared_arr

        dm_arr = np.where(np.isfinite(dm_arr), dm_arr, np.nan)
        g_squared_dxi_arr = np.where(np.isfinite(g_squared_dxi_arr), g_squared_dxi_arr, np.nan)

        fig, (ax_m, ax_xi) = plt.subplots(2, 1, sharex=True, figsize=(6.4, 7.2))

        ax_m.errorbar(
            g_squared_arr,
            m_arr,
            yerr=dm_arr,
            xerr=None,
            fmt="o-",
            markersize=3.0,
            capsize=3,
        )
        ax_m.set_ylabel("fitted m")
        ax_m.set_title("Finite-size fit parameters vs coupling g^2")
        ax_m.grid(True)

        ax_xi.errorbar(
            g_squared_arr,
            g_squared_xi_arr,
            yerr=g_squared_dxi_arr,
            xerr=None,
            fmt="o-",
            markersize=3.0,
            capsize=3,
        )
        ax_xi.set_xlabel("coupling g^2")
        ax_xi.set_ylabel("g^2 fitted xi")
        ax_xi.grid(True)

        fig.tight_layout()
        fig.savefig(params_png, dpi=200)
        print(f"Saved {params_png}")



def plot_energy_convergence(out_png, csv_path):
    g_vals, N_by_g, E0_by_g, E1_by_g, gaps_by_g = read_csv_gapConvergence(csv_path)

    plt.figure()
    for g in g_vals:
        N_vals = N_by_g[g]
        N_arr = np.asarray(N_vals, dtype=float)
        E0_density = np.asarray(E0_by_g[g], dtype=float) / N_arr
        E1_density = np.asarray(E1_by_g[g], dtype=float) / N_arr

        line, = plt.plot(
            N_vals,
            E0_density,
            marker="o",
            label=f"g={g} E0/N",
            markersize=2.5,
            linewidth=1.0,
        )
        plt.plot(
            N_vals,
            E1_density,
            marker="s",
            linestyle="--",
            color=line.get_color(),
            label=f"g={g} E1/N",
            markersize=2.5,
            linewidth=1.0,
        )

    plt.xlabel("N")
    plt.xticks(np.arange(3, 85, 5))
    plt.ylabel("energy density")
    plt.title("Energy density convergence")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

    raw_out_png = out_png.replace(".png", "_raw.png")
    plt.figure()
    for g in g_vals:
        N_vals = N_by_g[g]
        E0 = np.asarray(E0_by_g[g], dtype=float)
        E1 = np.asarray(E1_by_g[g], dtype=float)

        line, = plt.plot(
            N_vals,
            E0,
            marker="o",
            label=f"g={g} E0",
            markersize=2.5,
            linewidth=1.0,
        )
        plt.plot(
            N_vals,
            E1,
            marker="s",
            linestyle="--",
            color=line.get_color(),
            label=f"g={g} E1",
            markersize=2.5,
            linewidth=1.0,
        )

    plt.xlabel("N")
    plt.xticks(np.arange(3, 71, 3))
    plt.ylabel("energy")
    plt.title("Energy convergence")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(raw_out_png, dpi=200)


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
    plt.xlabel("steps k")
    plt.ylabel("Energy <X>")
    plt.title("Energy vs steps")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png.replace(".png", "_energy.png"), dpi=200)

    print(f"Saved {out_png.replace('.png','_energy.png')} and {out_png.replace('.png','_fidelity.png')}")
    plt.show()


def plot_multiple_fids(out_png, csv_files, mode="IF", theo=False ):
    plt.figure()

    for csv_path in csv_files:
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        ks = data[:,0]
        F  = data[:,2]
        #If = data[:,3]
        If = 1-F

        s_step = extract_s_step(csv_path)
        tau = ks #* s_step
        k = np.linspace(0, int(max(tau)), 100)
        f_exact = 0.5 * (1.0 + np.tanh(k*s_step))
        f_exact = f_exact**2

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
        if theo:
            plt.plot(k, f_theo, linestyle="--", linewidth=1, label="Analytic_"+f"{s_step}")

    # Analytic curve (use smallest step for smooth reference)
    tau_ref = np.linspace(0, max(tau), 500)
    F_exact = 0.5 * (1.0 + np.tanh(2.0 * tau_ref))
    F_exact = F_exact**2
    If_exact = 1.0 - F_exact


    plt.ylabel("Infidelity to Ground-Subspace")
    plt.title("Infidelity vs Steps")
    if mode == "IF":
        fid_theo = If_exact
    elif mode == "F":
        fid_theo = F_exact
        plt.ylabel("Fidelity to Ground-Subspace")
        plt.title("Fidelity vs Time")
    else:
        print("Unknown mode in Plot - default: Fidelity")
        fid_theo = F_exact
        plt.ylabel("Fidelity to Ground-Subspace")
        plt.title("Fidelity vs Time")

    #plt.plot(tau_ref, fid_theo, linestyle="--", linewidth=1, label="Analytic")
    #plt.xlabel("normalized time τ = s k")
    plt.xlabel("Number of steps k")
    plt.yscale("log")
    #plt.xscale("log")
    #plt.xlim(-0.1, 8.1)
    #plt.ylim(0.49, 1.01)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

def plot_multiple_energies(out_png, Elevels_file, csv_files, theo=False ):
    plt.figure()

    for i,csv_path in enumerate(csv_files):
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        ks = data[:,0]
        E  = data[:,1]

        s_step = extract_s_step(csv_path)
        tau = ks #* s_step
        k = np.linspace(0, int(max(tau)), 100)
        label = f"s={s_step}"

        plt.plot(tau, E, marker="o", label=label, markersize=1.5)
        if theo:
            plt.plot(k, E, linestyle="--", linewidth=1, label="Analytic_"+f"{s_step}")

        if i==0:
            data = np.loadtxt(Elevels_file, delimiter=",", skiprows=1)
            Elevels = data[:,0]
            #Plot energy levels
            for Etheo in Elevels:
                plt.hlines(Etheo, min(tau), max(tau), linestyle="--", linewidth=1)



    plt.ylabel("Energy <H>")
    plt.title("Energy vs Steps")


    #plt.xlabel("normalized time τ = s k")
    plt.xlabel("Number of steps k")
    #plt.yscale("log")
    #plt.xscale("log")
    #plt.xlim(-0.1, 8.1)
    #plt.ylim(-5, -3)
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
        plt.plot(s, k, marker="o", label=label, markersize=2.0, linewidth=1.0)

    plt.ylabel("Number of steps k")
    plt.title("Optimization for fixed target infidelity")
    plt.xlabel("Stepsize s")

    plt.legend()
    #plt.ylim(50, 220)
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
        plt.plot(e, k, marker="o", label=label, markersize=1.0)

    plt.ylabel("Number of steps k")
    plt.title("#Steps to reach infidelity target for step size s")
    plt.xlabel("Infidelity target")
    plt.xscale("log")
    plt.yticks(np.arange(min(k), max(k), 5))
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

def plot_infidelity_trace_fit(out_png, csv_files):
    plt.figure()
    st = [3, 3, 3, 4, 3, 3]
    for i,csv_path in enumerate(csv_files):
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        e = data[:,0]
        k = data[:,1]
        elog = np.log(e)

        s = extract_data_full(csv_path)
        label = f"s={s}"
        plt.plot(e, k, marker="o", label=label, markersize=2.0)

        pars, stdevs, xfit, yfit = daten_fitten(linear_funct, elog, k, 30, [30, -3])
        C = np.log(10) / np.log(3)
        m = un(pars[0], stdevs[0])
        A0 = un(pars[1], stdevs[1])
        B_tilde = -C / m
        #st_i = st[i]
        #B_tilde_val = np.round(B_tilde.n,st_i)
        #B_tilde_err = int(np.round(B_tilde.s,st_i)*1e3)
        #if B_tilde_err == 0:
            #B_tilde_err = 1
        plt.plot(np.exp(xfit), yfit, 'o', markersize = 0, label=f"fit: B={B_tilde}", linestyle='--', linewidth=0.5)

        #if B_tilde_err < 1000:
            #plt.plot(np.exp(xfit), yfit, 'o', markersize = 0, label=f"fit: B={B_tilde_val}({B_tilde_err})", linestyle='--', linewidth=0.5)

        print(
                f"fit: of step s={s}: "
                f"m={m} and "
                f"A0={A0}. "
                f"B_tilde={B_tilde}"
        )




    plt.ylabel("Number of steps k")
    plt.title("#Steps to reach infidelity target for stepsize s")
    plt.xlabel("Infidelity target")
    plt.xscale("log")
    plt.yticks(np.arange(5, 60, 5))
    plt.legend()
    #plt.ylim(30, 150)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)

# ------------------- Available Plot functions and command dict -----------------------
commands = {
    "single": plot_single_curves,
    "infid_trace": plot_infidelity_trace,
    "multi": plot_multiple_fids,
    "opti": plot_optimization,
    "infid_fit_trace": plot_infidelity_trace_fit
}
def main():
    command = sys.argv[1]

    if command == "single":
        plot_single_curves(sys.argv[2], sys.argv[3])

    elif command == "multi":
        plot_multiple_fids(sys.argv[2], sys.argv[3:])

    elif command == "multiE":
        plot_multiple_fids(sys.argv[2], sys.argv[5:])
        plot_multiple_energies(sys.argv[3], sys.argv[4], sys.argv[5:])

    elif command == "opti":
        plot_optimization(sys.argv[2], sys.argv[3:])

    elif command == "infid_trace":
        plot_infidelity_trace(sys.argv[2], sys.argv[3:])

    elif command == "infid_fit_trace":
        plot_infidelity_trace_fit(sys.argv[2], sys.argv[3:])

    elif command == "spectrum":
        plot_spectrum(sys.argv[2], sys.argv[3], sys.argv[4] )

    elif command == "spectrumPert":
        plot_spectrumPert(sys.argv[2], sys.argv[3] )

    elif command == "spectrumEgap":
        plot_spectrum_gap(sys.argv[2], sys.argv[3] )

    elif command == "gapConvergence":
        plot_gap_convergence(sys.argv[2], sys.argv[3])

    elif command == "gapConvergenceFit":
        plot_gap_convergence_fit(sys.argv[2], sys.argv[3])

    elif command == "energyConvergence":
        plot_energy_convergence(sys.argv[2], sys.argv[3])

    else:
        print("Unknown command")

if __name__ == "__main__":
    main()
