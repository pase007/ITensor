import sys
import csv
import os
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

def read_csv_correlation_profile(path):
    rows = []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            clean = {}
            for key, value in row.items():
                if key in {"N", "i0", "comp_i", "comp_j", "r", "rMin", "rMax"}:
                    clean[key] = int(value)
                else:
                    clean[key] = float(value)
            if "m_corr" not in clean:
                xi = clean.get("xi_corr", math.nan)
                clean["m_corr"] = 1.0 / xi if math.isfinite(xi) and xi != 0.0 else math.nan
            rows.append(clean)
    return rows

def read_mass_by_g_squared(path, mass_column):
    masses = {}

    with open(path, newline="") as f:
        r = csv.DictReader(f)
        if r.fieldnames is None:
            raise ValueError(f"{path} has no header")
        if "g" not in r.fieldnames:
            raise ValueError(f"{path} has no 'g' column")
        if mass_column not in r.fieldnames:
            raise ValueError(f"{path} has no '{mass_column}' column")

        for row in r:
            try:
                g = float(row["g"])
                mass = float(row[mass_column])
            except (TypeError, ValueError):
                continue

            if not math.isfinite(g) or not math.isfinite(mass):
                continue

            # Correlation profile files contain many rows per g, one per distance r.
            # The mass columns are constant for those rows, so keep the first value.
            g_squared_key = round(g * g, 12)
            masses.setdefault(g_squared_key, (g, mass))

    return masses

def read_mass_difference(left_path, right_path, left_column, right_column):
    left = read_mass_by_g_squared(left_path, left_column)
    right = read_mass_by_g_squared(right_path, right_column)
    common_g_squared = sorted(set(left).intersection(right))

    rows = []
    for g_squared in common_g_squared:
        g_left, mass_left = left[g_squared]
        g_right, mass_right = right[g_squared]
        rows.append({
            "g": 0.5 * (g_left + g_right),
            "g_squared": g_squared,
            left_column: mass_left,
            right_column: mass_right,
            "difference": mass_left - mass_right,
        })

    return rows

def merge_correlation_files(out_path, input_paths):
    if not input_paths:
        raise ValueError("mergeCorr needs at least one input file")

    header = None
    rows_by_key = {}
    duplicate_count = 0
    conflict_count = 0

    for path in input_paths:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise ValueError(f"{path} has no header")
            if header is None:
                header = reader.fieldnames
                missing = {"g", "N", "comp_i", "comp_j", "r"} - set(header)
                if missing:
                    raise ValueError(f"{path} is missing merge-key columns: {sorted(missing)}")
            elif reader.fieldnames != header:
                raise ValueError(f"{path} has a different header")

            for row in reader:
                try:
                    key = (
                        round(float(row["g"]), 12),
                        int(row["N"]),
                        int(row["comp_i"]),
                        int(row["comp_j"]),
                        int(row["r"]),
                    )
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"{path} has an invalid merge key row: {row}") from exc

                if key in rows_by_key:
                    duplicate_count += 1
                    if row != rows_by_key[key]:
                        conflict_count += 1
                        print(f"mergeCorr warning: duplicate key with different row, keeping first: {key}")
                    continue

                rows_by_key[key] = row

    sorted_rows = [
        rows_by_key[key]
        for key in sorted(rows_by_key)
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(sorted_rows)

    print(f"mergeCorr wrote {len(sorted_rows)} rows to {out_path}")
    print(f"mergeCorr skipped {duplicate_count} duplicate rows ({conflict_count} conflicting)")

# ----------------------------- Plot Data in different ways --------------
def plot_spectrum(out_png, csv_path, s_path):
    g, Espec, gap, gap2, gap3 = read_csv_spectrum(csv_path)
    gx, sy = read_csv(s_path)

    # --- Spectrum plot ---
    # Plot raw connected correlations C(r) from CSV columns r and corr for each g.
    # Plot raw connected correlations C(r) from CSV columns r and corr for each g.
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
    # Plot raw decay |C(r)| from CSV columns r and abs_corr for each g on a semilog y-axis.
    # Plot raw decay |C(r)| from CSV columns r and abs_corr for each g on a semilog y-axis.
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
        # Plot xi_fit * g^2 vs 1/g^2, where xi_fit = -1/slope from the local log-linear fit above.
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
        params_xi_gsq_png = out_png.replace(".png", "_fit_params_xi_gsq.png")
        fit_params_csv = csv_path.rsplit(".", 1)[0] + "_fit_params.csv"
        g_arr = np.asarray(fit_g_vals, dtype=float)
        g_squared_arr = g_arr**2
        m_arr = np.asarray(fit_m_vals, dtype=float)
        dm_arr = np.asarray(fit_dm_vals, dtype=float)
        xi_arr = np.asarray(fit_xi_vals, dtype=float)
        dxi_arr = np.asarray(fit_dxi_vals, dtype=float)

        dm_arr = np.where(np.isfinite(dm_arr), dm_arr, np.nan)

        with open(fit_params_csv, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["g", "g_squared", "m", "dm", "xi", "dxi"],
            )
            writer.writeheader()
            for g, g_squared, m, dm, xi, dxi in zip(
                g_arr,
                g_squared_arr,
                m_arr,
                dm_arr,
                xi_arr,
                dxi_arr,
            ):
                writer.writerow({
                    "g": g,
                    "g_squared": g_squared,
                    "m": m,
                    "dm": dm,
                    "xi": xi,
                    "dxi": dxi,
                })
        print(f"Saved {fit_params_csv}")

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
            xi_arr,
            yerr=dxi_arr,
            xerr=None,
            fmt="o-",
            markersize=3.0,
            capsize=3,
        )
        ax_xi.set_xlabel("coupling g^2")
        ax_xi.set_ylabel("fitted xi")
        ax_xi.grid(True)

        fig.tight_layout()
        fig.savefig(params_xi_gsq_png, dpi=200)
        print(f"Saved {params_xi_gsq_png}")



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

def plot_mass_difference(
    out_png,
    left_path,
    right_path,
    left_column="m_corr",
    right_column="gap",
    left_label=None,
    right_label=None,
):
    rows = read_mass_difference(left_path, right_path, left_column, right_column)
    if not rows:
        raise ValueError(
            f"No common finite g^2 values found between {left_path}:{left_column} "
            f"and {right_path}:{right_column}"
        )

    left_label = left_label or left_column
    right_label = right_label or right_column

    g_squared = [row["g_squared"] for row in rows]
    left_mass = [row[left_column] for row in rows]
    right_mass = [row[right_column] for row in rows]
    difference = [row["difference"] for row in rows]

    csv_name = os.path.basename(out_png).replace(".png", "_matched_difference.csv")
    csv_out = os.path.join(os.path.dirname(left_path) or ".", csv_name)
    with open(csv_out, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["g", "g_squared", left_column, right_column, "difference"],
        )
        writer.writeheader()
        writer.writerows(rows)

    fig, (ax_mass, ax_diff) = plt.subplots(2, 1, sharex=True, figsize=(6.4, 7.2))

    ax_mass.plot(g_squared, left_mass, marker="o", markersize=3.0, linewidth=1.0, label=left_label)
    ax_mass.plot(g_squared, right_mass, marker="s", markersize=3.0, linewidth=1.0, label=right_label)
    ax_mass.set_ylabel("mass")
    ax_mass.set_title("Matched masses")
    ax_mass.grid(True)
    ax_mass.legend(loc="best")

    ax_diff.axhline(0.0, color="black", linewidth=0.8)
    ax_diff.plot(
        g_squared,
        difference,
        marker="o",
        markersize=3.0,
        linewidth=1.0,
        label=f"{left_label} - {right_label}",
    )
    ax_diff.set_xlabel("coupling g^2")
    ax_diff.set_ylabel("mass difference")
    ax_diff.set_title("Mass difference")
    ax_diff.grid(True)
    ax_diff.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    print(f"Matched {len(rows)} common g^2 values")
    print(f"{'g':>10} {'g^2':>12} {left_label:>16} {right_label:>16} {'diff':>16}")
    for row in rows:
        print(
            f"{row['g']:10.6g} "
            f"{row['g_squared']:12.6g} "
            f"{row[left_column]:16.8g} "
            f"{row[right_column]:16.8g} "
            f"{row['difference']:16.8g}"
        )
    print(f"Saved {out_png} and {csv_out}")


def plot_correlation_profile(out_png, csv_path):
    rows = read_csv_correlation_profile(csv_path)
    if not rows:
        raise ValueError(f"No correlation data in {csv_path}")

    by_g = {}
    for row in rows:
        by_g.setdefault(row["g"], []).append(row)

    profile_png = out_png.replace(".png", "_profile.png")
    log_profile_png = out_png.replace(".png", "_log_profile.png")
    log_profile_fit_png = out_png.replace(".png", "_log_profile_fit_every_second.png")
    fit_xi_gsq_inv_gsq_png = out_png.replace(".png", "_fit_xi_times_gsq_vs_inv_gsq.png")
    fit_mass_gsq_png = out_png.replace(".png", "_fit_mass_vs_gsq.png")
    xi_png = out_png.replace(".png", "_xi_vs_g.png")
    mass_png = out_png.replace(".png", "_mass_vs_g.png")
    xi_gsq_png = out_png.replace(".png", "_xi_times_gsq_vs_g.png")
    xi_gsq_inv_gsq_png = out_png.replace(".png", "_xi_times_gsq_vs_inv_gsq.png")
    xi_gsq_inv_gsq_logy_png = out_png.replace(".png", "_xi_times_gsq_vs_inv_gsq_logy.png")

    # Plot raw connected correlations C(r) from CSV columns r and corr for each g.
    plt.figure()
    for g, grows in sorted(by_g.items()):
        grows = sorted(grows, key=lambda x: x["r"])
        r = [x["r"] for x in grows]
        corr = [x["corr"] for x in grows]
        plt.plot(r, corr, marker="o", markersize=2.5, linewidth=1.0, label=f"g={g:.3g}")
    plt.xlabel("distance r")
    plt.ylabel("connected correlation C(r)")
    plt.title("Ground-state connected correlations")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(profile_png, dpi=200)

    # Plot raw decay |C(r)| from CSV columns r and abs_corr for each g on a semilog y-axis.
    plt.figure()
    for g, grows in sorted(by_g.items()):
        grows = sorted(grows, key=lambda x: x["r"])
        r = [x["r"] for x in grows if math.isfinite(x["abs_corr"]) and x["abs_corr"] > 0.0]
        abs_corr = [x["abs_corr"] for x in grows if math.isfinite(x["abs_corr"]) and x["abs_corr"] > 0.0]
        if not r:
            continue
        plt.semilogy(r, abs_corr, marker="o", markersize=2.5, linewidth=1.0, label=f"g={g:.3g}")
    plt.xlabel("distance r")
    plt.ylabel("|C(r)|")
    plt.title("Correlation decay")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(log_profile_png, dpi=200)

    # Plot every second g curve as log(|C(r)|); fit data is computed here from CSV abs_corr via log(|C(r)|) = a + b*r.
    corr_fit_floor = 1e-16
    plt.figure()
    fit_g_vals = []
    fit_xi_vals = []
    fit_mass_vals = []
    for curve_idx, (g, grows) in enumerate(sorted(by_g.items())):
        if curve_idx % 2 != 0:
            continue

        grows = sorted(grows, key=lambda x: x["r"])
        fit_points = [
            (x["r"], x["abs_corr"])
            for x in grows
            if math.isfinite(x["abs_corr"]) and x["abs_corr"] >= corr_fit_floor
        ]
        if len(fit_points) < 2:
            continue

        r = np.array([x[0] for x in fit_points], dtype=float)
        abs_corr = np.array([x[1] for x in fit_points], dtype=float)
        log_abs_corr = np.log(abs_corr)
        slope, intercept = np.polyfit(r, log_abs_corr, 1)
        fit_log_abs_corr = intercept + slope * r
        xi_fit = -1.0 / slope if slope < 0.0 else math.nan
        mass_fit = -slope if slope < 0.0 else math.nan
        if math.isfinite(xi_fit) and math.isfinite(mass_fit):
            fit_g_vals.append(g)
            fit_xi_vals.append(xi_fit)
            fit_mass_vals.append(mass_fit)

        points = plt.scatter(r, log_abs_corr, s=6, label=f"g={g:.3g}")
        plt.plot(
            r,
            fit_log_abs_corr,
            linestyle="--",
            linewidth=1.0,
            color=points.get_facecolor()[0],
            label=f"fit g={g:.3g}, xi={xi_fit:.3g}, m={mass_fit:.3g}",
        )
    plt.xlabel("distance r")
    plt.ylabel("log(|C(r)|)")
    plt.title("Log correlation decay with linear fits")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(log_profile_fit_png, dpi=200)

    fit_inv_g_squared_vals = [
        1.0 / (g * g) if g != 0.0 else math.nan
        for g in fit_g_vals
    ]
    fit_g_squared_vals = [g * g for g in fit_g_vals]
    fit_xi_times_gsq = [
        xi * g * g if math.isfinite(xi) else math.nan
        for g, xi in zip(fit_g_vals, fit_xi_vals)
    ]

    if fit_g_vals:
        # Plot xi_fit*g^2 vs 1/g^2, where xi_fit = -1/slope from the local log-linear fit above.
        plt.figure()
        plt.plot(fit_inv_g_squared_vals, fit_xi_times_gsq, marker="o", markersize=3.0, linewidth=1.0, label="xi_fit * g^2")
        plt.xlabel("inverse coupling 1/g^2")
        plt.ylabel("xi_fit * g^2")
        plt.title("Slope-fit scaled correlation length")
        plt.grid(True)
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(fit_xi_gsq_inv_gsq_png, dpi=200)

        # Plot m_fit vs g^2, where m_fit = -slope from the local log-linear fit above.
        plt.figure()
        plt.plot(fit_g_squared_vals, fit_mass_vals, marker="o", markersize=3.0, linewidth=1.0, label="m_fit")
        plt.xlabel("coupling g^2")
        plt.ylabel("m_fit")
        plt.title("Slope-fit correlation mass")
        plt.grid(True)
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(fit_mass_gsq_png, dpi=200)

    g_vals = []
    xi_corr = []
    xi_gap = []
    for g, grows in sorted(by_g.items()):
        first = grows[0]
        g_vals.append(g)
        xi_corr.append(first["xi_corr"])
        xi_gap.append(first["xi_gap"])

    # Plot correlation lengths from CSV columns xi_corr and xi_gap = 1/gap against g^2.
    plt.figure()
    g_squared_vals = [g * g for g in g_vals]
    plt.plot(g_squared_vals, xi_corr, marker="o", markersize=3.0, linewidth=1.0, label="fit from |C(r)|")
    plt.plot(g_squared_vals, xi_gap, marker="s", markersize=3.0, linewidth=1.0, label="1/gap")
    plt.xlabel("coupling g^2")
    plt.ylabel("correlation length")
    plt.title("Correlation length scan")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(xi_png, dpi=200)

    m_corr = []
    gap = []
    for g, grows in sorted(by_g.items()):
        first = grows[0]
        m_corr.append(first["m_corr"])
        gap.append(first["gap"])

    # Plot masses from CSV columns m_corr and gap against g^2.
    plt.figure()
    plt.plot(g_squared_vals, m_corr, marker="o", markersize=3.0, linewidth=1.0, label="1/xi_corr")
    if any(math.isfinite(x) for x in gap):
        plt.plot(g_squared_vals, gap, marker="s", markersize=3.0, linewidth=1.0, label="DMRG gap")
    plt.xlabel("coupling g^2")
    plt.ylabel("mass")
    plt.title("Correlation mass scan")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(mass_png, dpi=200)

    xi_times_gsq = [
        xi * g * g if math.isfinite(xi) else math.nan
        for g, xi in zip(g_vals, xi_corr)
    ]

    # Plot scaled CSV correlation lengths xi_corr*g^2 against g^2.
    plt.figure()
    plt.plot(g_squared_vals, xi_times_gsq, marker="o", markersize=3.0, linewidth=1.0, label="xi_corr * g^2")
    plt.xlabel("coupling g^2")
    plt.ylabel("xi_corr * g^2")
    plt.title("Scaled correlation length")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(xi_gsq_png, dpi=200)

    inv_g_squared_vals = [
        1.0 / (g * g) if g != 0.0 else math.nan
        for g in g_vals
    ]

    # Plot scaled CSV correlation lengths xi_corr*g^2 against 1/g^2.
    plt.figure()
    plt.plot(inv_g_squared_vals, xi_times_gsq, marker="o", markersize=3.0, linewidth=1.0, label="xi_corr * g^2")
    plt.xlabel("inverse coupling 1/g^2")
    plt.ylabel("xi * g^2")
    plt.title("Scaled correlation length vs inverse coupling")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(xi_gsq_inv_gsq_png, dpi=200)

    print(f"Saved {profile_png}, {log_profile_png}, {log_profile_fit_png}, {fit_xi_gsq_inv_gsq_png}, {fit_mass_gsq_png}, {xi_png}, {mass_png}, {xi_gsq_png}, {xi_gsq_inv_gsq_png}, and {xi_gsq_inv_gsq_logy_png}")
    # Plot scaled CSV correlation lengths xi_corr*g^2 against 1/g^2 with logarithmic y-axis.
    plt.figure()
    plt.plot(inv_g_squared_vals, xi_times_gsq, marker="o", markersize=3.0, linewidth=1.0, label="xi_corr * g^2")
    plt.xlabel("inverse coupling 1/g^2")
    plt.ylabel("xi * g^2")
    plt.yscale("log")
    plt.title("Scaled correlation length vs inverse coupling")
    plt.grid(True)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(xi_gsq_inv_gsq_logy_png, dpi=200)


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

    elif command == "massDiff":
        left_column = sys.argv[5] if len(sys.argv) > 5 else "m_corr"
        right_column = sys.argv[6] if len(sys.argv) > 6 else "gap"
        plot_mass_difference(sys.argv[2], sys.argv[3], sys.argv[4], left_column, right_column)

    elif command == "mergeCorr":
        merge_correlation_files(sys.argv[2], sys.argv[3:])

    elif command == "corrProfile":
        plot_correlation_profile(sys.argv[2], sys.argv[3])

    else:
        print("Unknown command")

if __name__ == "__main__":
    main()
