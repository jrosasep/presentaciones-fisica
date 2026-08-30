#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproduce las dos figuras del doble péndulo usadas en la presentación.

Se integra el sistema exacto de un doble péndulo plano con masas y longitudes
iguales (m1=m2=L1=L2=1) y gravedad g=9.81. Para cada régimen se compara la
trayectoria de control con cinco perturbaciones de theta_2(0).

Condiciones de control documentadas en la presentación:
  región regular:  theta_1(0)=0.25, theta_2(0)=0.35 rad
  región caótica:  theta_1(0)=2.20, theta_2(0)=0.40 rad
  omega_1(0)=omega_2(0)=0 en ambos casos

La observable graficada es la posición horizontal de la segunda masa,
x_2(t)=L1 sin(theta_1)+L2 sin(theta_2).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "images" / "figures"
DATA = ROOT / "data"
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

M1 = 1.0
M2 = 1.0
L1 = 1.0
L2 = 1.0
G = 9.81
T = np.linspace(0.0, 25.0, 2501)
PERTURBATIONS = (0.0, 1e-8, 1e-6, 1e-4, 1e-2, 2.0)

COLORS = ("#223C6A", "#6F91C5", "#2E7D5B", "#B77A19", "#B24C4C", "#687080")
LINESTYLES = ("-", "--", ":", "--", "--", ":")
LINEWIDTHS = (2.6, 1.8, 1.8, 1.8, 1.8, 1.5)
LABELS = (
    "control",
    r"$\delta\theta_2(0)=10^{-8}$ rad",
    r"$\delta\theta_2(0)=10^{-6}$ rad",
    r"$\delta\theta_2(0)=10^{-4}$ rad",
    r"$\delta\theta_2(0)=10^{-2}$ rad",
    r"$\delta\theta_2(0)=2$ rad (referencia lejana)",
)


def rhs(_t: float, state: np.ndarray) -> np.ndarray:
    """Ecuaciones de movimiento en las variables (theta1, omega1, theta2, omega2)."""
    theta1, omega1, theta2, omega2 = state
    delta = theta1 - theta2
    denominator = 2.0 * M1 + M2 - M2 * np.cos(2.0 * delta)

    domega1 = (
        -G * (2.0 * M1 + M2) * np.sin(theta1)
        - M2 * G * np.sin(theta1 - 2.0 * theta2)
        - 2.0
        * M2
        * np.sin(delta)
        * (omega2**2 * L2 + omega1**2 * L1 * np.cos(delta))
    ) / (L1 * denominator)

    domega2 = (
        2.0
        * np.sin(delta)
        * (
            omega1**2 * L1 * (M1 + M2)
            + G * (M1 + M2) * np.cos(theta1)
            + omega2**2 * L2 * M2 * np.cos(delta)
        )
    ) / (L2 * denominator)

    return np.array((omega1, domega1, omega2, domega2))


def integrate(theta1_0: float, theta2_0: float) -> np.ndarray:
    state0 = np.array((theta1_0, 0.0, theta2_0, 0.0), dtype=float)
    solution = solve_ivp(
        rhs,
        (T[0], T[-1]),
        state0,
        t_eval=T,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    theta1 = solution.y[0]
    theta2 = solution.y[2]
    return L1 * np.sin(theta1) + L2 * np.sin(theta2)


def make_figure(stem: str, theta1_0: float, theta2_0: float, title: str) -> None:
    trajectories = np.column_stack(
        [integrate(theta1_0, theta2_0 + delta) for delta in PERTURBATIONS]
    )

    header = "t_s," + ",".join(
        ["x2_control_m"]
        + [f"x2_delta_theta2_{delta:g}_m" for delta in PERTURBATIONS[1:]]
    )
    np.savetxt(
        DATA / f"{stem}.csv",
        np.column_stack((T, trajectories)),
        delimiter=",",
        header=header,
        comments="",
    )

    fig = plt.figure(figsize=(10.0, 4.5), dpi=300, facecolor="white")
    ax = fig.add_axes((0.075, 0.17, 0.615, 0.74))
    for index, label in enumerate(LABELS):
        ax.plot(
            T,
            trajectories[:, index],
            color=COLORS[index],
            linestyle=LINESTYLES[index],
            linewidth=LINEWIDTHS[index],
            label=label,
            alpha=0.98 if index < 5 else 0.82,
        )

    ax.set_xlim(0.0, 25.0)
    ax.set_ylim(-2.1, 2.1)
    ax.set_xticks(np.arange(0.0, 25.1, 5.0))
    ax.set_yticks(np.arange(-2.0, 2.1, 1.0))
    ax.set_xlabel(r"tiempo $t$ [s]", fontsize=11)
    ax.set_ylabel(r"posición horizontal $x_2(t)$ [m]", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold", color="#223C6A", pad=8)
    ax.grid(True, color="#D9DEE7", linewidth=0.8, alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#2E2E2E")
    ax.tick_params(colors="#2E2E2E", labelsize=9)

    legend = ax.legend(
        title=r"Perturbación de $\theta_2(0)$",
        loc="center left",
        bbox_to_anchor=(1.04, 0.5),
        frameon=True,
        borderpad=0.9,
        labelspacing=1.0,
        fontsize=7.8,
        title_fontsize=8.8,
    )
    legend.get_frame().set_edgecolor("#D9DEE7")
    legend.get_frame().set_linewidth(0.9)

    fig.savefig(OUT / f"{stem}.png", dpi=300, facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", facecolor="white")
    plt.close(fig)


def main() -> None:
    make_figure(
        "doble_pendulo_region_regular",
        0.25,
        0.35,
        "Trayectorias cercanas en una región regular",
    )
    make_figure(
        "doble_pendulo_region_caotica",
        2.20,
        0.40,
        "Separación de trayectorias en una región caótica",
    )
    print(f"Figuras guardadas en {OUT}")
    print(f"Datos guardados en {DATA}")


if __name__ == "__main__":
    main()
