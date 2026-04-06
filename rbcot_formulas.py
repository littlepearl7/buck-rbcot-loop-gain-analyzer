# -*- coding: utf-8 -*-
"""
Lu et al., IEEE Trans. Power Electron., 2023 -- RBCOT 环路增益 Pade 简化式与近似设计式。
公式编号与论文一致：(12)(13) ESR；(24)-(27) 固定斜坡；(38)-(44) SW 滤波器。

修复记录：
1. esr_pade_loop: omega_c_approx 改为 sqrt(A0)*omega_p1（论文式13，功率级极点）
2. fixed_slope_pade_loop: Q_ps2 改为 sqrt(inner*(Ax+K))/(pi*den_q)（开根号）
3. sw_filter_pade_loop:
   - delta = LCo/tau - Ton/2（量纲 s，不除以 Ts）
   - A0r = LCo/(LCo+Axr)（<1，经零点提升后穿越0dB，与论文Fig.24吻合）
   - omega_pr2 = sqrt(LCo/(tau*Bxr))，Q_pr2 = sqrt(LCo*Bxr/(tau*Axr^2))
"""
from __future__ import annotations

import numpy as np
import control as ct


def esr_approx_loop(L, Cout, Rc, Rload, Ts, Ton):
    """
    ESR ripple compensation simplified loop gain model from the paper.

    This approximation keeps the ESR zero and the low-frequency power-stage
    double pole, while neglecting the switching-related high-frequency double
    pole p2 when generating the Bode plot.
    """
    s = ct.TransferFunction.s

    omega_z = 1.0 / (Rc * Cout)
    omega_p1 = 1.0 / np.sqrt(L * Cout * (1.0 + Rc / Rload))
    Q_p1 = 1.0 / (omega_p1 * (Cout * Rc + L / Rload))

    A0_num = L * Cout * (1.0 + Rc / Rload)
    term_1 = (Ts / np.pi) ** 2
    term_2 = (Ton / np.pi) ** 2
    term_3 = (Ton / 2.0) * (Rc * Cout - Ton / 2.0)
    A0_den = term_1 + term_2 + term_3
    if A0_den <= 0:
        A0_den = max(1e-30, term_1 + term_2)
    A0 = A0_num / A0_den

    G_z = 1.0 + s / omega_z
    G_p1 = 1.0 / (1.0 + s / (Q_p1 * omega_p1) + (s / omega_p1) ** 2)
    Floop = A0 * G_z * G_p1

    omega_c_approx = np.sqrt(max(A0, 0.0)) * omega_p1
    PM_approx_rad = np.arctan(omega_c_approx / omega_z)

    meta = dict(
        A0=A0,
        omega_z=omega_z,
        omega_p1=omega_p1,
        Q_p1=Q_p1,
        omega_c_approx=omega_c_approx,
        PM_approx_rad=PM_approx_rad,
    )
    return Floop, meta


def esr_approx_bode(L, Cout, Rc, Rload, Ts, Ton, w):
    """
    Approximate ESR model Bode data that matches the paper's design equations.

    Magnitude ignores the ESR zero's gain boost around crossover as assumed in
    Eq. (13), while phase includes the zero contribution as in Eq. (14).
    """
    _, meta = esr_approx_loop(L, Cout, Rc, Rload, Ts, Ton)

    omega_p1 = meta["omega_p1"]
    Q_p1 = meta["Q_p1"]
    omega_z = meta["omega_z"]
    A0 = meta["A0"]

    den_mag = np.sqrt((1.0 - (w / omega_p1) ** 2) ** 2 + (w / (Q_p1 * omega_p1)) ** 2)
    mag = A0 / np.maximum(den_mag, 1e-30)

    phase_p1 = -np.degrees(np.arctan2(w / (Q_p1 * omega_p1), 1.0 - (w / omega_p1) ** 2))
    phase_z = np.degrees(np.arctan(w / omega_z))
    phase_deg = np.unwrap(np.radians(phase_p1 + phase_z)) * 180.0 / np.pi
    freq = w / (2.0 * np.pi)

    return freq, mag, phase_deg, meta


def esr_pade_loop(L, Cout, Rc, Rload, Ts, Ton):
    """
    式(12)(13)：Floop(s) = A0(1+s/ωz) / [(二阶 p1)(二阶 p2)]
    A0 分母：(Tsw/π)²+(Ton/π)² + (Ton/2)·Rc·Co − (Ton/2)²
    
    修复：交换 p1 和 p2 的定义
    - p1 = 采样相关极点（高频）
    - p2 = 功率级极点（低频）
    """
    s = ct.TransferFunction.s
    omega_z = 1.0 / (Rc * Cout)
    
    # 功率级极点（低频）
    omega_p1 = 1.0 / np.sqrt(L * Cout * (1.0 + Rc / Rload))
    Q_p1 = 1.0 / (omega_p1 * (Cout * Rc + L / Rload))

    # 采样相关极点（高频）
    term = 1.0 + np.pi**2 * (Rc * Cout / (2.0 * Ton) - 0.25) + (Ts / Ton) ** 2
    omega_p2 = np.pi * np.sqrt(term) / Ts
    denom_q2 = Rc * Cout / Ts - Ton / (2.0 * Ts) + Ts / (2.0 * Ton)
    if abs(denom_q2) < 1e-18:
        denom_q2 = 1e-18 * np.sign(denom_q2 if denom_q2 != 0 else 1.0)
    Q_p2 = np.sqrt(term) / (np.pi * denom_q2)
    Q_p2 = float(np.clip(Q_p2, 0.05, 50.0))

    A0_num = L * Cout * (1.0 + Rc / Rload)
    
    # 改进：按论文式(12)重新组织计算，避免小项被淹没
    term_1 = (Ts / np.pi) ** 2
    term_2 = (Ton / np.pi) ** 2
    term_3 = (Ton / 2.0) * (Rc * Cout - Ton / 2.0)
    A0_den = term_1 + term_2 + term_3
    
    if A0_den <= 0:
        A0_den = max(1e-30, term_1 + term_2)
    A0 = A0_num / A0_den

    G_z = 1.0 + s / omega_z
    G_p1 = 1.0 / (1.0 + s / (Q_p1 * omega_p1) + (s / omega_p1) ** 2)
    G_p2 = 1.0 / (1.0 + s / (Q_p2 * omega_p2) + (s / omega_p2) ** 2)
    Floop = A0 * G_z * G_p1 * G_p2

    # 论文式(13): omega_c ≈ sqrt(A0) * omega_p1（功率级极点，低频）
    # PM ≈ arctan(omega_c / omega_z)，式(14)
    omega_c_approx = np.sqrt(max(A0, 0.0)) * omega_p1
    PM_approx_rad = np.arctan(omega_c_approx / omega_z)

    meta = dict(
        A0=A0, omega_z=omega_z, omega_p1=omega_p1, Q_p1=Q_p1,
        omega_p2=omega_p2, Q_p2=Q_p2,
        omega_c_approx=omega_c_approx,
        PM_approx_rad=PM_approx_rad,
    )
    return Floop, meta


def fixed_slope_pade_loop(L, Cout, Rload, Vout, Ts, Ton, Vslope):
    """
    式(24)-(27)。Sr = Vslope/Tsw，K = Sr·Tsw·Co·L/Vout = Vslope·Co·L/Vout。
    ωps2、Qps2 为论文完整式（Qps2 分母含 π·(Tsw·Ax/2 + Ton·K/2)）。
    """
    s = ct.TransferFunction.s
    Sr = Vslope / Ts
    K = Sr * Ts * Cout * L / Vout

    Ax = (Ts / np.pi) ** 2 + (Ton / np.pi) ** 2 - (Ton / 2.0) ** 2
    Ax = max(Ax, 1e-30)
    A0_s = L * Cout / (Ax + K)

    omega_xs = np.pi / Ts
    Q_xs = 2.0 / np.pi
    omega_ps1 = 1.0 / np.sqrt(L * Cout)
    Q_ps1 = Rload * np.sqrt(Cout / L)

    # 论文式(24)-(27)：omega_ps2 和 Q_ps2
    # inner = Ts^2*Ax + Ton^2*K  (量纲 s^4，用于 omega_ps2 的根号内)
    # Q_ps2 = sqrt(inner*(Ax+K)) / [pi*(Ts/2*Ax+Ton/2*K)]
    inner_ps2 = Ts**2 * Ax + Ton**2 * K
    inner_ps2 = max(inner_ps2, 1e-30)
    omega_ps2 = np.pi * np.sqrt((Ax + K) / inner_ps2)
    den_q = (Ts / 2.0) * Ax + (Ton / 2.0) * K
    if abs(den_q) < 1e-30:
        den_q = 1e-30
    # 修复：Q_ps2 = sqrt(inner*(Ax+K)) / (pi*den_q)，而不是 inner*(Ax+K)/(pi*den_q)^2
    Q_ps2 = np.sqrt(inner_ps2 * (Ax + K)) / (np.pi * den_q)
    Q_ps2 = float(np.clip(Q_ps2, 0.05, 50.0))

    num = A0_s * (1.0 + s / (Q_xs * omega_xs) + (s / omega_xs) ** 2)
    den_ps2 = 1.0 + s / (Q_ps2 * omega_ps2) + (s / omega_ps2) ** 2
    den_ps1 = 1.0 + s / (Q_ps1 * omega_ps1) + (s / omega_ps1) ** 2
    Floop = num / (den_ps2 * den_ps1)

    omega_c_s_approx = np.sqrt(max(A0_s, 0)) * omega_ps1
    meta = dict(
        Ax=Ax, K=K, A0_s=A0_s, omega_xs=omega_xs, Q_xs=Q_xs,
        omega_ps1=omega_ps1, Q_ps1=Q_ps1, omega_ps2=omega_ps2, Q_ps2=Q_ps2,
        omega_c_s_approx=omega_c_s_approx,
        Vslope_best=2.3 * Ts**2 * Vout / (np.pi**2 * L * Cout),
    )
    return Floop, meta


def sw_filter_pade_loop(L, Cout, Rload, Ts, Ton, Rx, Cx):
    """
    论文式(38)-(44)：SW 滤波器 IL AC 斜坡补偿环路增益。

    符号与论文一致：
      tau   = Rx*Cx
      delta = LCo/tau - Ton/2          [单位 s]
      Axr   = (Ts/pi)^2 + (Ton/pi)^2 + delta*(Ton/2)   [s^2]
      Bxr   = (Ton/pi)^2*delta + Ton*Ts^2/(2*pi^2)      [s^3]
      A0r   = LCo / (LCo + Axr)        [无量纲，< 1]
      omega_zr = 1/tau
      omega_pr1 = 1/sqrt(LCo),  Q_pr1 = RL*sqrt(Co/L)
      omega_pr2 = sqrt(LCo/(tau*Bxr)),  Q_pr2 = sqrt(LCo*Bxr/(tau*Axr^2))
    """
    s = ct.TransferFunction.s
    tau = Rx * Cx
    lc  = L * Cout

    # delta [s]: 对应论文中 LCo/tau - Ton/2
    delta = lc / tau - Ton / 2.0
    # delta 必须 > 0 才能在 Axr 里给出正贡献；物理上 tau << 2*LCo/Ton 时才有意义
    delta = float(max(delta, 1e-30))

    # Axr [s^2]
    Ax_r = (Ts / np.pi) ** 2 + (Ton / np.pi) ** 2 + delta * (Ton / 2.0)
    Ax_r = float(max(Ax_r, 1e-30))

    # Bxr [s^3]
    Bx_r = (Ton / np.pi) ** 2 * delta + Ton * (Ts ** 2) / (2.0 * np.pi ** 2)
    Bx_r = float(max(Bx_r, 1e-30))

    # A0r = LCo / (LCo + Axr)  论文式(39)，验证与Fig.24一致
    # A0r < 1，零点 omega_zr << omega_pr1，增益从 A0r 出发经零点提升后穿越 0dB
    A0_r = lc / (lc + Ax_r)

    # ωzr = 1/tau  (SW滤波器零点)
    omega_zr = 1.0 / tau

    # 功率级极点
    omega_pr1 = 1.0 / np.sqrt(lc)
    Q_pr1 = Rload * np.sqrt(Cout / L)

    # 采样相关极点
    # ωpr2 = sqrt(LCo / (tau * Bxr))  [sqrt(s^2/(s*s^3)) = sqrt(1/s^2) = 1/s ✓]
    omega_pr2 = np.sqrt(lc / (tau * Bx_r))
    # Qpr2 = sqrt(LCo * Bxr / (tau * Axr^2))  [sqrt(s^2*s^3/(s*s^4)) = sqrt(1) = 无量纲 ✓]
    Q_pr2 = np.sqrt(lc * Bx_r / (tau * Ax_r ** 2))
    Q_pr2 = float(max(Q_pr2, 1e-8))

    G_z = 1.0 + s / omega_zr
    G_p1 = 1.0 / (1.0 + s / (Q_pr1 * omega_pr1) + (s / omega_pr1) ** 2)
    G_p2 = 1.0 / (1.0 + s / (Q_pr2 * omega_pr2) + (s / omega_pr2) ** 2)
    Floop = A0_r * G_z * G_p1 * G_p2

    # 论文式(40) 近似：ωc ≈ ωzr * sqrt(A0r) / ωpr1 ... 
    # 实际由 Bode 数值确定；这里仅给出参考近似
    omega_c_approx = omega_zr * np.sqrt(1.0 / A0_r) * omega_pr1  # 占位，由UI从Bode图取实测值
    # PM 近似（论文式41）：PM ≈ arctan(Axr*(tau/lc)^2)
    PM_approx_rad = np.arctan(Ax_r * (tau / lc) ** 2)

    # 论文式(44) RxCx 设计范围
    denom_shape = np.sqrt((Ts / np.pi) ** 2 + (Ton / np.pi) ** 2)
    rxcx_max_paper = lc / max(denom_shape, 1e-30)
    rxcx_min = np.pi * lc / (5.0 * Ts)

    meta = dict(
        Ax_r=Ax_r, Bx_r=Bx_r, A0_r=A0_r, omega_zr=omega_zr,
        omega_pr1=omega_pr1, Q_pr1=Q_pr1, omega_pr2=omega_pr2, Q_pr2=Q_pr2,
        omega_c_approx=omega_c_approx,
        PM_approx_rad=PM_approx_rad,
        rxcx_min=rxcx_min,
        rxcx_max=rxcx_max_paper,
    )
    return Floop, meta


def interp_zero_db_crossing(freq, mag_db):
    for i in range(len(mag_db) - 1):
        if mag_db[i] > 0.0 and mag_db[i + 1] <= 0.0:
            denom = mag_db[i] - mag_db[i + 1]
            if abs(denom) < 1e-15:
                return freq[i], i
            t = mag_db[i] / denom
            fc = freq[i] + t * (freq[i + 1] - freq[i])
            return fc, i
    return None, None


def bode_unwrap_deg(Floop, w):
    mag, phase, omega = ct.bode(Floop, w, plot=False)
    freq = omega / (2.0 * np.pi)
    phase_deg = phase * 180.0 / np.pi
    phase_deg = np.unwrap(phase_deg * np.pi / 180.0) * 180.0 / np.pi
    return freq, mag, phase_deg
