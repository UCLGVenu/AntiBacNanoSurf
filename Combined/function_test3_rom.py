import numpy as np


def calculate_e(pitch_in, height):
    h = height
    force = 25e-9
    pitch = pitch_in
    coat = 10e-9
    E_ti = 157e9
    nu_ti = 0.35

    rad = 30e-9

    E_si_d = [12, 38.5, 132, 165]
    E_si_val = [53, 68, 132, 165]

    d = 2 * (rad * 1e9)

    v_si = (rad ** 2) / ((rad + coat) ** 2)
    v_ti = 1 - v_si

    E_si = -1
    nu_si = 0.27
    nu = (v_si * nu_si) + (v_ti * nu_ti)

    for i in range(3):
        if d >= 165:
            E_si = 165 * 1e9
            break
        if E_si_d[i] < d < E_si_d[i + 1]:
            dif = (d - E_si_d[i]) / (E_si_d[i + 1] - E_si_d[i])
            E_si = (E_si_val[i] + (dif * (E_si_val[i + 1] - E_si_val[i]))) * 1e9
            break

    E_si = 1 / ((v_si / E_si) + (v_ti / E_ti))

    T_tilt = (1.3 / (2 * np.pi)) * (1 + nu) * ((2 * (1 - nu)) + (1 - (1 / (4 * (1 - nu)))))

    ld = h / (2 * rad)

    k_base = ((16 / 3) * (ld ** 3)) + (((7 + 6 * nu) / 3) * ld) + (8 * T_tilt * ld ** 3)
    k_top = (np.pi / 4) * (2 * rad) * E_si
    k = k_top / k_base

    defl = abs(force / k)
    energy = 0.5 * k * (defl ** 2)
    avg_defl = 0.3752 * defl

    separation = pitch - ((2 * defl) + (2 * rad))  # Worst Case Scenario
    separation_flat = pitch - (defl + (2 * rad))  # Flat Case Scenario
    max_deflection = (pitch - (2 * rad)) / 2
    max_energy = 0.5 * k * (max_deflection ** 2)

    return energy
