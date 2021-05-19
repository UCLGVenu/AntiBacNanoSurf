import numpy as np
import matplotlib.pyplot as plt

# Fixed Constants

lamda = 365e-9
k0 = (2 * np.pi) / lamda
w0 = 3e8 * k0


def get_kx_layer(n_layer):
    kx_layer = k0 * np.sqrt((n_layer**2) - (np.sin(angle) ** 2))
    return kx_layer


def get_fs(p, r, coating):
    a_in = np.pi * (r ** 2)
    a_out = np.pi * ((r + coating) ** 2)
    box = p ** 2
    f_in = a_in / box
    f_coat = (a_out - a_in) / box
    return f_in, f_coat


def get_effective_n(n_host, n1, n2, f1, f2):
    e_host, e1, e2 = n_host ** 2, n1 ** 2, n2 ** 2
    q1 = f1 * (e1 - e_host) / (e1 + (2 * e_host))
    q2 = f2 * (e2 - e_host) / (e2 + (2 * e_host))
    q = q1 + q2
    e_eff = e_host * ((1 + (2 * q)) / (1 - q))
    return np.sqrt(e_eff)


def get_rt(k1, k2):
    r = (k1 - k2) / (k1 + k2)
    t = (2 * k1) / (k1 + k2)
    return r, t


def get_rt_pn(k1, k2, n1, n2):
    r = ((n2 ** 2 * k1) - (n1 ** 2 * k2)) / ((n1 ** 2 * k2) + (n2 ** 2 * k1)) * -1
    t = (2 * (n1 * n2) * k1) / ((n1 ** 2 * k2) + (n2 ** 2 * k1))
    #t = (n1 / n2) * (r + 1)
    return r, t


# Set variable here. Requires changing of 'var' and variable in loop.
# BE CAREFUL for power of 10 factors. Will result in weird looking graphs otherwise.

angles = np.arange(0, 90, 1)
coatings = np.arange(0, 15.5, 0.5) * 1e-9
heights = np.arange(300, 1020, 1) * 1e-9
pitches = np.arange(80, 505, 5) * 1e-9
top_fraction = np.arange(0.02, 1.02, 0.02)
bot_rs = np.arange(20, 92, 2) * 1e-9
var = bot_rs

Rs = []
Ts = []
As = []

for q in range(len(var)):

    # Initialise non-fixed variables.

    coat = 10e-9 # Thickness of APPLIED COATING
    angle = (np.pi / 180) * 0

    pitch = 200e-9
    bot_r = bot_rs[q]  # Radius of INNER PILLAR (Substrate only)
    top_r = bot_r
    height = 500e-9
    kz = k0 * np.sin(angle)

    n_air = 1
    n_si = np.complex(6.52, -2.67)
    #n_si = np.complex(1.49, 0)
    n_ti = np.complex(2.88, -0.017)
    kx_air = k0 * np.sqrt((1 ** 2) - (np.sin(angle) ** 2))
    kx_si = k0 * np.sqrt((n_si**2) - (np.sin(angle) ** 2))
    kx_si2 = np.sqrt(((n_si**2 - 1)*(k0**2)) + kz**2)

    sections = 100  # Set number of sections to iterate over - note #iterations is multiplied by len(var). Can get slow.
    depth = height / sections
    step = (bot_r - top_r) / sections
    r_val = bot_r

    f1, f2 = get_fs(pitch, r_val, coat)
    n_bot = get_effective_n(n_air, n_si, n_ti, f1, f2)
    kx_layer1 = get_kx_layer(n_bot)
    kx_layer2 = get_kx_layer(n_si)
    r1, t1 = get_rt(kx_layer1, kx_layer2)
    rp1, tp1 = get_rt_pn(kx_layer1, kx_layer2, n_bot, n_si)

    matrix = np.array([[1, r1], [r1, 1]]) * (1 / t1)
    matrix2 = np.array([[1, rp1], [rp1, 1]]) * (1 / tp1)

    kx_layer2 = kx_layer1
    r_val -= step

    for i in range(sections - 1):
        exponent = kx_layer2 * depth * 1j
        p_matrix = np.array([[np.exp(exponent), 0], [0, np.exp(-1 * exponent)]])

        matrix = np.matmul(p_matrix, matrix)
        matrix2 = np.matmul(p_matrix, matrix2)

        n_previous = get_effective_n(n_air, n_si, n_ti, f1, f2)
        f1, f2 = get_fs(pitch, r_val, coat)
        n_layer = get_effective_n(n_air, n_si, n_ti, f1, f2)
        kx_layer1 = get_kx_layer(n_layer)

        r2, t2 = get_rt(kx_layer1, kx_layer2)
        rp2, tp2 = get_rt_pn(kx_layer1, kx_layer2, n_layer, n_previous)

        t_matrix = np.array([[1, r2], [r2, 1]]) * (1 / t2)
        t_matrix2 = np.array([[1, rp2], [rp2, 1]]) * (1 / tp2)

        matrix = np.matmul(t_matrix, matrix)
        matrix2 = np.matmul(t_matrix2, matrix2)

        kx_layer2 = kx_layer1
        r_val -= step

    exponent = kx_layer2 * depth * 1j
    p_matrix = np.array([[np.exp(exponent), 0], [0, np.exp(-1 * exponent)]])
    matrix = np.matmul(p_matrix, matrix)
    matrix2 = np.matmul(p_matrix, matrix2)

    f1, f2 = get_fs(pitch, r_val, coat)
    n_top = get_effective_n(n_air, n_si, n_ti, f1, f2)
    kx_layer = get_kx_layer(n_top)
    r3, t3 = get_rt(kx_air, kx_layer)
    rp3, tp3 = get_rt_pn(kx_air, kx_layer, 1, n_top)

    t_matrix = np.array([[1, r3], [r3, 1]]) * (1 / t3)
    t_matrix2 = np.array([[1, rp3], [rp3, 1]]) * (1 / tp3)

    matrix = np.matmul(t_matrix, matrix)
    matrix2 = np.matmul(t_matrix2, matrix2)

    r = matrix[1][0] / matrix[0][0]
    t = 1 / matrix[0][0]

    rp = matrix2[1][0] / matrix2[0][0]
    tp = 1 / matrix2[0][0]

    s_in = 0.5 * (kx_air / w0)
    s_r = 0.5 * (kx_air / w0) * (r.real ** 2 + r.imag ** 2)
    #s_t = 0.5 * (kx_si / w0) * (t.real ** 2 + t.imag ** 2)
    s_t = 0.5 * ((kx_si / w0) * t * np.conj(t)).real

    # Poynting for p

    theta_out = np.arctan(kz/kx_si)

    #rfac = rp * np.conj(rp)
    #tfac = tp * np.conj(tp)

    '''
    p_in = 0.5 * (1/w0) * ((kz * np.cos(angle) * np.conj(np.sin(angle))) - (np.conj(kx_air) * np.cos(angle) * np.conj(np.cos(angle))))
    p_r = 0.5 * (rfac/w0) * ((kz * np.cos(angle) * np.conj(np.sin(angle))) - (kx_air * np.cos(angle) * np.conj(np.cos(angle))))
    p_t = 0.5 * (tfac/w0) * ((kz * np.cos(theta_out) * np.conj(np.sin(theta_out))) - (np.conj(kx_air) * np.cos(theta_out) * np.conj(np.cos(theta_out)))).real
    '''

    p_in2 = 0.5 * (kx_air)
    p_r2 = 0.5 * (kx_air) * (rp.real ** 2 + rp.imag ** 2)
    p_t2 = 0.5 * (kx_si * tp * np.conj(tp)).real
    #p_t2 = 0.5 * ((kx_air / w0) * (tp*np.conj(tp))).real

    R = s_r / s_in
    T = s_t / s_in

    #Rp = p_r / p_in
    #Tp = p_t / p_in

    Rp2 = p_r2 / p_in2
    Tp2 = p_t2 / p_in2

    #print(R, T, 1 - (R + T))
    #print(Rp2, Tp2, 1 - (Rp2 + Tp2))

    R_avg = (R + Rp2) / 2
    T_avg = (T + Tp2) / 2
    A_avg = ((1-(R+T)) + (1-(Rp2+Tp2))) / 2

    Rs.append(R_avg)
    Ts.append(T_avg)
    As.append(A_avg)

# Get value of maxima for oscillatory graphs

for i in range(1, len(Rs) - 1, 1):
    if Rs[i - 1] < Rs[i] and Rs[i + 1] < Rs[i]:
        print(var[i])

plt.plot(var, Rs)
plt.plot(var, Ts)
plt.plot(var, As)
plt.legend(['Reflection', 'Transmission', 'Absorption'])
plt.title('R,T,A varying Radius - Angle 0, Height 500, Coat 10, Pitch 200')
plt.ylabel('Intensity')
plt.xlabel('Radius / m')
plt.show()

