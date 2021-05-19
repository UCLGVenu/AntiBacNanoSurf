import numpy as np
import matplotlib.pyplot as plt

# Define Constants

pitch = 200e-9
radius = 30e-9
#coating = 5e-9
force = 25e-9
height = 500e-9

# E_si
# @300nm = ~165
# @188nm = ~132
# @38.5nm = 68
# @12nm = 53

poisson_si = 0.27
poisson_ti = 0.35

E_si_d = [12, 38.5, 132, 165]
E_si_val = [53, 68, 132, 165]

d = 2 * (radius * 1e9)

for i in range(3):
    if E_si_d[i] <= d < E_si_d[i+1]:
        dif = (d - E_si_d[i])/(E_si_d[i+1] - E_si_d[i])
        E_si = (E_si_val[i] + (dif * (E_si_val[i+i] - E_si_val[i])))* 1e9

E_ti = 157e9

coatings = np.arange(0,20,1) * 1e-9
deflections = []
energies = []

for i in range(20):

    G_si = E_si / (2*(1+poisson_si))
    G_ti = E_ti / (2*(1+poisson_ti))

    v_si = (radius ** 2)/((radius+coatings[i])**2)
    v_ti = 1 - v_si

    G_comp = (G_si * G_ti) / ((v_si * G_ti) + (v_ti * G_si))
    E_comp = (E_si * v_si) + (E_ti * v_ti)
    poisson_comp = (E_comp/(2*G_comp)) - 1
    rom_poisson = (v_si*poisson_si) + (v_ti*poisson_ti)

    sma = np.pi * 0.25 * (radius+coatings[i])**4
    area = np.pi * (radius+coatings[i])**2

    '''# Method 1 - My Bad Method
    
    k_bend = (3 * E_comp * sma) / height**3
    K_timo = (6 + (6*poisson_comp)) / (7 + (6*poisson_comp))
    k_shear = (K_timo * G_comp * area) / height
    
    d_pillar = force * ((1/k_bend) + (1/k_shear))'''

    nu = poisson_comp

    T_tilt = 1.3 * ((1+nu)/(2*np.pi)) * ((2*(1-nu)) + (1 - (1/(4*(1-nu)))))
    d_tilt = 8 * T_tilt * ((height/(2*radius+coatings[i]))**2) * (4/np.pi) * (force/(E_comp * (2*radius+coatings[i])))

    #d_total = d_pillar + d_tilt


    # Version 2

    coeff_1 = 4 * force * 1/(2 * radius * np.pi)
    c1 = (16 / (3*E_comp)) * ((height/(radius*2))**3)
    c2 = ((7+(6*poisson_comp))/(3*E_comp)) * (height / (radius * 2))
    c3 = 8 * T_tilt * (1/E_comp) * ((height/(radius*2))**2)

    d_total2 = coeff_1 * (c1 + c2 + c3)

    k = force / d_total2

    U = 0.5 * k * (d_total2**2)

    deflections.append(d_total2)
    energies.append(U)

print(deflections)
plt.plot(deflections)
plt.ylabel('Maximum Deflection / m')
plt.xlabel('Coating Thickness / nm')
plt.title('Deflection of Si Pillar H500, R30, P 200, F 25nN with TiO2 Coating')
plt.show()

plt.plot(energies)
plt.ylabel('Energy stored / J')
plt.xlabel('Coating Thickness, nm')
plt.title('Energy of Si Pillar H500, R30, P 200, F 25nN with TiO2 Coating')
plt.show()