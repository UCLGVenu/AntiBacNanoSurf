import numpy as np
import matplotlib.pyplot as plt
from Bending_Limits2 import func1
import Light_Absorptions_Test2 as LAT
import function_test2 as FT

# Func 1 return Radius, Max Heights, Energies
# Accepts no values (at the moment) - edit Bending_Limits.py to change variables

# Func 2 returns Dictionary - Key = Radius, Value = Arr[Heights, Absorption]
# Accepts arr(Radius), arr(Max Heights)

Rs, Hs, Es = func1()
results = LAT.func2(Rs, Hs)

legend = []
for k, v in results.items():
    #plt.plot(v[0], v[1])
    legend.append(str(k))

#plt.legend(legend)
#plt.show()

for k, v in results.items():
    heights = v[0]
    absorbs = v[1]
    rad = float(k)
    E = [FT.calculate_e(rad, p) for p in heights]
    normalised = [a * b for a,b in zip(absorbs, E)]
    plt.plot(heights, normalised)

plt.legend(legend)
plt.xlabel('Height / nm')
plt.ylabel('Efficiency Figure of Merit / AU')
plt.title('Bactericidal Figure of Merit for Radii/Heights')
plt.show()
