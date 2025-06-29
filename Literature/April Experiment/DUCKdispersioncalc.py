import numpy as np
import matplotlib.pyplot as plt

a0 = 549.706 # crystal to source
L = 117.584*2 # spectrometer length
n = 1 # first order
twod = 10.64/n # lattice spacing of ADP in angstrom
Ecentral = 1580 # in eV
C = 12398 # Umrechnungsfaktor
theta0 = np.arcsin(12398/(Ecentral*twod))

def func(E):
    d = L*(1/np.sqrt((twod/C)**2*E**2-1) - 1/np.sqrt((twod/C)**2*Ecentral**2-1))
    return d

t = np.linspace(1541, 1620, 300)
plt.plot(t,func(t), color = 'red', linewidth = 2.5, label = 'Analytical Calculation')

# from simulation
def d(E):
    a = 0.01402
    b = -2.79
    c = 1580
    d = -b/(2*a) - np.sqrt(b**2-4*a*(c-E))/(2*a)
    return d
#plt.plot(t, d(t), color = 'blue', linewidth = 2.5, label = 'mmpxrt Simulation')
plt.legend()
plt.ylabel("d [mm]")
plt.xlabel("E [eV]")

fit = np.poly1d(np.polyfit(t, func(t), 2))
print(fit)
plt.plot(t, fit(t))
# plt.show()

print(func(1598.4))