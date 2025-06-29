import numpy as np
import matplotlib.pyplot as plt

# a0 = 549.706 # crystal to source
# l = 50 # length of crystal
# R = 155.04 # radius of curvature
# n = 2 # second order
# # twod = 2.024/n
# twod = 19.84 # lattice spacing of mica in angstrom
# Ecentral = 1600 # in eV
# nRays = 200 # number of rays to map
# C = 12398
# theta0 = np.arcsin(12398/(Ecentral*twod/n))

Ecentral = 1600
R = 150
n = 2 # second order
twod = 19.84/n # lattice spacing of mica in angstrom
theta0 = np.arcsin(12398/(Ecentral*twod)) # calculate central angle from bragg equation
phi0 = np.pi/2-theta0
a0 = R*np.cos(phi0)/np.cos(2*phi0)
b0 = ((2*np.cos(phi0))/R - 1/a0)**(-1) # follows from focusing condition on saggital plane
print(a0, b0)
ld = 27.6 # detector chip length
l = 50 # length of crystal in mm
C = 12398
nRays = 1000


def r(E):
    s = np.sqrt(a0**2+R**2-2*a0*R*n*C/twod*1/Ecentral)
    r = R*n*C/twod*1/E + np.sqrt((R*n*C/twod)**2/E**2 + s**2 - R**2)
    return r
x0 = np.sqrt(r(Ecentral)**2 + R**2 - 2*R*n*C/twod*r(Ecentral)/Ecentral)/(1-R*twod/(2*n*C)*Ecentral/r(Ecentral))
print(x0)
print(a0*np.sin(2*theta0))
def func(E):
    d = np.sqrt(r(E)**2 + R**2 - 2*R*n*C/twod*r(E)/E)/(1-R*twod/(2*n*C)*E/r(E)) - x0
    return d


t = np.linspace(1450, 1805, nRays)
fit = np.poly1d(np.polyfit(t, func(t), 2))
fit = np.poly1d(np.polyfit(func(t), t, 2))
print(fit)
plt.plot(t,func(t), color = 'red', label = 'Analytical Calculation', zorder = 10)
plt.show()
# from simulation
def d(E):
    a = -0.02194
    b = 10.57
    c = 1600
    d = -b/(2*a) + np.sqrt(b**2-4*a*(c-E))/(2*a)
    return d

t = np.linspace(1430, 1820, nRays)
plt.plot(t, d(t), color = 'blue', label = 'mmpxrt Simulation')
plt.ylabel("d [mm]")
plt.xlabel("E [eV]")
#plt.show()

#----------------------------------------------------------------------------------
#From ray tracing
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sc

# known quantities
Ecentral = 1600
R = 150
n = 2 # second order
twod = 19.84/n # lattice spacing of mica in angstrom
theta0 = np.arcsin(12398/(Ecentral*twod)) # calculate central angle from bragg equation
phi0 = np.pi/2-theta0
b0 = ((2*np.cos(phi0))/R - 1/a0)**(-1) # follows from focusing condition on saggital plane
ld = 27.6 # detector chip length
l = 50 # length of crystal in mm
diff = 2 # how close the rays should be when determining central ray for detector shifting calc. [mm]
# for FSSR 1D
a0 = R*np.cos(phi0)/np.cos(2*phi0)
b0 = R*np.cos(phi0)
# ax = plt.axes()
# ax.set_xlabel("x [mm]")
# ax.set_ylabel("y [mm]")
# line calculations
def ycalc (x, m, b):
    y = m*x + b
    return y



# for initial path from source to crystal and for slope after reflection
def bcalc (m):
    b = a0*(np.sin(theta0)+m*np.cos(theta0))-R
    return b
def xscalc(m,b):
    xs = -m*b/(m**2+1) + np.sqrt((m*b/(m**2+1))**2 + (R**2-b**2)/(m**2+1))
    return xs
def mcalc (E):
    theta = np.arcsin(12398/(E*twod))
    func = lambda m : np.tan(theta) - np.abs((m-xscalc(m,bcalc(m))/np.sqrt(R**2-xscalc(m,bcalc(m))**2))/(1+m*xscalc(m,bcalc(m))/np.sqrt(R**2-xscalc(m,bcalc(m))**2)))
    m = sc.fsolve(func, -1)
    mr = np.tan(np.arctan(m)+2*theta)
    # x = np.linspace(-10, 10, 201)
    # plt.plot(x, func(x))
    return m,mr

# now do the full calculation and mapping. First initialize detector line info for graph. Also define Schnittpunkt for central wavelength. Comes from choice of coordinate system
dectX = []
dectY = []
EonCrystal = []
mrs = []
brs = []
xs0 = 0
ys0 = -R
# source position first, and plot it after everything else so that it covers
x0 = -a0*np.cos(theta0); y0 = a0*np.sin(theta0)-R

def fullcalc(E, colorplt):
    # make sure rays around central energy are red
    if (np.abs(E - Ecentral) <= 0.3):
        colorplt = "red"
    # next do calculations for intial ray and reflected ray
    m, mr = mcalc(E)
    ys = ycalc(xscalc(m,bcalc(m)), m, bcalc(m))
    xs = xscalc(m, bcalc(m))
    # break calculation if ray doesn't fall on crystal anymore
    if (np.abs(xs - xs0) > l/2):
        return
    br = ys - mr*xs
    # now look at focusing point
    theta = np.arcsin(12398/(E*twod))
    ai = np.sqrt((xs-x0)**2+(ys-y0)**2)
    bs = (2*np.cos(np.pi/2-theta)/R - 1/ai)**(-1)
    xf = bs/np.sqrt(1+mr**2) + xs
    yf = mr*xf + br
    # store focusing points for later use. Condition of if statement for later calculation need mr and br for detector deviation calculation.
    if (E != Ecentral):
        dectX.append(xf[0])
        dectY.append(yf[0])
        mrs.append(mr[0])
        brs.append(br[0])
        EonCrystal.append(E)
    # plot the lines, where reflected line only runs until the saggital focusing condition
    if (E == Ecentral):
        return xf, yf

# now run through all the energies
Erange = t
xf0, yf0 = fullcalc(Ecentral, "red")
for E in Erange:
    fullcalc(E, 'purple')


# calculate dispersion on the detector
EonCrystal.append(Ecentral)
dectX.append(xf0[0])
dectY.append(yf0[0])
# d as norm of diff between central E vector and other E vector. E lower then Ecentral makes d neg. This is 
# to bring it in line with simulation convention
d = []
for i in range(len(dectX)):
    xdi = dectX[i]
    ydi = dectY[i]
    dvalue = np.sqrt((xdi-xf0)**2 + (ydi-yf0)**2)
    if (xdi < xf0):
        dvalue *= -1
    d.append(dvalue[0])
# fit quadratically
fit = np.poly1d(np.polyfit(d, EonCrystal, 2))
plt.plot(fit(d), d, color = "green", label = 'Simple Ray Tracing', zorder = 20)
print(EonCrystal[0], EonCrystal[-2])

plt.legend()
# plt.show()