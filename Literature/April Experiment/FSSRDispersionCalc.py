import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sc
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# known quantities
a0 = 750
Ecentral = 1600
R = 155.04
n = 2 # second order
twod = 19.84/n # lattice spacing of mica in angstrom
theta0 = np.arcsin(12398/(Ecentral*twod)) # calculate central angle from bragg equation
phi0 = np.pi/2-theta0
b0 = ((2*np.cos(phi0))/R - 1/a0)**(-1) # follows from focusing condition on saggital plane
ld = 27.6 # detector chip length
l = 50 # length of crystal in mm
nRays = 4 # number of rays to map
diff = 2 # how close the rays should be when determining central ray for detector shifting calc. [mm]
# for FSSR 1D
a0 = R*np.cos(phi0)/np.cos(2*phi0)
b0 = R*np.cos(phi0)
print(a0, b0, theta0*180/np.pi)
fig = plt.figure()
plt.axis('off')
# ax = plt.axes()
# ax.set_xlabel("x [mm]")
# ax.set_ylabel("y [mm]")
# line calculations
def ycalc (x, m, b):
    y = m*x + b
    return y

# for circles

t = np.linspace(0, 2*np.pi, 200)
plt.plot(R*np.cos(t), R*np.sin(t), linestyle = "dashed", color = "black", alpha = 0.5)
plt.plot(R/2*np.cos(t), -R/2 + R/2*np.sin(t), linestyle = "dashed", color = "black", alpha = 0.5)

# and plotting crystal
xcrystal = np.linspace(-l/2,l/2,500)
plt.plot(xcrystal, -np.sqrt(R**2-xcrystal**2), color = "blue", linewidth = 2)


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
plt.scatter(x0, y0, color = "red", alpha = 1, zorder = 100)

def fullcalc(E, colorplt):
    # make sure rays around central energy are red
    if (np.abs(E - Ecentral) <= 0.3):
        colorplt = "red"
    # next do calculations for intial ray and reflected ray
    m, mr = mcalc(E)
    ys = ycalc(xscalc(m,bcalc(m)), m, bcalc(m))
    xs = xscalc(m, bcalc(m))
    # break calculation if ray doesn't fall on crystal anymore
    if (np.abs(xs - xs0) > l/2 or len(EonCrystal) > 0):
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
    xr = np.linspace(xs, xf, 200)
    plt.plot(xr, ycalc(xr,mr,br), color = colorplt, linewidth = 0.6)

    x = np.linspace(x0, xs, 200)
    plt.plot(x, ycalc(x, m, bcalc(m)), color = colorplt, linewidth = 0.6)

    # # zoomed part
    # axins = zoomed_inset_axes(ax, 2, loc=1) # zoom = 2
    # axins.plot(xr, ycalc(xr,mr,br), color = colorplt, linewidth = 0.6)
    # axins.plot(x, ycalc(x, m, bcalc(m)), color = colorplt, linewidth = 0.6)
    # axins.set_xlim(-50, 50)
    # axins.set_ylim(-150, -50)
    # plt.xticks(visible=False)
    # plt.yticks(visible=False)
    # mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    # return xf0 and yf0 to use in detector calculation
    if (E == Ecentral):
        return xf, yf

# now run through all the energies
Erange = np.linspace(Ecentral-250, Ecentral+250, nRays)
xf0, yf0 = fullcalc(Ecentral, "red")
for E in Erange:
    fullcalc(E, 'purple')


# graph the dectector line as well
# md, bd = np.polyfit(dectX, dectY, 1)
# xgraphdect = np.linspace(x0, dectX[-1], 100)
# plt.plot(xgraphdect, md*xgraphdect + bd, linestyle = "dashed", color = "gray", alpha = 0.5)

# now show how far detector goes, xd is location of bottom end of detector chip. xt is top of detector chip.
# func = lambda xd : (xf0 - xd)**2 + (yf0 - md*xd - bd)**2 - (ld/2)**2
# xd = sc.fsolve(func, 80)
# xt = xd - 2*(xd - xf0)
# # want to offset chip so that equal numbers of enegies hit on both sides of central energy
# numAbove = 0
# numBelow = 0
# for i in dectX:
#     if (i > xf0):
#         numBelow += 1
#     if (i < xf0):
#         numAbove += 1

# topIndex = 0
# botIndex = 0
# for j in dectX: # effectively scans from top to bottom of detector, since Emin on top
#     if (j >= xt and j != xf0):
#         topIndex = dectX.index(j)
#         # print(topIndex)
#         break
# botIndex doesn't matter anymore because the number of rays above and below central energy will never be equal.
# Just set the shift in detector so that as many higher E's are recorded as possible.
# for j in dectX:
#     if (j > xd and j != xf0):
#         botIndex = dectX.index(j) - 1
#         print(botIndex)
#         break
# topIndex += round(np.abs((numAbove - numBelow)/2))
# botIndex += round(np.abs((numAbove - numBelow)/2))
# print (topIndex, botIndex)

# print("Number of rays above E central: ", numAbove,"; Number of rays below E central: ", numBelow, "; Number of rays shown in graph: ", len(dectX))
# # print(xf0, dectX[numAbove])
# offset = dectX[topIndex] - xt
# xdectline = np.linspace(xt + offset, xd + offset, 100)
# # plots detector
# plt.plot(xdectline, md*xdectline + bd, color = 'black', linewidth = 2.5)
# check if detector length is correct
# print(np.sqrt((xd - xt)**2 + (ycalc(xd, md, bd) - ycalc(xt, md, bd))**2))
# print("X offset of detector: ", offset)





# print(EonCrystal[0], EonCrystal[-1])
# # this is plot for diff of focus point to detector line
# diff = []
# for i in range (len(EonCrystal)):
#     xb = (bd-brs[i])/(mrs[i]-md)
#     thisdiff = np.sqrt((xb-dectX[i])**2+(ycalc(xb, md, bd)-dectY[i])**2)
#     diff.append(thisdiff)
# # plt.plot(EonCrystal, diff)

# # calculate dispersion on the detector
# EonCrystal.append(Ecentral)
# dectX.append(xf0[0])
# dectY.append(yf0[0])
# # d as norm of diff between central E vector and other E vector. E lower then Ecentral makes d neg. This is 
# # to bring it in line with simulation convention
# d = []
# for i in range(len(dectX)):
#     xdi = dectX[i]
#     ydi = dectY[i]
#     dvalue = np.sqrt((xdi-xf0)**2 + (ydi-yf0)**2)
#     if (xdi < xf0):
#         dvalue *= -1
#     d.append(dvalue[0])
# fit quadratically
# fit = np.poly1d(np.polyfit(d, EonCrystal, 2))
# plt.figure()
# plt.scatter(d, EonCrystal, s = 1/2)
# plt.plot(d, fit(d), color = "r")
# print(fit)

# plt.axis([-R/2-10, R/2+10, -R-10, 10])
plt.show()