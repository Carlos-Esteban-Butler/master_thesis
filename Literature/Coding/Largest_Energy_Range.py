from numpy import cos, sin, arcsin, arccos, pi, sqrt

a = 1000 # crystal to source
l = 50 # length of crystal
R = 150 # radius of curvature
n = 1
twod = 2.024/n
Ecentral = 9537
thetaC = arcsin(12398/(Ecentral*twod)) # calculate central angle from bragg equation
print(thetaC*180/pi)
phi = pi/2-thetaC
x = sqrt((a)**2+(R)**2-2*a*R*cos(phi))
s1min = cos(thetaC)*a-l/2
s1max = cos(thetaC)*a+l/2
s2 = a*sin(thetaC)-(R-sqrt(R**2-(l/2)**2))
amin = sqrt(s1min**2+s2**2)
amax = sqrt(s1max**2+s2**2)

phimin = arccos(1/(2*amin*R)*(amin**2+R**2-x**2))
phimax = arccos(1/(2*amax*R)*(amax**2+R**2-x**2))

Emax = 12398/(sin(pi/2-phimin)*twod)
Emin = 12398/(sin(pi/2-phimax)*twod)

Bandwidth = Emax-Emin
print('Mathematical central energy is ',Bandwidth/2+Emin)
print('Bandwidth is: ', Bandwidth)
print('Min Angle is: ', (pi/2-phimin)*180/pi)
print('Max Angle is: ', (pi/2-phimax)*180/pi)
print('Min Energy is: ', Emin)
print('Max Energy is: ', Emax)

# define function to calculate saggital and meridional focus point. Different because of astigmatism
def saggital(phi, a):
    return (2*cos(phi)/R - 1/a)**(-1)
def meridional (phi, a):
    return (2/(R*cos(phi)) - 1/a)**(-1)
# calculate sag and mer focus points for central and min angles
bs = saggital(phi, a)
bm = meridional(phi, a)
bsmin = saggital(phimin, amin)
bmmin = meridional(phimin, amin)
bsmax = saggital(phimax, amax)
bmmax = meridional(phimax, amax)

print(bs, bsmin, bsmax, bm, bmmin, bmmax)