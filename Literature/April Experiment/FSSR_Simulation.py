from numpy import cos, sin, arcsin, arccos, pi, sqrt, arctan

runSim = True # decide whether or not to run simulation
# simulation parameters
SimName = "FSSR-1D Unfocused pure reflection"
SimLength = 5 # amount of rays, exponent of exponential
a = 549.706 # crystal to source
l = 50 # length of crystal
R = 155.04 # radius of curvature
n = 2 # second order
# twod = 2.024/n
twod = 19.84/n # lattice spacing of mica in angstrom
Ecentral = 1600 # in eV


# from here is calculation

thetaC = arcsin(12398/(Ecentral*twod)) # calculate central angle from bragg equation
phi = pi/2-thetaC
b = ((2*cos(phi))/R - 1/a)**(-1) # follows from focusing condition on saggital plane
x = sqrt((a)**2+(R)**2-2*a*R*cos(phi))
s1min = cos(thetaC)*a-l/2
s1max = cos(thetaC)*a+l/2
s2 = a*sin(thetaC)-(R-sqrt(R**2-(l/2)**2))
amin = sqrt(s1min**2+s2**2)
amax = sqrt(s1max**2+s2**2)

phimin = arccos(1/(2*amin*R)*(amin**2+R**2-x**2))
phimax = arccos(1/(2*amax*R)*(amax**2+R**2-x**2))
phimin = (90-45.4)*pi/180
phimax = (90-58.53)*pi/180
Emax = 12398/(sin(pi/2-phimin)*twod)
Emin = 12398/(sin(pi/2-phimax)*twod)

Bandwidth = Emax-Emin
Emath = Bandwidth/2+Emin

# calculate angle of detector, for which the detector is tangential to roland circle (hopefully gives optimal imaging and spectra resolution)
# beta = arctan(sin(phi)/((2*b)/R - cos(phi)))
# beta = beta*180/pi
# result: this doesnt seem to help much, as it leads to an even more intense curve in the dispersion direction
# next: try to align the meridional focussing conditions instead. Previous calculation with saggital focusing condition should be valid
# reasoning: basically got a point source, ie saggital focusing should not play a huge role, considering distance differences due to detector tilt.
# --> best spectral resolution when focusing along a line in dispersive plane.


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

# print(bs, bm, bsmin, bmmin, bsmax, bmmax)

# print(bs, bsmin, bsmax, bm, bmmin, bmmax)

# simulation

if (runSim):

    from mmpxrt import mmpxrt

    p=mmpxrt.init()
    p['simulation']['numraysE']=SimLength
    p['simulation']['name']=SimName
    p['simulation']['numerical_intersect_halving']=1 # prevents infinite loop due to stong crystal curve
    p['simulation']['num_processes']=1

    p['source']['EcentralRay']=Ecentral
    p['source']['EmaxBandwidth']= Bandwidth + 500 # make sure to get max bandwidth
    p['source']['size']= 100e-3 # set source size to cube with length of 100 micrometer

    p['crystal']['d2']=twod
    p['crystal']['width']=10
    p['crystal']['length']=l
    p['crystal']['maxThickness']=3.5
    p['crystal']['thickness']=0 # use exponential penetration
    p['crystal']['radius_w']=R
    p['crystal']['radius_l']=R
    p['crystal']['mosaicity']=0 # 0 bc monocrystal
    # p['crystal']['crystalliteRockingcurveWidth']=rctest # this in radians, values taken from paper # leave out bc monocrystal
    p['crystal']['rockingCurveFWHM']=0.002322 # include bc monocrystal
    p['crystal']['integrated_reflectivity']=0.0000289 # both taken from Hölzer, in FSSR data folder. Is an old paper, so value should be higher than
    # now, which means if resolution is acceptable in simulation, should be great in experiment. 

    p['geometry']['detRot']=0
    p['geometry']['CrystalSource']=a
    p['geometry']['CrystalDetector']=b+5 # move camera 5 mm away from crystal
    p['geometry']['detectorLength']=27.6
    p['geometry']['detectorWidth']=6.9
    #p['geometry']['detectorPxSize']=200
    # p['geometry']['ThBragg0']=thetaC # force central bragg angle, since doing from Emath otherwise
    rrrs = mmpxrt.spectrometer(p)
    mmpxrt.spectrometer_evaluate(p,rrrs)


print('Central angle is: ', thetaC*180/pi)
print('Mathematical central energy is ',Emath)
print('Bandwidth is: ', Bandwidth)
print('Min Angle is: ', (pi/2-phimin)*180/pi)
print('Max Angle is: ', (pi/2-phimax)*180/pi)
print('Min Energy is: ', Emin)
print('Max Energy is: ', Emax)
# print('Detector rotation angle in degrees is ', -beta)