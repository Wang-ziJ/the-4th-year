# Interactive peak static overpressure plot using glasstone and matplotlib


import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Cursor
# Overpressure functions of U.S. and Soviet origin.

import numpy as np
from scipy.integrate import quad

'''the following code was copied from glasstone.utilities.py'''
# functions reused throughout glasstone module
class ValueOutsideGraphError(Exception):
    """This exception indicates that the requested input value falls outside the
graphs found in the original source."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# Unfortunately, historical NWE models emply a bewildering array of non-Si units.
# Even worse, many of the models utilize arbitrary combinations of units. This unit
# conversion function is provided to provide standard inputs in SI units for all
# parts of glasstone.

class UnknownUnitError(Exception):
    """This exception indicates that the requested unit to convert is unknown."""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
def convert_units(v, unitsfrom, unitsto):
    if unitsfrom == unitsto:
        return v
    # yield
    elif unitsfrom == 'kT' and unitsto== 'MT':
        return v / 1000.0
    elif unitsfrom == 'MT' and unitsto== 'kT':
        return v * 1000.0
    # distance
    elif unitsfrom == 'm' and unitsto== 'kilofeet':
        return v / 304.8
    elif unitsfrom == 'm' and unitsto== 'km':
        return v / 1000.0
    elif unitsfrom == 'km' and unitsto== 'm':
        return v * 1000.0
    elif unitsfrom == 'kilofeet' and unitsto== 'm':
        return 304.8 * v
    elif unitsfrom == 'yards' and unitsto== 'm':
        return v / 1.09361
    elif unitsfrom == 'm' and unitsto== 'yards':
        return v * 1.09361
    elif unitsfrom == 'ft' and unitsto== 'm':
        return v * 0.3048
    elif unitsfrom == 'm' and unitsto== 'ft':
        return v / 0.3048
    elif unitsfrom == 'kilofeet' and unitsto== 'km':
        return convert_units(v, 'kilofeet', 'm') / 1000.0
    elif unitsfrom == 'kilofeet' and unitsto== 'mi':
        return v / 5.28
    elif unitsfrom== 'mi' and unitsto== 'km':
        return v * 1.60934
    elif unitsfrom== 'km' and unitsto== 'mi':
        return v / 1.60934
    elif unitsfrom== 'km' and unitsto== 'kilofeet':
        return v / 0.3048
    elif unitsfrom== 'yards' and unitsto== 'meters':
        return v * 0.9144
    elif unitsfrom== 'yards' and unitsto== 'km':
        return v * 0.0009144
    elif unitsfrom== 'meters' and unitsto== 'yards':
        return v / 0.9144
    elif unitsfrom== 'km' and unitsto== 'yards':
        return v / 0.0009144
    #pressure
    elif unitsfrom == 'psi' and unitsto== 'kg/cm^2':
        return v * 0.070307
    elif unitsfrom == 'kg/cm^2' and unitsto== 'psi':
        return v / 0.070307
    elif unitsfrom == 'MPa' and unitsto== 'psi':
        return v * 145.037738
    elif unitsfrom == 'psi' and unitsto== 'MPa':
        return v / 145.037738
    elif unitsfrom == 'kg/cm^2' and unitsto== 'MPa':
        return convert_units(convert_units(v, 'kg/cm^2', 'psi'), 'psi', 'MPa')
    elif unitsfrom == 'MPa' and unitsto== 'kg/cm^2':
        return convert_units(convert_units(v, 'psi', 'kg/cm^2'), 'MPa', 'psi')
    elif unitsfrom =='Pa':
        return convert_units(v, 'MPa', unitsto) / 1e6
    elif unitsto == 'Pa':
        return convert_units(v, unitsfrom, 'MPa') * 1e6
    # speed
    elif unitsfrom == 'm/s' and unitsto== 'mph':
        return v * 2.23694
    elif unitsfrom == 'mph' and unitsto== 'm/s':
        return v / 2.23694
    elif unitsfrom == 'm/s' and unitsto== 'km/h':
        return v * 3.6
    elif unitsfrom == 'km/h' and unitsto== 'm/s':
        return v / 3.6
    elif unitsfrom == 'mph' and unitsto== 'km/h':
        return v * 1.60934
    elif unitsfrom == 'km/h' and unitsto== 'mph':
        return v / 1.60934
    # wind shear
    elif unitsfrom == 'm/s-km' and unitsto == 'mph/kilofoot':
        return v * 0.13625756613945836
    # dose
    # under normal circumstances this isn't quite right, as Roentgens were
    # usually employed as a unit of exposure rather than dose. However, WSEG-10
    # used an unusual unit, the Equivilent Residual Dose, which does convert
    # directly into Sv:
    elif unitsfrom == 'Roentgen' and unitsto == 'Sv':
        return v / 100.0
    else:
        raise UnknownUnitError((unitsfrom, unitsto))

def dict_reverse(d):
    new_dict = {}
    for k in d:
        new_dict[k] = d[k][::-1]
    return new_dict







'''the following part was modified from glassrone.overpressure.py'''
# First, some utility functions; these should go in a dedicated file eventually
def scale_range(bomb_yield, ground_range):
    return ground_range / (bomb_yield**(1.0 / 3))

def scale_height(bomb_yield, burst_height):
    return burst_height / (bomb_yield**(1.0 / 3))

# Overpressure function from pp. 60-71 of H.L. Brode, _Airblast From Nuclear Bursts-
# Analytic Approximations_ (Pacific-Sierra Research Corporation, 1986).
# Appears to have been made using Fourier transforms on data from an model of a blast
# wave over an ideal surface.  The many local functions it contains are of no general
# physical significance.

def _brode(z, r, y):
    """Brode equation for approximating peak static overpressure for a 1 kT burst.
    Units are kT (yield) and kilofeet (ground range, burst height).
    Caveats: This is accurate to ~10%, and presumes sea level ambient pressure.
    """
    def a(z):
        return 1.22 - ((3.908 * z**2) / (1 + 810.2 * z**5))
    def b(z):
        return 2.321 + ((6.195 * z**18) / (1 + 1.113 * z**18)) - ((0.03831 * z**17) / (1 + 0.02415 * z**17)) + (0.6692 / (1 + 4164 * z**8))
    def c(z):
        return 4.153 - ((1.149 *  z**18) / (1 + 1.641 * z **18)) - (1.1 / (1 + 2.771 * z**2.5))
    def d(z):
        return -4.166 + ((25.76 * z**1.75) / (1 + 1.382 * z**18)) + ((8.257 * z) / (1 + 3.219 * z))
    def e(z):
        return 1 - ((0.004642 * z**18) / (1 + 0.003886 * z**18))
    def f(z):
        return 0.6096 + ((2.879 * z**9.25) / (1 + 2.359 * z**14.5)) - ((17.5 * z**2) / (1 + 71.66 * z**3))
    def g(z):
        return 1.83 + ((5.361 * z**2) / (1 + 0.3139 * z**6))
    def h(z, r, y):
        return ((8.808 * z**1.5) / (1 + 154.5 * z**3.5)) - ((0.2905 + 64.67 * z**5) / (1 + 441.5 * z**5)) - ((1.389 * z) / (1 + 49.03 * z**5)) + ((1.094 * r**2) / ((781.2 - (123.4 * r) + (37.98 * r**1.5) + r**2) * (1 + (2 * y))))
    def j(y):
        return ((0.000629 * y**4) / (3.493e-9 + y**4)) - ((2.67 * y**2) / (1 + (1e7 * y**4.3)))
    def k(y):
        return 5.18 + ((0.2803 * y**3.5) / (3.788e-6 + y**4))
    return (10.47 / r**a(z)) + (b(z) / r**c(z)) + ((d(z) * e(z)) / (1 + (f(z) * r**g(z)))) + h(z, r, y) + (j(y) / r**k(y))

def _brodeop(bomb_yield, ground_range, burst_height):
    """Calculate overpressure for arbitrary-height airbursts using Brode 
    equation.

    Units: kT, kilofeet
    
    Warning: ground_range = 0 results in divide-by-zero error.
    """
    z = (burst_height / ground_range)
    y = scale_height(bomb_yield, burst_height)
    x = scale_range(bomb_yield, ground_range)
    r = (x**2 + y**2)**0.5
    return _brode(z, r, y)

def brode_overpressure(y, r, h, yunits='kT', dunits='m', opunits='kg/cm^2'):
    """Estimate peak static overpressure at radius r from the epicenter of a burst
with yield y and a burst height of h using the Brode equation."""
    yld = convert_units(y, yunits, 'kT')
    ground_range = convert_units(r, dunits, 'kilofeet')
    height = convert_units(h, dunits, 'kilofeet')
    op = _brodeop(yld, ground_range, height)
    return convert_units(op, 'psi', opunits)


# Functions from NRDC, _The U.S. Nuclear War Plan: A Time for Change_
# These were, in turn, sourced from the help files of the Reagan-era
# Defense Nuclear Agency DOS programs BLAST and WE


# Altitude scaling factors (SP, SD, and ST)

def _altitude_t(alt):
    if 0 <= alt < 11000:
        return 1 - (2 * 10**9)**-0.5 * alt
    if 11000 <= alt < 20000:
        return 0.7535 * (1 + (2.09 * 10**-7) * alt)
    if alt >= 20000:
        return 0.684 *  (1 + (5.16 * 10**-6) * alt)

def _altitude_p(alt):
    if 0 <= alt < 11000:
        return _altitude_t(alt)**5.3
    if 11000 <= alt < 20000:
        return 1.6**0.5 * (1 + (2.09 * 10**-7) * alt)**-754 
    if alt >= 20000:
        return 1.4762 *  (1 + (5.16 * 10**-6) * alt)**-33.6

def _altitude_sp(alt):
    return _altitude_p(alt)

def _altitude_sd(alt):
    return _altitude_sp(alt)**(-1.0/3)

def _altitude_st(alt):
    return _altitude_sd(alt) * _altitude_t(alt)**-0.5

# 'Altitude-dependent speed of sound'
# 'rule of thumb: C increases 1.8% for each 10 deg. C rise above 15 deg. C'

def _altitude_speed_of_sound(alt):
    return (340.5 * _altitude_sd(alt)) / _altitude_st(alt)

# The Defense Nuclear Agency 1kT Free-Air Overpressure Standard
def _DNA1kTfreeairop(r):
    return (3.04 * 10**11)/r**3 + (1.13 * 10**9)/r**2 + (7.9 * 10**6)/(r * (np.log(r / 445.42 + 3 * np.exp(np.sqrt(r / 445.42) / -3.0)))**0.5)

# PFREE
def _DNAfreeairpeakop(r, y, alt):
    r1 = r / (_altitude_sd(alt) * y**(1.0/3))
    return _DNA1kTfreeairop(r1) * _altitide_sp(alt)

# r in these functions is not scaled for yield to facilitate their re-use for the
# airburst functions. 
def _shock_strength(op):
   return op / 101325 + 1

def _shock_gamma(op):
  xi = _shock_strength(op)
#  xi[xi <= 0.0] = 0.01

  t = 10**-12 * xi**6
  z = np.log(xi) - (0.47 * t) / (100 + t)
  return 1.402 - (3.4 * 10**-4 * z**4) / (1+ 2.22 * 10**-5 * z**6)

def _shock_mu(g):
    return (g + 1) / (g - 1)

def _mass_density_ratio(op):
    xi = _shock_strength(op)
    mu = _shock_mu(_shock_gamma(op))
    return (1 + mu * xi) / (5.975 + xi)

def _DNA1kTfreeairdyn(r):
    op = _DNA1kTfreeairop(r)
    return 0.5 * op * (_mass_density_ratio(op) - 1)

# QFREE
def _DNAfreeairpeakdyn(r, y, alt):
    r1 = r / (_altitude_sd(alt) * y**(1.0/3))
    return _DNA1kTfreeairdyn(r1) * _altitide_sp(alt)

def _scaledfreeairblastwavetoa(r):
    r2 = r * r
    return (r2 * (6.7 + r)) / (7.12e6 + 7.32e4 * r + 340.5 * r2)

def _freeairblastwavetoa(r, y, alt):
    return _scaledfreeairblastwavetoa(r) * _altitude_st(alt) * y**(1.0/3)

# Rankine-Hugoniot Factors
def _normal_reflection_factor(op):
    g = _shock_gamma(op)
    n = _mass_density_ratio(op)
    return 2 + ((g + 1) * (n - 1)) / 2

def _peak_particle_mach_number(pfree):
    n = _mass_density_ratio(pfree)
    return ((pfree * (1 - (1 /n))) / 142000)**0.5

def _shock_front_mach_number(pfree):
    n = _mass_density_ratio(pfree)
    vc = _peak_particle_mach_number(pfree)
    return vc / (1 - 1 / n)

def _scale_slant_range(r, y, alt):
    sgr = r / y**(1.0/3)
    shob = alt / y**(1.0/3)
    return np.sqrt(sgr**2 + shob**2)

def _regular_mach_merge_angle(r, y, alt):
    pfree = _DNA1kTfreeairop(_scale_slant_range(r, y, alt))
    t = 340 / pfree**0.55
    u = 1 / (7782 * pfree**0.7 + 0.9)
    return np.arctan(1 / (t + u))

def _merge_region_width(r, y, alt):
    pfree = _DNA1kTfreeairop(_scale_slant_range(r, y, alt))
    t = 340 / pfree**0.55
    w = 1 / (7473 * pfree**0.5 + 6.6)
    v = 1 / (647 * pfree**0.8 + w)
    return np.arctan(1 / (t + v))

def _regular_mach_switching_parameter(r, y, alt):
    sgr = r / y**(1.0/3)
    shob = alt / y**(1.0/3)
    alpha = np.arctan(shob / sgr)
    s = (alpha - _regular_mach_merge_angle(r, y, alt)) / _merge_region_width(r, y, alt)
    s0 = max(min (s, 1), -1)
    return 0.5 * (np.sin(0.5 * np.pi * s0) + 1)

def _p_mach(r, y, alt):
    sgr = r / y**(1.0/3)
    shob = alt / y**(1.0/3)
    alpha = np.arctan(shob / sgr)
    a = min(3.7 - 0.94 * np.log(sgr), 0.7)
    b = 0.77 * np.log(sgr) - 3.8 - 18 / sgr
    c = max(a, b)
    return _DNA1kTfreeairop(sgr / 2**(1.0/3)) / (1 - c * np.sin(alpha))

def _p_reg(r, y, alt):
    sgr = r / y**(1.0/3)
    pfree = _DNA1kTfreeairop(_scale_slant_range(r, y, alt))
    shob = alt / y**(1.0/3)
    alpha = np.arctan(shob / sgr)
    r_n = 2 + ((_shock_gamma(pfree) + 1) * (_mass_density_ratio(pfree) - 1)) / 2
    f = pfree / 75842
    d = (f**6 * (1.2 + 0.07 * f**0.5) ) / (f**6 + 1)
    return pfree * ((r_n - 2) * np.sin(alpha)**d + 2)

# PAIR
def _DNAairburstpeakop(r, y, alt):
    sigma = _regular_mach_switching_parameter(r, y, alt)
    if sigma == 0:
        return _p_mach(r, y, alt)
    elif 0 < sigma < 1:
        return _p_reg(r, y, alt) * sigma + _p_mach(r, y, alt) * (1 - sigma)
    elif sigma == 1:
        return _p_reg(r, y, alt)

#QAIR
def _DNAairburstpeakdyn(r, y, alt):
    pair = _DNAairburstpeakop(r, y, alt)
    sigma = _regular_mach_switching_parameter(r, y, alt)
    sgr = r / y**(1.0/3)
    shob = alt / y**(1.0/3)
    alpha = np.arctan(shob / sgr)
    n_q = _mass_density_ratio(pair)
    return 0.5 * pair * (n_q - 1) * (1 - (sigma * np.sin(alpha)**2))

# Airburst blast wave time-of-arrival
def _scaledmachstemformationrange(y, alt):
    shob = alt / y**(1.0/3)
    return shob**2.5 / 5822 + 2.09 * shob**0.75

def _slantrangescalingfactor(r, y, alt):
    sgr = r / y**(1.0/3)
    x_m = _scaledmachstemformationrange(y, alt)
    if sgr <= x_m:
        return 1
    else:
        return 1.26 - 0.26 * (x_m / sgr)

# TAAIR
def _airburstblastwavetoa(r, y, alt):
    v = _slantrangescalingfactor(r, y, alt)
    r1 = _scale_slant_range(r, y, alt) / v
    ta_air = _scaledfreeairblastwavetoa(r1)
    return ta_air * y**(1.0/3) * v

# Overpressure total impulse

def _scaledopposphasedursurf(r, y, alt):
    v = _slantrangescalingfactor(r, y, alt)
    r1 = _scale_slant_range(r, y, alt) / v
    ta_air = _scaledfreeairblastwavetoa(r1)
    t_0 = np.log(1000 * ta_air) / 3.77
    return (155 * np.exp(-20.8 * ta_air) + np.exp(-t_0**2 + 4.86 * t_0 + 0.25)) / 1000

def _scaledopposphasedur(r, y, alt):
    shob = alt / y**(1.0/3)
    v = _slantrangescalingfactor(r, y, alt)
    r1 = _scale_slant_range(r, y, alt) / v
    ta_air = _scaledfreeairblastwavetoa(r1) * v
    t_0 = np.log(1000 * ta_air) / 3.77
    dp_surf = (155 * np.exp(-20.8 * ta_air) + np.exp(-t_0**2 + 4.86 * t_0 + 0.25)) / 1000
    dp_unmod = dp_surf * (1 - (1 - 1 / (1 + 4.5e-8 * shob**7)) * (0.04 + 0.61 / (1 + ta_air**1.5 / 0.027)))
    return dp_unmod * 1.16 * np.exp(-abs(shob / 0.3048 - 156) / 1062)

# Functions for double-peak overpressure waveforms

def _peaksequalityapprox(shob):
    return 138.3 / (1 + 45.5 / shob)

def _peakstimeseperationapprox(shob, sgr, x_m):
    return max(shob / 8.186e5 * (sgr - x_m)**1.25, 1e-12)

def _DNA_b(sgr, shob, ta_air, dp, t):
    """A lot of physically meaningless internal functions used in the DNA overpressure and dynamic pressure calculations."""
    s = 1 - 1 / (1 + (1 / 4.5e-8 * shob**7)) - ((5.958e-3 * shob**2) / (1 + 3.682e-7 * shob**7)) / (1 + sgr**10 / 3.052e14)
    f = s * ((2.627 * ta_air**0.75) / (1 + 5.836 * ta_air) + (2341 * ta_air**2.5) / (1 + 2.541e6 * ta_air**4.75  - 0.216)) + 0.7076 - 3.077 / (1e-4 * ta_air**-3 + 4.367)
    g = 10 + s * (77.58 - 154 * ta_air**0.125 / (1 + 1.375 * np.sqrt(ta_air)))
    h = s * ((17.69 * ta_air) / (1 + 1803 * ta_air**4.25) - (180.5 * ta_air**1.25) / (1 + 99140 * ta_air**4) - 1.6) + 2.753 + 56 * ta_air /(1 + 1.473e6 * ta_air**5)
    return (f * (ta_air / t)**g + (1 - f) * (ta_air / t)**h) * (1 - (t - ta_air) / dp)

def _opatscaledtime(r, y, alt, sgr, shob, x_m, ta_air, dp, t):
    
    s = 1 - 1 / (1 + (1 / 4.5e-8 * shob**7)) - ((5.958e-3 * shob**2) / (1 + 3.682e-7 * shob**7)) / (1 + sgr**10 / 3.052e14)
    f = s * ((2.627 * ta_air**0.75) / (1 + 5.836 * ta_air) + (2341 * ta_air**2.5) / (1 + 2.541e6 * ta_air**4.75  - 0.216)) + 0.7076 - 3.077 / (1e-4 * ta_air**-3 + 4.367)
    g = 10 + s * (77.58 - 154 * ta_air**0.125 / (1 + 1.375 * np.sqrt(ta_air)))
    h = s * ((17.69 * ta_air) / (1 + 1803 * ta_air**4.25) - (180.5 * ta_air**1.25) / (1 + 99140 * ta_air**4) - 1.6) + 2.753 + 56 * ta_air /(1 + 1.473e6 * ta_air**5)
    b = (f * (ta_air / t)**g + (1 - f) * (ta_air / t)**h) * (1 - (t - ta_air) / dp)
    if x_m > sgr or shob > 116:
        return _DNAairburstpeakop(r, y, alt) * b
    else:
        x_e = _peaksequalityapprox(shob)
        e = max(min(abs((sgr - x_m) / (x_e - sgr)), 50), 0.02)
        w = 0.583 / (1 + 2477 / shob**2)
        d = 0.23 + w + 0.27 * e + e**5 * (0.5 - w)
        a = (d - 1) * (1 - 1 / (1 + e**-20))
        dt = _peakstimeseperationapprox(shob, sgr, x_m)
        v_0 = shob**6 / (2445 * (1 + shob**6.75 / 3.9e4) * (1 + 9.23 * e**2))
        c_0 = (1.04 - 1.04 / (1 + 3.725e7 / sgr**4)) / ((a + 1) * (1 + 9.872e8 / shob**9))
        g_a = max(min((t - ta_air) / dt, 400), 0.0001)
        v = 1 + v_0 * g_a**3 / (g_a**3 + 6.13)
        c = c_0 * (1 / (g_a**-7 + 0.923 * g_a**1.5)) * (1 - ((t - ta_air) / dp)**8)
        return _DNAairburstpeakop(r, y, alt) * (1 + a) * (b * v + c)

def _overpressureatscaledtime(r, y, alt, t):
    sgr = r / y**(1.0/3)
    shob = alt / y**(1.0/3)
    x_m = _scaledmachstemformationrange(y, alt)
    v = _slantrangescalingfactor(r, y, alt)
    r1 = _scale_slant_range(r, y, alt) / v
    ta_air = _scaledfreeairblastwavetoa(r1) * v
    dp = _scaledopposphasedur(r, y, alt)
    ss = _opatscaledtime(r, y, alt, sgr, shob, x_m, ta_air, dp, t)
    return np.append( np.zeros(np.sum(t <= ta_air)) , ss[-np.sum(t > ta_air)::])

# In lieu of the 20-point Gauss-Legendre quadrature used in the original
# BLAST.EXE, this fuction uses scipy.integrate.quad to call the venerable FORTRAN
# library QUADPACK and perform a Gauss-Kronod quadrature. This appears
# to be more accurate than the BLAST.EXE help file claims for the original
# approach, which is not surprising as it uses an adaptive algorithm that
# attempts to reduce error to within a particular tolerance.

# IPTOTAL
def _overpressuretotalimpulse(r, y, alt):
    sgr = r / y**(1.0/3)
    shob = alt / y**(1.0/3)
    x_m = _scaledmachstemformationrange(y, alt)
    v = _slantrangescalingfactor(r, y, alt)
    r1 = _scale_slant_range(r, y, alt) / v
    ta_air = _scaledfreeairblastwavetoa(r1) * v
    dp = _scaledopposphasedur(r, y, alt)
    t_p = 13 * (ta_air + dp) / 14
    scaled_impulse, _ = quad(lambda t: _opatscaledtime(r, y, alt, sgr, shob, x_m, ta_air, dp, t), ta_air, ta_air + dp)
    return scaled_impulse * y**(1.0/3)

def _dpt(shob, sgr, x_m, ta_air, dp_q, pair, b, t):
    if x_m > sgr or shob > 116:
        return pair * b
    else:
        dt = _peakstimeseperationapprox(shob, sgr, x_m)
        g_a = max(min((t - ta_air) / dt, 400), 0.0001)
        x_e = _peaksequalityapprox(shob)
        e = max(min(abs((sgr - x_m) / (x_e - sgr)), 50), 0.02)
        w = 0.583 / (1 + 2477 / shob**2)
        d = 0.23 + w + 0.27 * e + e**5 * (0.5 - w)
        a = (d - 1) * (1 - 1 / (1 + e**-20))
        v_0 = shob**6 / (2445 * (1 + shob**6.75 / 3.9e4) * (1 + 9.23 * e**2))
        c_0 = (1.04 - 1.04 / (1 + 3.725e7 / sgr**4)) / ((a + 1) * (1 + 9.872e8 / shob**9))
        v = 1 + v_0 * g_a**3 / (g_a**3 + 6.13)
        c = c_0 * (1 / (g_a**-7 + 0.923 * g_a**1.5)) * (1 - ((t- ta_air) / dp_q)**8)
        return pair * (1 + a) * (b * v + c)
        
# Dynamic pressure total impulse

def _dynamicpressureatscaledtime(r, y, alt, t):
    pair = _DNAairburstpeakop(r, y, alt)
    sgr = r / y**(1.0/3)
    shob = alt / y**(1.0/3)
    x_m = _scaledmachstemformationrange(y, alt)
    v = _slantrangescalingfactor(r, y, alt)
    sr = _scale_slant_range(r, y, alt)
    ta_air = _scaledfreeairblastwavetoa(sr / v) * v
    shob_0 = shob / 0.3048
    sgr_0 = sgr / 0.3048
    shob_x = abs(shob_0 - 200) + 200
    sgr_x = sgr_0 - 200
    dp_0 = 0.3 + 0.42 * np.exp(-shob_x / 131)
    if sgr_x > 0:
        dp_x = dp_0 + 4.4e-5 * sgr_x
    else:
        dp_x = dp_0 + sgr_x * (1.0 / 2361 - abs(shob_x -533)**2 / 7.88e7)
    if shob_0 >= 200:
        dp_q = dp_x
    else:
        dp_q = dp_x * (1 + 0.2 * np.sin(shob_0 * np.pi /200))
    delta_0 = max(shob_0**1.52 / 16330 - 0.29, 0)
    delta = 2.38 * np.exp(-7e-7 * abs(shob_0 - 750)**2.7 - 4e-7 * sgr_0**2) + delta_0
    b = _DNA_b(sgr, shob, ta_air, dp_q, t)
    dpt = _dpt(shob, sgr, x_m, ta_air, dp_q, pair, b, t)
    dpt = np.abs(dpt)
    n_q = _mass_density_ratio(dpt)
    ss = 0.5 * dpt * (n_q - 1) * (dpt / pair)**delta
    return np.append( np.zeros(np.sum(t <= ta_air)) , ss[-np.sum(t > ta_air)::])

def DNA_static_overpressure(y, r, h, yunits='kT', dunits='m', opunits='kg/cm^2'):
    """Estimate peak static overpressure at range r from a burst of yield y using the
the Defense Nuclear Agency 1kT standard free airburst overpressure. This assumes a
thermally ideal surface."""
    yld = convert_units(y, yunits, 'kT')
    gr = convert_units(r, dunits, 'm')
    height = convert_units(h, dunits, 'm')
    op = _DNAairburstpeakop(gr, yld, height)
    return convert_units(op, 'Pa', opunits)

def DNA_dynamic_pressure(y, r, h, yunits='kT', dunits='m', opunits='kg/cm^2'):
    """Estimate peak pynamic overpressure at range r from a burst of yield y using the
the Defense Nuclear Agency 1kT standard free airburst overpressure, assuming an ideal
surface. Many real-world surfaces are not ideal (most, in the opinion of Soviet 
analysts), meaning that this function has only limited predictove capability."""
    yld = convert_units(y, yunits, 'kT')
    gr = convert_units(r, dunits, 'm')
    height = convert_units(h, dunits, 'm')
    dyn = _DNAairburstpeakdyn(gr, yld, height)
    return convert_units(dyn, 'Pa', opunits)




'''the following part was modified from examples.interactive_overpressure.py'''
plt.style.use('dark_background')
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)


t = np.arange(0.01, 10, 0.0025)

r0 = 500
y0 = 20
h0 = 500

s1 = convert_units(_overpressureatscaledtime(r0,y0,h0,t), 'Pa', 'kg/cm^2')
s2 = convert_units(_dynamicpressureatscaledtime(r0,y0,h0,t), 'Pa', 'kg/cm^2')
l1, = plt.plot(t, s1, lw=2, color='magenta', label='static overpressure')

l2, = plt.plot(t, s2, lw=2, color='blue', label='dynamic pressure')
plt.axis([0, 10, -5, 15])
plt.grid()

ax.set_title('interactive overpressure and dynamic pressure calculator')

ax.set_xlabel('time ($s$)')
ax.set_ylabel('pressure ($kg/cm^2$)')
ax.legend()


axyield = plt.axes([0.25, 0.09, 0.65, 0.03])
axdistance = plt.axes([0.25, 0.05, 0.65, 0.03])
axheight = plt.axes([0.25, 0.01, 0.65, 0.03])

syield = Slider(axyield, 'Yield ($kT$):', 1.0, 1000.0, valinit=y0)
sdistance = Slider(axdistance, 'Distance ($m$):', 0.1, 2000.0, valinit=r0)
sheight = Slider(axheight, 'Burst altitude ($m$):', 0.1, 2000.0, valinit=h0)

def update(val):
    burst_distance = sdistance.val
    burst_height = sheight.val
    bomb_yield = syield.val
    l1.set_ydata(convert_units(_overpressureatscaledtime(burst_distance,bomb_yield,burst_height,t), 'Pa', 'kg/cm^2'))
    l2.set_ydata(convert_units(_dynamicpressureatscaledtime(burst_distance,bomb_yield,burst_height,t), 'Pa', 'kg/cm^2'))
    fig.canvas.draw_idle()

syield.on_changed(update)
sdistance.on_changed(update)
sheight.on_changed(update)

cursor = Cursor(ax, useblit=False, color='red', linewidth=2 )
plt.show()
