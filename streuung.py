# Molecular scattering via quadrature, Beddard pp 611-617
# Scattering angles with Lennard-Jones potential
from numpy import *
from scipy.optimize import fsolve, bisect
from scipy.integrate import quad
from matplotlib.pyplot import *

E0 = 0.2
sig = 1.0
n = 550  # how many b values
bm = 3.0  # max b value

U = lambda r, eps: 4*eps*((sig/r)**12 - (sig/r)**6)  # LJ potential
dU = lambda r, eps: 24*eps*(-2*(sig/r)**13 + (sig/r)**7)
g = lambda r, b, eps: b/(r*r*sqrt(1.0 - U(r, eps)/E0 - (b/r)**2))

""""
def F(a,b,c):
    return (
        lambda x: (x-a)*(x-b)*(x-c),
        lambda x: a*b + a*c + b*c - 2*(a + b + c)*x + 3*x**2,
        lambda x: -2*(a + b + c - 3*x)
    )
"""

def find_largest_root(f, fp, fpp, a=1e-2, b=1e2):
    """
    f: f(x)
    fp: df(x)/dx
    fp: d^2f(x)/dx^2
    a, b: Intervallgrenzen
    """
    r = 1.0
    #######################################################
    #                                                     #
    # TODO: Implementieren Sie hier die Nullstellen Suche #
    #                                                     #
    #######################################################
    r = []
    
    try:
        zeta = bisect(fpp, a, b)
        
        try:
            mu1 = bisect(fp, a, zeta)
        except:
            r.append(bisect(f, a, zeta))
        try:
            mu2 = bisect(fp, zeta, b)
        except:
            r.append(bisect(f, zeta, b))
        
        try:
            r.append(bisect(f, a, mu1))
        except:
            pass
        try:
            r.append(bisect(f, mu1, m2))
        except:
            pass
        try:
            r.append(bisect(f, mu2, b))
        except:
            pass
    except:
        try:
            mu = bisect(fp, a, b)
            try:
                r.append(bisect(f, a, mu))
            except:
                pass
            try:
                r.append(bisect(f, a, mu))
            except:
                pass
        except:
            try:
                r.append(bisect(f, a, b))
            except:
                pass
    
    #if len(r) == 0:
    #    return bisect(f, a, b)
        
    return max(r)


########################
#                      #
# TODO: Select epsilon #
#                      #
########################
#eps_range = array([0.01,1.0])
#eps_range = linspace(0.01, 2, 20)
eps_range = linspace(0.01, 2, 20)

# Radien und Winkel fuer alle eps und b
# Ein array pro fixem eps Wert
# Siehe auch Plots unten.
r0E0 = []
angle_eps = []

for eps in eps_range:
    bval = linspace(0.001, bm, n)
    # Radien und Winkel fuer aktuelles eps
    r0s = []
    angle = []

    for b in bval:
        ur0 = lambda r0: r0*r0*(1 - U(r0, eps)/E0) - b**2
        dur0 = lambda r0: 2*r0*(1 - U(r0, eps)/E0) - r0*r0*dU(r0, eps)/E0
        d2ur0 = lambda r: 2+(40*(-11 + 2*r**6)*eps)/(E0*r**12)
        ##############################################################
        #                                                            #
        # TODO: Finden Sie die groesste Nullstelle und berechnen Sie #
        #       den minimalen Abstand r0 und Streuwinkel theta       #
        #       fuer das gegebene eps und b.                         #
        #                                                            #
        ##############################################################
        #r0 = fsolve(ur0, 10, fprime=dur0)
        r0 = find_largest_root(ur0, dur0, d2ur0)
        dthetadr = lambda r: g(r, b, eps)
        theta = pi - 2. * quad(dthetadr, r0, np.inf)[0]
        # TODO: try here or correct b interval?
        #try:
        #    theta = pi - 2. * quad(dthetadr, r0, np.inf)[0]
        #except ZeroDivisionError:
        #    theta = float('nan')
        #    print("encountered division by zero")    
        r0s.append(r0)
        angle.append(theta)

    r0E0.append(array(r0s))
    angle_eps.append(array(angle))


figure(figsize=(10, 10))
for eps, angle in zip(eps_range, angle_eps):
    plot(bval, angle, label=r'$\varepsilon = %.3f$' % eps)
ylim(-2*pi, pi)
xlabel('$b$')
ylabel(r'$\theta_0$')
legend(loc='lower left')
grid(True)
savefig('scattering_angle.pdf')

figure(figsize=(10, 10))
for eps, r0 in zip(eps_range, r0E0):
    plot(bval, r0, label=r'$\varepsilon = %.3f$' % eps)
grid(True)
xlabel('$b$')
ylabel('$r_0$')
legend(loc='upper left')
savefig('minimal_distance.pdf')
