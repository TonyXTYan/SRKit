"""Check: Hawking lifetime of an M_*-scale hole vs alpha_G, in proton-Compton (atomic-like) units and in Planck units.
Also: the hole mass whose lifetime exceeds the U384-ms range in our universe."""
import math
hbar=1.054571817e-34; c=2.99792458e8; mp=1.67262192e-27; Msun=1.3271244e20/6.67430e-11; yr=3.15576e7
aG0=6.67430e-11*mp**2/(hbar*c)
def t_ev(aG, lam=1.0):
    G=aG*hbar*c/mp**2; M=lam*aG**-1.5*mp
    return 5120*math.pi*G**2*M**3/(hbar*c**4)
t0, t1 = t_ev(aG0), t_ev(aG0/10)
tP=lambda aG: math.sqrt(hbar*(aG*hbar*c/mp**2)/c**5)
print('alpha_G0 = %.3e'%aG0)
print('exponent in proton-Compton units: %.3f (expect -2.5)'%(math.log(t1/t0)/math.log(0.1)))
print('exponent in Planck units:         %.3f (expect -3.0)'%(math.log((t1/tP(aG0/10))/(t0/tP(aG0)))/math.log(0.1)))
U384_yr=2**384/1e3/yr
M_over=(U384_yr/(t_ev(aG0, 1)/yr))**(1/3)*(aG0**-1.5*mp)/Msun  # scale mass from M_* reference
# simpler: t = 2.1e67 yr (M/Msun)^3
G=6.674e-11; tsun=5120*math.pi*G**2*Msun**3/(hbar*c**4)/yr
print('t_ev(1 Msun) = %.2e yr; U384-ms range %.2e yr -> overflow above %.2e Msun'%(tsun,U384_yr,(U384_yr/tsun)**(1/3)))
