import math
import numpy as np
# Planck 2018, the game's cosmology (R71.10): H0 = 67.4 km/s/Mpc, Omega_m = 0.315, flat
H0=67.4/3.0856775814913673e19; Om,OL=0.315,0.685
Gyr=3.15576e16; c=299792458.0; Gly=9.4607304725808e24
lna=np.linspace(math.log(1e-3),math.log(1e12),400001)
a=np.exp(lna); Ha=H0*np.sqrt(Om/a**3+OL); dl=lna[1]-lna[0]
t=np.cumsum(1/Ha)*dl; t-=np.interp(0,lna,t)
eta=np.cumsum(c/(a*Ha))*dl
eta_inf=eta[-1]+c/(a[-1]*Ha[-1])
chi_eh=eta_inf-eta
print('event horizon now: %.1f Gly'%(np.interp(0,lna,chi_eh)/Gly))
for name,dMpc in [('Virgo cluster',16.5),('Coma cluster',100)]:
    chi=dMpc*3.0856775814913673e22; i=np.argmax(chi_eh<chi)
    print('%s: signals sent later than %.0f Gyr from now never arrive'%(name,t[i]/Gyr))
a1=9.80665; d=16.5*3.0856775814913673e22
print('1 g flip-and-burn to Virgo: tau = %.1f yr'%(2*c/a1*math.acosh(1+a1*d/(2*c*c))/3.15576e7))
kB=1.380649e-23; G=6.67430e-11; hbar=1.054571817e-34; Ms=1.3271244e20/G
S=lambda M:4*math.pi*G*M**2/(hbar*c)
print('S_BH 1 Msun %.1e kB; 4e6 Msun %.1e kB'%(S(Ms),S(4e6*Ms)))
print('Milky Way starlight entropy rate L/T0: %.1e kB/s'%(4e36/2.725/kB))
print('Sun over 10 Gyr into 2.7 K sky: %.1e kB'%(3.828e26*1e10*3.15576e7/2.725/kB))
