import math
import numpy as np
# Planck 2018, the game's cosmology (R71.10): H0 = 67.4 km/s/Mpc, Omega_m = 0.315, flat, with radiation (9e-5)
# kept here because this script integrates numerically. SRKit's closed forms (srkit-flrw) are matter and Lambda only:
# without radiation the age is 13.80 Gyr, the particle horizon 47.0 Gly and the reachable share 4.5%.
H0=67.4/3.0856775814913673e19; Om,Or=0.315,9e-5; OL=1-Om-Or
Gyr=3.15576e16; yr=3.15576e7; c=299792458.0; Gly=9.4607304725808e24
lna=np.linspace(math.log(1e-8),math.log(1e12),800001)
a=np.exp(lna); Ha=H0*np.sqrt(Or/a**4+Om/a**3+OL); dl=lna[1]-lna[0]
t=np.cumsum(1/Ha)*dl; t_now=np.interp(0,lna,t); tr=t-t_now
eta=np.cumsum(c/(a*Ha))*dl
eta_inf=eta[-1]+c/(a[-1]*Ha[-1])
eta_now=np.interp(0,lna,eta)
chi_eh_now=eta_inf-eta_now; chi_ph=eta_now
print('age %.2f Gyr; particle horizon %.1f Gly; event horizon %.1f Gly'%(t_now/Gyr,chi_ph/Gly,chi_eh_now/Gly))
f=(chi_eh_now/chi_ph)**3
print('reachable share of observable comoving volume: %.1f%%'%(100*f))
# galaxies: ~2e11 (conservative) to 2e12 (Conselice 2016) in observable universe
for N in (2e11,2e12):
    print('  of %.0e observable galaxies: %.1e reachable; lost to horizon per yr now: %.0f'%(N,N*f,N*f*3*c*yr/chi_eh_now))
# launch later: reachable comoving radius
for tl in (1,10,50,100):
    e=np.interp(tl*Gyr,tr,eta); r=(eta_inf-e)
    print('launch at +%d Gyr: reach %.2f Gly comoving (%.1f%% of today\'s reach volume)'%(tl,r/Gly,100*(r/chi_eh_now)**3))
# light arrival time for targets at comoving chi, launched now
for chiG in (0.054,1,5,10,15,16):
    e=eta_now+chiG*Gly
    print('target %5.3g Gly comoving: arrival at +%.3g Gyr'%(chiG,np.interp(e,eta,tr)/Gyr))
