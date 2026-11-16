import math
from scipy.special import hyp2f1, gamma
from scipy.integrate import quad
c=299792458.0; G=6.67430e-11; Msun=1.3271244e20/G; Gyr=3.15576e16
# Planck 2018 (R71.10): H0 = 67.4 km/s/Mpc, Omega_m = 0.315 give 1/H_Lambda = 17.53 Gyr
HL=1/(17.53*Gyr)
def F(tau):
    W=math.sinh(1.5*tau)**-2
    return W**(1/3)*hyp2f1(0.5,1/3,4/3,-W)
print('F(0+) ~ %.4f, closed form %.4f'%(F(1e-6), gamma(1/6)*gamma(1/3)/(3*math.sqrt(math.pi))))
# independent: matter-only-normalised flat LCDM, a = sinh^(2/3)(1.5 HL t); comoving horizon chi_eh = c * int_t^inf dt/a ; matter mass = rho_m a^3 (4/3 pi chi^3), rho_m a^3 = const
# rho_m a^3 with a=sinh^(2/3): Friedmann H^2 = HL^2 (1 + sinh^-2) => rho_m a^3 = 3 HL^2/(8 pi G) (in units where a^3 sinh^-2 ... ) 
rho_m_a3=3*HL**2/(8*math.pi*G)
def Mreach(t):
    chi=c*quad(lambda s: math.sinh(1.5*HL*s)**(-2/3), t, t+60/HL, limit=200)[0]
    return rho_m_a3*4/3*math.pi*chi**3
for tG in (1,13.79,50):
    t=tG*Gyr
    print('t=%5.2f Gyr: numeric %.3e Msun, closed form %.3e Msun'%(tG, Mreach(t)/Msun, c**3/(2*G*HL)*F(HL*t)**3/Msun))
print('crunch time 2pi/3H: %.1f Gyr'%(2*math.pi/3*17.53))
print('Hawking=dS mass c^3/(4 G H): %.2e Msun'%(c**3/(4*G*HL)/Msun))
