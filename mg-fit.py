import scipy.special as sp
import numpy as np
import matplotlib.pylab as plt
from scipy.misc import derivative
from scipy.optimize import curve_fit

uv_data = np.loadtxt('uv_samp.tsv')

# set matrice/solvent refraction index
nmatrice = 1.334

L = 3

xdata = []
ydata = []

spr = 0
sprabs = 0

for line in uv_data:
    if line[0]%5==0 and line[0]>= 400:
        ydata.append(line[1])
        xdata.append(line[0])
    if line[1] > sprabs and line[0]>450:
        spr = line[0]
        sprabs = line[1]


        
ydata = np.array(ydata)/sprabs
xdata = np.array(xdata)
# set the A parameter coeficients

A1 = 1.08883
A2 = 0.38984
A3 = -0.03956
A4 = 0

n_m = 10
c = 299792458 #speed of light

#input gold dielectric values
Epsilon_Bulk = np.loadtxt('J&CAu.txt')

#input normalized data
data = np.loadtxt('Normalized Sample.Dat')
    
gammabulk = 1/2.94e-14
v_F = 1.4e6
ne = 5.9e28
mEFF = 1.1*9.1093897e-31
e = 1.60217733e-19
E_0 = 8.854e-12
w_p = np.sqrt(ne*e**2/(E_0*mEFF))


#Ricatti-Bessel functions of first and second kind
def nu_L(L,v):
    #return sp.riccati_jn(L,v)
    nu_list = []
    for l in range(1,L+1):

        nu_list.append(v*np.sqrt(np.pi/(2*v))*(sp.jv(l+0.5,v) + 1j*sp.yv(l+0.5,v)))   
    return np.array(nu_list)


def dif(x, fn,L):
    h = x/float(1000)
    """This function takes the arguments:
    x: The point at which the gradient needs to be found.
    fn: The function to be differentiated.
    h: The step size of the function which defaults to h = 0.0000000001."""
    
    return (fn(L,x + h) - fn(L,x))/h


def psi_L(L,v):
    
    #return sp.riccati_yn(L,v)
    psi_list = []
    
    for l in range(1,L+1):
        psi_list.append(v*np.sqrt(np.pi/(2*v))*(sp.jv(l+0.5,v))) 
    
    return np.array(psi_list)

def difnu_L(L,v):
    return dif(v,nu_L,L)

def difpsi_L(L,v):
    return dif(v,psi_L,L)

def a_n(L,m,x):
    return -(m*psi_L(L,m*x)*difpsi_L(L,x)-difpsi_L(L,m*x)*psi_L(L,x))/(m*nu_L(L,m*x)*difpsi_L(L,x)- difnu_L(L,m*x)*psi_L(L,x))

def b_n(L,m,x):
     return -(psi_L(L,m*x)*difpsi_L(L,x)-m*difpsi_L(L,m*x)*psi_L(L,x))/(nu_L(L,m*x)*difpsi_L(L,x)- m*difnu_L(L,m*x)*psi_L(L,x))



def gammareduced(R):
    A = A1 + A2*R*1e9 + A3*(R*1e9)**2 + A4*(R*1e9)**3 
    
    if A<0: 
        gamma = gammabulk 
    else:
        gamma = gammabulk + A*v_F/R

    return gamma

def extinction(q,R):

    gamma = gammareduced(R)
    
    lmda = Epsilon_Bulk[:,0][q]
    
    E = Epsilon_Bulk[:,1][q]
    
    Ei = Epsilon_Bulk[:,2][q]
    
    k = 2*np.pi/(lmda*1e-9)

    w = c*k

    enne = np.sqrt(E + w_p**2*((w**2+gammabulk**2)**-1-(w**2+gamma**2)**-1)
            + 1j*(Ei + w_p**2/w*(gamma/(w**2+gamma**2)-gammabulk/(w**2+gammabulk**2))))

    x = k*R*enne  
    
    m = nmatrice/enne

    a = a_n(3,m,x)

    b = b_n(3,m,x)
    
    asum = 0
    bsum = 0

    for i in range(0,3):
        asum += (2*(i+1)+1)*a[i]
        bsum += (2*(i+1)+1)*b[i]
    
    return -2*np.pi/(nmatrice*k)**2*(np.real(asum + bsum))


def absorb(q,N,R,B):
    d = 10e-3
    return -B*np.log10(np.exp(-1*d*N*extinction(q,R)))

#popt, pcov = curve_fit(absorb, np.linspace(40,130,91,dtype = int),
#    ydata,p0=[1E10,10e-9,0.5], bounds = (0,[1e20,50e-9,1000000]))

#plt.plot(xdata,absorb(np.linspace(40,130,91,dtype
#    = int),popt[0],popt[1],popt[2]), label = 'Fit')
#plt.plot(xdata,ydata, label = 'data')
#plt.legend(loc=0)
#plt.show()

def mie_extinction(q,a,b):
    e = np.sqrt(1.-(float(b)/a)**2.)
    E = Epsilon_Bulk[:,1][q]
    Ei = Epsilon_Bulk[:,2][q]
    lmda = np.array(Epsilon_Bulk[:,0][q])
    print(lmda)
    Pa = (1.-e**2)/e**2*(1./(2*e)*np.log((1.+e)/(1.-e))-1.)
    Pb = (1-Pa)/2

    V = (4*np.pi/3)*a*b**2
    k = 2*np.pi/(lmda*1e-9)
    w = c*k
    
    R = (a*b**2)**(1/3)

    gamma = gammareduced(R)

    Em = nmatrice**2

    E1 = E #+ w_p**2*((w**2+gammabulk**2)**(-1)-(w**2+gamma**2)**(-1))

    E2 = Ei #+ w_p**2/w*(gamma/(w**2+gamma**2)-gammabulk/(w**2+gammabulk**2))
    

    return 2*np.pi*V*Em**(3/2.)/(3*lmda*1e-9)*(E2/Pa**2/((E1+(1-Pa)/Pa*Em)**2
        + E2**2) + 2*E2/Pb**2/((E1+(1-Pb)/Pb*Em)**2 + E2**2))

lmda_list = np.linspace(0,160,161, dtype = int)

xarray = Epsilon_Bulk[:,0][lmda_list]

def absorbmg(q,N,a,b,B):
    d = 1e-3
    return -B*np.log10(np.exp(-1*d*N*mie_extinction(q,a,b)))

def axis(R,p):
    v = R**3
    b = (v/p)**(1/3)
    a = (p*b)
    return a,b
    

# Demonstrating exctinction with varying ar
mie_sp = absorbmg(lmda_list,10000000000,50e-9,10e-9,1)

mie_g1 = mie_extinction(lmda_list,*axis(5e-9,1.2))

mie_g2 = mie_extinction(lmda_list,*axis(5e-9,2))

mie_g3 = mie_extinction(lmda_list,15e-9,5e-9)

plt.plot(xarray, mie_sp, label = 'Mie particle with ar = 1, V = 125 nm^3')

#plt.plot(xarray,mie_g1, label = 'MG particle with ar = 1.2, V = 125 nm^3')

#plt.plot(xarray,mie_g2, label = 'MG particle with ar = 1.5, V = 125 nm^3')

#plt.plot(xarray,mie_g3, label = 'MG particle with ar = 2, V = 125 nm^3')

plt.legend(loc=0)

plt.show()


# Demonstrating exctinction with varying size




def gauss(a,b,Sg):
    return 1/(Sg*np.sqrt(2*np.pi))*np.exp(-(a/b-1)**2/(s*Sg**2))

def total(q,a,b,Sg, Ns, Ne):
    return 1
