import numpy as np
from lmfit import minimize, Parameters, Model
import mg
import matplotlib.pylab as plt



#we need to define a function, normalize data

#set to +-1 for +-5 nm offset

data = np.loadtxt("Sample.dat")

x = []
y = []

spr =0
sprabs = 0

for line in data:
    if line[0]%5==0 and line[0]>= 300 and line[0]<=1000:
        y.append(line[1])
        x.append(int(line[0])
    if line[1] > sprabs and line[0]>450:
        spr = line[0]
        sprabs = line[1]

q_list = np.linspace(20,160,141, dtype = int)[::-1]

y = y

#define fitting model according to MG function
miemodel = Model(mg.absorbtot)

#initialized fitting parameters
params = miemodel.make_params(Ns = 1e12, Ne = 1e12, B = 1.,
        R = 4.2e-9, Sg = 0.75)

#Set constant (fixed) parameters
params['Ns'].set(vary = True)
params['Sg'].set(vary = True)
params['B'].set(vary = False)
params['R'].set(vary = False)

result = miemodel.fit(y,params, q=q_list)

print(result.fit_report())

plt.plot(x,y, 'bo')
#plt.plot(x, result.init_fit, 'k--')
plt.plot(x, result.best_fit, 'r-')
plt.show()



