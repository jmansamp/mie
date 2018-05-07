import numpy as np
from lmfit import minimize, Parameters, Model
import mg 
import matplotlib.pylab as plt



#we need to define a function, normalize data 

data = np.loadtxt("uv_samp.tsv")

x = []
y = []

spr =0
sprabs = 0

for line in data:
    if line[0]%5==0 and line[0]>= 200:
        y.append(line[1])
        x.append(int(line[0]))
    if line[1] > sprabs and line[0]>450:
        spr = line[0]
        sprabs = line[1]

q_list = np.linspace(0,130,131, dtype = int)



miemodel = Model(mg.absorbtot)

params = miemodel.make_params(Ns = 4e14, Ne = 1e5, B = 6542.15582,
        R = 7.9977e-9, vary=False, Sg = 0.5)

params['Ns'].set(vary = False)
params['B'].set(vary = False)
params['R'].set(vary = False)

result = miemodel.fit(y,params, q=q_list)

print(result.fit_report())

plt.plot(x,y, 'bo')
#plt.plot(x, result.init_fit, 'k--')
plt.plot(x, result.best_fit, 'r-')
plt.show()


        
