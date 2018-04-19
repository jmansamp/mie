


def mie_extinction(w,ei):
    e = np.sqrt(1-(b/a)**2 
    Pa = (1-e**2)/e**2*[1/(2*e)*np.ln((1+e)/(1-e))-1]
    Pb = (1-Pa)/2
    Pc = (1-Pa)/2
    V = (4*np.pi/3)*a*b**2

    E1 = E + w_p**2*((w**2+gammabulk**2)**-1-(w**2+gamma**2)**-1)
    
    E2 = Ei + w_p**2/w*(gamma/(w**2+gamma**2)-gammabulk/(w**2+gammabulk**2)
    return 2*np.p*V*Em**(3/2.)/(3*lmda)* E2/Pa**2/((E1+(1-Pa)/Pa*Em)**2 + E2**2)
    + E2/Pb**2/((E1+(1-Pb)/Pa*Em)**2 + E2**2) + E2/Pc**2/((E1+(1-Pc)/Pc*Em)**2 + E2**2) 
