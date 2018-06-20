
import numpy 
import random

def exp(x,tauinv) :
    return tauinv*numpy.exp(-x*tauinv)



def aleatoire(x) :
    liste=[]
    [liste.append(sum(x[:i+1])) for i in range(0,6)]
    nb=random.random()
    liste.append(nb)
    inde=list(sorted(liste)).index(nb)
    return inde 