import numpy 
import matplotlib
matplotlib.use('Agg')
import random  
import function 
from variable import active_time
from variable import activation_time
from variable import end_time


def generatif(trans, transbase , OMEGA ,probabehavior,larva,version):

    changement=[[],[],[],[],[],[]]
    probaS=numpy.zeros((6,end_time))
    gene=numpy.zeros((larva , end_time+1))
    comp=[]
    compbase=[]
    temps_comp=[]
    [comp.append([]) for i in range(0,larva)]
    [compbase.append([]) for i in range(0,larva)]
    [temps_comp.append([]) for i in range(0,larva)]
    
    if version == "version1":


        for i in range(0,larva):
            Behav=function.aleatoire(probabehavior)
            compbase[i].extend([[Behav,0.0]])
            deltatemps=0.0
            for j in range(0,activation_time):
                DT=0
                old=Behav
                gene[i,j]=Behav
                nbb=random.random()
                if nbb>(numpy.exp(-OMEGA[0][Behav])):
                    dt=-numpy.log(nbb)/OMEGA[0][Behav]
                    compbase[i][-1][1]+=dt
                    Behav=function.aleatoire(transbase[Behav])
                    compbase[i].extend([[Behav,1-dt]])    
                else :
                    compbase[i][-1][1]+=1   
            del compbase[i][-1]
            comp[i].extend([[Behav,0.0]])
            for j in range(1,active_time):
                old=Behav
                gene[i,j+activation_time-1]=Behav
                nbb=random.random()
                if nbb>(numpy.exp(-OMEGA[j][Behav])):
                    dt=-numpy.log(nbb)/OMEGA[j][Behav]
                    comp[i][-1][1]+=dt
                    a=j
                    while sum(trans[j][Behav])==0 :
                        trans[j][Behav]=trans[a][Behav]
                        a-=1
                    Behav=function.aleatoire(trans[j][Behav])
                    comp[i].extend([[Behav,1-dt]])  
                else :
                    comp[i][-1][1]+=1
            del comp[i][-1]     




    if version == "version2":

        for i in range(0,larva):

            Behav=function.aleatoire(probabehavior)
            #compbase[i].extend([[Behav,0.0]])
 
            time=0
            time=0
            while time <activation_time :
                if time>0 :
                    Behav=function.aleatoire(transbase[Behav])
                old=time
                nbb=random.random()
                dt=-numpy.log(nbb)/OMEGA[0][Behav]
                compbase[i].extend([[Behav,dt]]) 
                #print(Behav , dt)
                time=time+dt
                if time>activation_time :
                    del compbase[i][-1]
                    compbase[i].extend([[Behav,activation_time-old]]) 
                    time=activation_time
                for e in range(int(old),int(time)-1): 
                    #print(time,e) 
                    gene[i,e]=Behav 
            #del compbase[i][-1]
            time=1
            #print(trans[1])
            while time<active_time:
                #print(Behav,j,len(trans),sum(trans[j][Behav]))
                old=time
                #if j>0: 
                a=int(time)
                b=int(time)
                A=trans[int(time)][Behav]
                #print(A,trans[int(time)][Behav])
                while sum(A)==0 :
                    a-=1
                    if a<0:
                        a=0
                    b+=1
                    if b>15:
                        b=15
                    A=[]
                    [A.append((trans[a][Behav][k]+trans[b][Behav][k])/2) for k in range(0,6)]
                SUM=sum(A)
                A = [k/SUM for k in A]
                #print(A)
                Behav=function.aleatoire(A)   
                nbb=random.random()
                dt=-numpy.log(nbb)/OMEGA[int(time)][Behav]
                comp[i].extend([[Behav,dt]])
                time=time+dt
                if time>active_time :
                    time=active_time
                for e in range(int(old),int(time)-1): 
                    gene[i,e+activation_time-1]=Behav 
               
        print('version2')
        #print(comp)



    if version== "version3":

        for i in range(0,larva):

            Behav=function.aleatoire(probabehavior)
            compbase[i].extend([[Behav,0.0]])
            temps_comp[i].extend([[Behav,0.0]])

            for j in range(0,activation_time):
                    Dt=1.0
                    gene[i,j]=Behav
                    nbb=random.random()
                    while nbb>(numpy.exp(-OMEGA[0][Behav]*Dt)):
                        dt=-numpy.log(nbb)/OMEGA[0][Behav]
                        compbase[i][-1][1]+=dt
                        Dt-=dt
                        Behav=function.aleatoire(transbase[Behav])
                        compbase[i].extend([[Behav,0]])  
                        temps_comp[i].extend([[Behav,j]])
                        nbb=random.random()  
                    compbase[i][-1][1]+=Dt
            
            comp[i].extend([[Behav,0]])
            #del compbase[i][-1]

            for j in range(1,active_time+1):
                Dt=1
                old=Behav
                gene[i,j+activation_time-1]=Behav
                nbb=random.random()
                while nbb>(numpy.exp(-OMEGA[j][Behav]*Dt)):
                    dt=-numpy.log(nbb)/OMEGA[j][Behav]
                    Dt-=dt
                    comp[i][-1][1]+=dt
                    a=j
                    b=j
                    A=trans[j][Behav]
                #print(A,trans[int(time)][Behav])
                    # while sum(A)==0 :
                    #     a-=1
                    #     if a<0:
                    #         a=0
                    #     b+=1
                    #     if b>15:
                    #         b=15
                    #     A=[]
                    #     [A.append((trans[a][Behav][k]+trans[b][Behav][k])/2) for k in range(0,6)]
                    SUM=sum(A)
                    A = [k/SUM for k in A]
                    Behav=function.aleatoire(A)
                    comp[i].extend([[Behav,0]])
                    temps_comp[i].extend([[Behav,j+activation_time]])
                    nbb=random.random()
                comp[i][-1][1]+=Dt

    return gene, comp, compbase, temps_comp



    if version== "version4":

        for i in range(0,larva):

            Behav=function.aleatoire(probabehavior)
            compbase[i].extend([[Behav,0.0]])

            for j in range(0,activation_time):
                    Dt=1.0
                    gene[i,j]=Behav
                    nbb=random.random()
                    while nbb>(numpy.exp(-OMEGA[0][Behav]*Dt)):
                        dt=-numpy.log(nbb)/OMEGA[0][Behav]
                        compbase[i][-1][1]+=dt
                        Dt-=dt
                        Behav=function.aleatoire(transbase[Behav])
                        compbase[i].extend([[Behav,0]])  
                        nbb=random.random()  
                    compbase[i][-1][1]+=Dt  
            
            comp[i].extend([[Behav,0]])
            #del compbase[i][-1]

            for j in range(1,active_time+1):
                Dt=1
                old=Behav
                gene[i,j+activation_time-1]=Behav
                nbb=random.random()
                while nbb>(numpy.exp(-OMEGA[j][Behav]*Dt)):
                    dt=-numpy.log(nbb)/OMEGA[j][Behav]
                    Dt-=dt
                    comp[i][-1][1]+=dt
                    a=j
                    b=j
                    A=trans[j][Behav]
                #print(A,trans[int(time)][Behav])
                    while sum(A)==0 :
                        a-=1
                        b+=1
                        A=[]
                        [A.append((trans[a][Behav][k]+trans[b][Behav][k])/2) for k in range(0,6)]
                    SUM=sum(A)
                    A = [k/SUM for k in A]
                    Behav=function.aleatoire(A)
                    comp[i].extend([[Behav,0]])
                    nbb=random.random()
                comp[i][-1][1]+=Dt



