import numpy 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors  
import random    
from variable import active_time
from variable import activation_time
from variable import end_time



def proba(gene) :

    proba=numpy.zeros((6,end_time))

    for j in range(0,end_time) :
        for i in range(0,6) :
            proba[i,j]=list(gene[:,j]).count(i)/len(gene[:,j])  
    suM=numpy.zeros(end_time)
    for j in range(0,end_time):        
        suM[j]=numpy.sum(proba[:,j])         

    return proba, suM


def matrix2(temps_comp) : 

    matrix2=[]
    [matrix2.append([]) for i in range(0,int(active_time+1))]
    for j in range(0,int(active_time+1)):
        [matrix2[j].append([]) for i in range(0,6)]
        [matrix2[j][i].extend([0,0,0,0,0,0]) for i in range(0,6)]   

    for i in range(0,len(temps_comp)) :
        for j in range(0,len(temps_comp[i])-2) :

            if temps_comp[i][j][1]>activation_time :
                matrix2[int(temps_comp[i][j][1])-activation_time+1][int(temps_comp[i][j][0])][int(temps_comp[i][j-1][0])]+=1
                #print(int(temps_comp[i][j][1])-29)

            else :
                matrix2[0][int(temps_comp[i][j][0])][int(temps_comp[i][j-1][0])]+=1

    for j in range(0,len(matrix2)):   
        for i in range(0,6) :
            SUM=[]
            [SUM.append(sum(matrix2[j][i])) for k in range(0,6)]
            if SUM[0]!=0 :
                matrix2[j][i]=list(matrix2[j][i]/numpy.asarray(SUM))

    return matrix2


def matrix34(comp) : 

    matrix3=numpy.zeros((6,36))
    matrix4=numpy.zeros((36,36))
    for i in range(0,len(comp)) :
        for j in range(0,len(comp[i])-3) :
            matrix3[int(comp[i][j][0]),int(comp[i][j+1][0])*6+int(comp[i][j+2][0])]+=1
            if j<len(comp[i])-4 :
                matrix4[int(comp[i][j][0])*6+int(comp[i][j+1][0]),int(comp[i][j+2][0])*6+int(comp[i][j+3][0])]+=1
    for i in range(0,6):
        #print(i,numpy.sum(matrix3[i,:]))
        if numpy.sum(matrix3[i,:])!=0 :
            matrix3[i,:]=matrix3[i,:]/numpy.sum(matrix3[i,:])  
        if numpy.sum(matrix4[i*6:i*6+6,:])!=0 :
            for j in range(0,6):
                matrix4[i*6+j,:]=matrix4[i*6+j,:]/numpy.sum(matrix4[i*6:i*6+6,:]) 

    return matrix3 , matrix4    

def MatrixCorre(comp,compbase) :
    MC=numpy.zeros((36,36))





def FonctionCorrelationMOD(comp,compbase) :
    Fonction=numpy.zeros((24))
    Fonctionbase=numpy.zeros((40))
    for i in range(0,len(comp)) :
        Correbase=[]
        [Correbase.append(compbase[i][j][1]) for j in range(0,len(compbase[i]))]
        delta=1
        #print(len(Correbase),i,len(compbase[i]))
        while   sum(Correbase[0:delta])<10 :
            for j in range(0,len(Correbase)):
                if int((sum(Correbase[j:j+delta]))*4)<40 :
                    Fonctionbase[int((sum(Correbase[j:j+delta]))*4)]+=1/(len(Correbase)+len(comp))
            delta+=1
        
        Corre=[]
        [Corre.append(comp[i][j][1]) for j in range(0,len(comp[i]))]
        delta=1
        while sum(Corre[0:0+delta])<6  :
            for j in range(0,len(Corre)):
                if int((sum(Corre[j:j+delta]))*4)<24:
                    Fonction[int((sum(Corre[j:j+delta]))*4)]+=1/(len(Corre)+len(comp))
            delta+=1
        #print(i)
    plt.plot(numpy.linspace(0,6,24),Fonction,color= '#60ACB7')
    plt.plot(numpy.linspace(0,10,40),Fonctionbase,color='#e46a4f')
    return    


def distrib_time(comp,compbase) :
    Histo=[]
    for toto in [comp] :
        proba_time=[]
        [proba_time.append([]) for i in range(0,6)]
        for i in range(0,len(toto)):
            [proba_time[int(toto[i][j][0])].append(toto[i][j][1]) for j in range(0,len(toto[i]))]
        Histo.extend([proba_time])

    Histobase=[]
    for toto in [compbase] :
        proba_time=[]
        [proba_time.append([]) for i in range(0,6)]
        for i in range(0,len(toto)):
            [proba_time[int(toto[i][j][0])].append(toto[i][j][1]) for j in range(0,len(toto[i]))]
        Histobase.extend([proba_time])

    return Histo, Histobase 




def CovarianceTemporelle2(comp,compbase):  
    old=0
    lon=0
    for j in range(0,len(compbase)):
        lon+=len(compbase[j])/len(compbase)
        lenght=len(compbase[j])
        if lenght >old:
            old=lenght
    baselim=lon+(lenght-lon)/2

    lon=0
    old=0
    for j in range(0,len(comp)):
        lon+=len(comp[j])/len(comp)
        lenght=len(comp[j])
        if lenght >old:
            old=lenght

    aclim=lon+(lenght-lon)/2
    All2=numpy.zeros((int(aclim)))
    allbi2=numpy.zeros((int(aclim)))
    MOY=numpy.zeros((6))
    MOYt1=numpy.zeros((6))
    nbmoy=numpy.zeros((6))
    Covbase=numpy.zeros((6,int(baselim)))
    Cov=numpy.zeros((6,int(aclim)))
    nbbase=numpy.zeros((6,int(baselim)))
    nb=numpy.zeros((6,int(aclim)))


    for i in range(0,len(compbase)) :
        for j in range(0,len(compbase[i])-1) :
            MOY[int(compbase[i][j][0])]+=compbase[i][j][1]
            nbmoy[int(compbase[i][j][0])]+=1
            MOYt1[int(compbase[i][j][0])]+=compbase[i][j+1][1]

    for i in range(0,len(compbase)) :
        for k in range(0,int(baselim)) :
            for j in range(0,len(compbase[i])-k) :
                Covbase[int(compbase[i][j][0]),k]+=(compbase[i][j][1]-MOY[int(compbase[i][j][0])]/nbmoy[int(compbase[i][j][0])])*(compbase[i][j+k][1]-MOYt1[int(compbase[i][j][0])]/nbmoy[int(compbase[i][j][0])])
                nbbase[int(compbase[i][j][0]),k]+=1
    Covbase=Covbase/nbbase


    MOY=numpy.zeros(6)
    MOYt1=numpy.zeros((6,int(aclim)))
    nbmoy=numpy.zeros(6)
    for i in range(0,len(comp)) :
        for j in range(0,len(comp[i])-1) :
            MOY[int(comp[i][j][0])]+=comp[i][j][1]
            nbmoy[int(comp[i][j][0])]+=1
            for k in range(0,min(int(aclim)-j-1,len(comp[i])-1-j)):
                MOYt1[int(comp[i][j][0]),k]+=comp[i][j+k][1]

    for i in range(0,len(comp)) :
        #print(k,j,len(comp[i])-k)
        for j in range(0,len(comp[i])) :
            for k in range(0,min(int(aclim),len(comp[i])-j)):
                Cov[int(comp[i][j][0]),k]+=(comp[i][j][1]-MOY[int(comp[i][j][0])]/nbmoy[int(comp[i][j][0])])*(comp[i][j+k][1]-MOYt1[int(comp[i][j][0]),k]/nbmoy[int(comp[i][j][0])])
                All2[k]+=(comp[i][j][1]-MOY[int(comp[i][j][0])]/nbmoy[int(comp[i][j][0])])*(comp[i][j+k][1]-MOYt1[int(comp[i][j][0]),k]/nbmoy[int(comp[i][j][0])])
                allbi2[k]+=1
                nb[int(comp[i][j][0]),k]+=1

    Cov=Cov/nb
    numpy.nan_to_num(Cov)
    numpy.nan_to_num(Covbase)

    return Cov, Covbase




def CovarianceTemporelle(LARVA,LARVAbase):  
    old=0
    lon=0
    for j in range(0,len(LARVAbase)):
        lon+=len(LARVAbase[j])/len(LARVAbase)
        lenght=len(LARVAbase[j])
        if lenght >old:
            old=lenght
    baselim=lon+(lenght-lon)/2

    lon=0
    old=0
    for j in range(0,len(LARVA)):
        lon+=len(LARVA[j])/len(LARVA)
        lenght=len(LARVA[j])
        if lenght >old:
            old=lenght

    aclim=lon+(lenght-lon)/2
    All=numpy.zeros((int(aclim)))
    allbi=numpy.zeros((int(aclim)))
    bibu1=numpy.zeros((int(aclim),int(aclim)))
    Cov=numpy.zeros((int(aclim),int(aclim)))
    #print(active_time)
    moy=numpy.zeros((int(aclim)))
    nb=numpy.zeros((int(aclim)))
    for i in range(len(LARVA)):
        MIN=min((int(aclim),len(LARVA[i])))
        #print(sum(LARVA[i]))
        for j in range(0,MIN):
            moy[j]+=LARVA[i][j]
            nb[j]+=1
        for j in range(0,MIN) :
            for k in range(0,MIN):
                if j<=k :
                    Cov[j,k]+=((LARVA[i][j]-moy[j]/nb[j])*(LARVA[i][k]-moy[k]/nb[k]))
                    All[abs(k-j)]+=((LARVA[i][j]-moy[j]/nb[j])*(LARVA[i][k]-moy[k]/nb[k]))
                    allbi[abs(k-j)]+=1
                #print(((LARVA[i][j]-moy[j]/nb[j])*(LARVA[i][k]-moy[k]/nb[k])))
                    bibu1[j,k]+=1

    Cov=Cov/bibu1

    Covbase=numpy.zeros((int(baselim),int(baselim)))
    bibu=numpy.zeros((int(baselim),int(baselim)))

    moy=numpy.zeros((int(baselim)))
    nb=numpy.zeros((int(baselim)))
    for i in range(len(LARVAbase)):
        #print(sum(LARVA[i]))
        MIN=min((int(baselim),len(LARVAbase[i])))
        for j in range(0,MIN):
            moy[j]+=LARVAbase[i][j]
            nb[j]+=1
        for j in range(0,MIN) :
            for k in range(0,MIN):
                Covbase[j,k]+=((LARVAbase[i][j]-moy[j]/nb[j])*(LARVAbase[i][k]-moy[k]/nb[k]))
                bibu[j,k]+=1

    Covbase=Covbase/bibu


    return Cov, Covbase




def CovarianceTemporelleRoll(comp):  

    moy=numpy.zeros((6,9))
    nb=numpy.zeros((6,9))


    for i in range(len(comp)):
        length = len(sorted(comp[i],key=len, reverse=True)[0])
        y=numpy.array([xi+[None]*(length-len(xi)) for xi in comp[i]])
        for behav in range(0,6):
            if float(behav) in y[:,0]:
                for W in numpy.where(y==float(behav))[0]:
                    for k in range(-4,5):
                        if W+k< len(y) and W+k>0 :
                            moy[behav,k+4]+=y[W+k,1]
                            nb[behav,k+4]+=1  

    
    Cov=numpy.zeros((6,9))
    nbcov=numpy.zeros((6,9))


    for i in range(len(comp)):
        length = len(sorted(comp[i],key=len, reverse=True)[0])
        y=numpy.array([xi+[None]*(length-len(xi)) for xi in comp[i]])
        for behav in range(0,6):
            if float(behav) in y[:,0]:
                for W in numpy.where(y==float(behav))[0]:
                    for k in range(-4,5):
                        if W+k< len(y) and W+k>0 :
                            Cov[behav,k+4]+=((y[W+k,1]-moy[behav,k+4]/nb[behav,k+4])*(y[W,1]-moy[behav,4]/nb[behav,4]))
                            nbcov[behav,k+4]+=1

    return Cov/nbcov




