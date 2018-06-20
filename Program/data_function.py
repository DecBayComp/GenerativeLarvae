import numpy 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors  
from variable import activation_time
from variable import active_time
from variable import end_time



def data(TAB):
    compbase=[]
    LARVAbase=[]
    B=numpy.zeros((1,2))
    TABB=TAB[(TAB[:,0]<activation_time)]
    for i in range(int(TABB[0,3]),int(TABB[-1,3])+1):
        T=TABB[(TABB[:,3])==i]
        if len(T)>0 :
            TABLEAU=numpy.zeros((len(T),2))    
            a=0
            for j in range(1,len(T)):
                if T[j-1,2]!=T[j,2]:
                    TABLEAU[a,0]=int(T[j-1,1])
                    TABLEAU[a,1]=T[j-1,2]
                    a+=1
            TABLEAU[a,0]=int(T[-1,1])
            TABLEAU[a,1]=T[-1,2]
            TABLEAU=TABLEAU[(TABLEAU[:,1]!=0.0)]
            #B=numpy.concatenate((B,TABLEAU))
            compbase.extend([TABLEAU.tolist()])
            LARVAbase.extend([list(TABLEAU[:,1])])
    #compbase=B[1:,:].tolist()
    #LARVAbase=list(B[1:,1])

    comp=[]
    LARVA=[]
    B=numpy.zeros((1,2))
    TABB=TAB[(TAB[:,0]>activation_time)]
    for i in range(int(TABB[0,3]),int(TABB[-1,3])+1):
        T=TABB[(TABB[:,3])==i]
        if len(T)>0 :
            TABLEAU=numpy.zeros((len(T),2))    
            a=0
            for j in range(1,len(T)):
                if T[j-1,2]!=T[j,2]:
                    TABLEAU[a,0]=int(T[j-1,1])
                    TABLEAU[a,1]=T[j-1,2]
                    a+=1
            TABLEAU[a,0]=int(T[-1,1])
            TABLEAU[a,1]=T[-1,2]
            TABLEAU=TABLEAU[(TABLEAU[:,1]!=0.0)]
            comp.extend([TABLEAU.tolist()])
            LARVA.extend([list(TABLEAU[:,1])])
   # comp=B[1:,:].tolist()
   # LARVA=list(B[1:,1])

    return comp, compbase, LARVA, LARVAbase

def probadata(TAB) :
    PROBA=numpy.zeros(((end_time)*2+1,6))
    NB=numpy.zeros(((end_time)*2+1))

    for i in range(0,len(TAB)):
        if int(TAB[i,0]*2)<91 :
            PROBA[int(TAB[i,0]*2),int(TAB[i,1])]+=1
            NB[int(TAB[i,0]*2)]+=1

    for i in range(0,6) :
        PROBA[:,i]=PROBA[:,i]/NB

    numpy.nan_to_num(PROBA)

    return PROBA


# def EXPcomtime(TAB):

#     LARVAbase=[[]]
#     LARVA=[[]]

#     nblarva=0
#     TABB=TAB[(TAB[:,0]<activation_time)]
#     for i in range(int(TABB[-1,4])) :
#         T=TABB[(TABB[:,4]==float(i))]
#         if len(T)>0 :
#             for j in range(int(T[-1,3])):
#                 T2=T[(T[:,3]==float(j))]
#                 if len(T2)>0:
#                     nblarva+=1
#                     LARVAbase.extend([[[int(T2[0,1]),T2[0,2]]]])
#                     old_time=T2[0,2]
#                     for k in range(0,len(T2)):
#                         if old_time!=T2[k,2]:
#                             LARVAbase[nblarva].extend([[int(T2[k,1]),T2[k,2]]])
#                             old_time=T2[k,2]


#     nblarva=0
#     TABB=TAB[(TAB[:,0]>activation_time)]
#     old_time=TABB[0,2]
#     for i in range(int(TABB[-1,4])) :
#         T=TABB[(TABB[:,4]==float(i))]   
#         if len(T)>0 :
#             for j in range(int(T[-1,3])):
#                 T2=T[(T[:,3]==float(j))]
#                 if len(T2)>0:
#                     nblarva+=1
#                     LARVA.extend([[[int(T2[0,1]),T2[0,2]]]])
#                     old_time=T2[0,2]
#                     for k in range(0,len(T2)):
#                         if old_time!=T2[k,2]:
#                             LARVA[nblarva].extend([[int(T2[k,1]),T2[k,2]]])
#                             old_time=T2[k,2]

#     return LARVA, LARVAbase


# def larvaEXP(TAB):
#     TABB=TAB[(TAB[:,0]<activation_time)]
#     old_time=TABB[0,2]
#     old_larva=TABB[0,3]
#     LARVAbase=[[]]
#     LARVA=[[]]
#     for i in range(int(TABB[-1,4])) :
#         T=TABB[(TABB[:,4]==float(i))]   
#         for j in range(int(TABB[-1,3])):
#             T2=T[(T[:,3]==float(j))]
#             if len(T2)>0:
#                 LARVAbase.extend(([list(set(list(T2[:,2])))]))

#     TABB=TAB[(TAB[:,0]>activation_time)]
#     TABB=TABB[(TABB[:,0]<end_time)]
#     old_time=TABB[0,2]
#     old_larva=TABB[0,3]
#     for i in range(int(TABB[-1,4])) :
#         T=TABB[(TABB[:,4]==float(i))]   
#         for j in range(int(TABB[-1,3])):
#             T2=T[(T[:,3]==float(j))]
#             if len(T2)>0:
#                 LARVA.extend(([list(set(list(T2[:,2])))]))
#           #  if sum(LARVA[-1])>active_time :
#               #  print(sum(LARVA[-1]),len(LARVA)-1,j,i)
    

#     print('bloblobase')
#     return LARVAbase, LARVA



# def FonctionCorrelationEXP(LARVA,LARVAbase) :
#     Fonction=np.zeros((40))
#     Fonctionbase=np.zeros((40))

#     #print('bloblo',len(LARVA))        
#     for i in range(0,len(LARVAbase)):
#         delta=1
#         while sum(LARVAbase[i][0:delta])<10 : 
#             if delta>len(LARVAbase[i]) :
#                 break
#             for j in range(0,len(LARVAbase[i])):
#                 if int((sum(LARVAbase[i][j:j+delta])))*4<40 :
#                     Fonctionbase[int((sum(LARVAbase[i][j:j+delta]))*4)]+=1/(len(LARVAbase)+len(LARVAbase[i]))
#             delta+=1
    
#     for i in range(0,len( LARVA)):
#         delta=1
#         while sum(LARVA[i][0:delta])<10 :
#             if delta>len(LARVA[i]) :
#                 break
#             for j in range(0,len(LARVA[i])):
#                 if int((sum(LARVA[i][j:j+delta])))*4<40 :
#                     #print(int((sum(LARVA[i][j:j+delta]))*2),delta)
#                     Fonction[int((sum(LARVA[i][j:j+delta]))*4)]+=1/(len(LARVA)+len(LARVA[i]))
#             delta+=1
#     plt.plot(np.linspace(0,10,40),Fonction, color = '#2998C9')
#     plt.plot(np.linspace(0,10,40),Fonctionbase,color='#E97106')
#     return


# def CovarianceTemporelleEXT(LARVA,LARVAbase):  
#     fig = plt.figure(figsize=(40,20))
#     ax=fig.subplots((2,1))

#     Cov=np.zeros((activ_time*4,activ_time*4))

#     for i in range(LARVA):
#         moy=0.0
#         for j in range(len(LARVA[i])):
#             moy+=LARVA[i][j]/len(LARVA[i])
#         for j in range(1,len(LARVA[i])) :
#             Cov(int(LARVA[i][j-1]*4),int(LARVA[i][j]*4))=((LARVA[i][j-1]-moy)*(LARVA[i][j]*4-moy))/(len(LARVA[i])-1)/len(LARVA)

#     ax[1].imshow(Cov,map=viridis)

#     Cov=np.zeros((activ_time*4,activ_time*4))

#     for i in range(LARVA):
#         moy=0.0
#         for j in range(len(LARVA[i])):
#             moy+=LARVA[i][j]/len(LARVA[i])
#         for j in range(1,len(LARVA[i])) :
#             Cov(int(LARVA[i][j-1]*4),int(LARVA[i][j]*4))=((LARVA[i][j-1]-moy)*(LARVA[i][j]*4-moy))/(len(LARVA[i])-1)/len(LARVA)

#     ax[0].imshow(Cov,map=viridis)

#     return fig







