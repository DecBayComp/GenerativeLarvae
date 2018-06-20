import numpy 
import glob
import re 
import argparse
import pickle
from variable import active_time
from variable import activation_time
from variable import end_time
from variable import processus_nb


def moyenneall(files):

    numbers = re.compile(r'(\d+)')

    def numericalSort(value):
        parts = numbers.split(value)
        parts[1::2] = map(int, parts[1::2])
        return parts

    PROBAS=numpy.zeros(((6,end_time)))
    HISTO=[[],[],[],[],[],[]]
    Color=['#17202a','#c0392b','#8bc34a','#2e86c1','#26c6da','#f1c40f']
    #print(files)

    nb=0
    for name  in sorted(files,key=numericalSort): 
        with open(name, 'rb') as fifi :
            (probaS,suM ,transgene ,matrix3,matrix4,matrixBase3 ,matrixBase4 ,Histo, Histobase,Cov,comp,compbase,Covbase ,Cov2, Covbase2,Covroll)=pickle.load(fifi)
        #print(transgene)
        if nb==0 :
            TRANSGENE=numpy.asarray(transgene)
            MATRIX3=numpy.asarray(matrix3)
            MATRIX3BASE=numpy.asarray(matrixBase3)
            MATRIX4=numpy.asarray(matrix4)
            MATRIX4BASE=numpy.asarray(matrixBase4)
            #print(len(Histo[0]))

            HISTO=Histo
            HISTOBASE=Histobase
            #print(len(Histo),len(Histo[0]))
            COV=numpy.asarray(Cov)
            
            COVBASE=numpy.asarray(Covbase)
            # COV2=Cov2
            # COVBASE2=Covbase2
            COVR=Covroll

        else:

            for i in range(0,6):
                PROBAS[i,:]+=probaS[i,:]/suM
            TRANSGENE=TRANSGENE+numpy.asarray(transgene)
            MATRIX3=MATRIX3+numpy.asarray(matrix3)
            MATRIX3BASE=MATRIX3BASE+numpy.asarray(matrixBase3)
            MATRIX4=MATRIX4+numpy.asarray(matrix4)
            MATRIX4BASE=MATRIX4BASE+numpy.asarray(matrixBase4)
            #print(len(HISTO))
            [HISTO[i].extend(Histo[i]) for i in range(0,6)]
            [HISTOBASE[i].extend(Histobase[i]) for i in range(0,6)]
            
            COV=COV+numpy.asarray(Cov)
            # print(len(Cov))
            COVBASE=COVBASE+numpy.asarray(Covbase)
            # COV2=COV2+Cov2
            # COVBASE2=COVBASE2+Covbase2
            COVR=COVR+Covroll

        nb+=1

    PROBAS=[elemt/nb for elemt in PROBAS]
    for j in range(0,len(TRANSGENE)): 
        for i in range(0,len(TRANSGENE[j])):
            TRANSGENE[j][i]=[elemt/nb for elemt in TRANSGENE[j][i]]
    MATRIX3=[elemt/nb for elemt in list(MATRIX3)]
    MATRIX3BASE=[elemt/nb for elemt in list(MATRIX3BASE)]
    MATRIX4=[elemt/nb for elemt in list(MATRIX4)]
    COVR=[elemt/nb for elemt in COVR]
    MATRIX4BASE=[elemt/nb for elemt in list(MATRIX4BASE)]
    COVBASE=[elemt/nb for elemt in list(COVBASE)]
    COV=[elemt/nb for elemt in list(COV)]
    #print(HISTO)

    #print(Histo[0])

    DATAGENE=(HISTOBASE, HISTO, PROBAS,TRANSGENE,MATRIX3,MATRIX3BASE, MATRIX4, MATRIX4BASE, COVR,comp,compbase, COV,COVBASE)

    return DATAGENE


