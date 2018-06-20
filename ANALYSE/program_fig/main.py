if __name__ == "__main__":
    # setting the hyper parameters
    import numpy 
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors  
    import result  
    import sequence 
    import variable as vble
    import glob
    import re 
    import argparse
    from variable import active_time
    from variable import activation_time
    from variable import end_time
    import resultsfinal
    import pickle
    #from variable import files

    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--matrix2", action="store_true")
    parser.add_argument("-mm","--matrix34", action="store_true")
    parser.add_argument("-s","--sequence", action="store_true")
    parser.add_argument("-d","--distribution", action="store_true")
    parser.add_argument("-c","--Corre", action="store_true")
    parser.add_argument("-v","--Cov", action="store_true")
    parser.add_argument("-vr","--CovRoll", action="store_true")
    args = parser.parse_args()

    numbers = re.compile(r'(\d+)')
    args.version='version3' 

    def numericalSort(value):
        parts = numbers.split(value)
        parts[1::2] = map(int, parts[1::2])
        return parts
    Color=['#17202a','#c0392b','#8bc34a','#2e86c1','#26c6da','#f1c40f']
    DOSSIER='GMR_11A07_AE_01_UAS_TNT_2_0003_p_8_45s1x30s0s#p_8_105s10x2s10s'
    files=glob.glob('../../PROGRAMME_TARS/Liste_fig/'+DOSSIER+'/*')
    #name_0 in sorted(glob.glob('../DATA/TAB'+DOSSIER+'*'),key=numericalSort): 
    name_0='../../PROGRAMME_TARS/DATA/TABt5'+DOSSIER+'.pkl'
    name_files=open(name_0,'rb')
    print(name_files)
    (trans,transbase,compdata, compbasedata, LARVEdata, LARVAbasedata, PROBA , matrix3data , matrix4data ,matrixBase3data, matrixBase4data,histodata, histodatabase, Covdata, Covbasedata)=pickle.load(name_files)

    # for i in range(0,6):
    #     plt.plot(numpy.linspace(0,active_time+1,active_time+1),OMEGA[:,i],Color[i])
    # plt.savefig('../proba/OMEGA'+args.version+name+'.png') 
    # plt.clf()
    (histobase, histo, probaS,transgene ,matrix3,matrixBase3,matrix4,matrixBase4,Covroll, comp,compbase, Cov, Covbase)=resultsfinal.moyenneall(files)
    #print(transgene)
    #print(len(histo),len(histo[0]))
    #print(len(PROBA),PROBA)
    fig = plt.figure(figsize=(20,20))
    ax=fig.subplots(1,1)
    for i in range(0,6):
        
        plt.plot(numpy.linspace(0,end_time-1,end_time),probaS[i],color=Color[i],linewidth=5.0)
        plt.plot(numpy.linspace(0,end_time,end_time*2+1),PROBA[i,:],color=Color[i],linewidth=1.5)
    plt.show()
    plt.ylim(0,1)
    plt.xlim(0, end_time-1)
    plt.savefig('../proba/Proba'+args.version+DOSSIER+'.png')  
    plt.clf()

    if args.matrix2 :

        #transdata = result.matrix2(temps_comp)
        pcm3=[]
        [pcm3.append([]) for i in range(17)]
        pcm4=[]
        [pcm4.append([]) for i in range(17)]
        pcminter=[]
        [pcminter.append([]) for i in range(17)]

        fig = plt.figure(figsize=(110,12))
        ax=fig.subplots(3,16)
        pcm3base= ax[0,0].imshow(transbase, cmap='jet', aspect = 'auto')
        pcm4base=ax[1,0].imshow(transgene[0], cmap='jet', aspect = 'auto')
        pcmint=ax[2,0].imshow(transbase-transgene[0], cmap='jet', aspect = 'auto')
        plt.colorbar(pcm3base , ax=ax[0,0])
        plt.colorbar(pcm4base, ax=ax[1,0])
        plt.colorbar(pcmint, ax=ax[2,0])
        ax[0,0].set_title('Base, t='+'[0'+str(activation_time)+']')
        for i in range(1,16):
            ax[0,i].set_title('t='+str(activation_time+i-1))
            pcm3[i]=ax[0,i].imshow(trans[i], cmap='jet', aspect = 'auto')
            pcm4[i]=ax[1,i].imshow(transgene[i+1], cmap='jet', aspect = 'auto')
            pcminter[i]=ax[2,i].imshow(numpy.array(trans[i])-numpy.array(transgene[1+1]), cmap='jet', aspect = 'auto')
            plt.colorbar(pcm3[i] , ax=ax[0,i])
            plt.colorbar(pcm4[i], ax=ax[1,i])
            plt.colorbar(pcminter[i], ax=ax[2,i])
        plt.savefig('../proba/trans'+args.version+DOSSIER+'.png') 
        plt.clf()

    if args.matrix34 :

        # matrix3 , matrix4 = result.matrix34(comp)
        # matrixBase3, matrixBase4=result.matrix34(compbase)
        # matrix3data , matrix4data = result.matrix34(compdata)
        # matrixBase3data, matrixBase4data=result.matrix34(compbasedata)

        fig = plt.figure(figsize=(20,20))
        ax=fig.subplots(3,2)
        pcm=ax[0,0].imshow(matrix3,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
        pcm1=ax[0,1].imshow(matrixBase3,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
        pcm2=ax[1,0].imshow(matrix3data,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
        pcm3=ax[1,1].imshow(matrixBase3data,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
        pcm4=ax[2,0].imshow(matrix3data-matrix3, cmap='jet', aspect = 'auto')
        pcm5=ax[2,1].imshow(matrixBase3data-matrixBase3, cmap='jet', aspect = 'auto')
        for tot in ax :
            for axes in tot:
                print(axes)
                axes.set_xticks(numpy.arange(-0.5, 35.5, 6));
                axes.set_yticks(numpy.arange(-0.5, 5.5, 1));
                axes.set_xticklabels(numpy.arange(0, 36, 6));
                axes.set_yticklabels(numpy.arange(0, 6, 1));

                axes.grid(color='w', linestyle='-', linewidth=2)
        #ax[0].grid(color='w', linestyle='-', linewidth=2)

        for pp in [[pcm,ax[0,0]], [pcm1,ax[0,1]], [pcm2,ax[1,0]], [pcm3,ax[1,1]], [pcm4,ax[2,0]], [pcm5,ax[2,1]]]:
            plt.colorbar(pp[0] , ax=pp[1])
        #plt.colorbar(pcm4, ax=ax[1])
        plt.savefig('../matrix/matrice_3_'+args.version+DOSSIER+'.png')  
        plt.clf()

        fig = plt.figure(figsize=(20,20))
        ax=fig.subplots(3,2)
        pcm=ax[0,0].imshow(matrix4,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
        pcm1=ax[0,1].imshow(matrixBase4,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
        pcm2=ax[1,0].imshow(matrix4data,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
        pcm3=ax[1,1].imshow(matrixBase4data,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
        pcm4=ax[2,0].imshow(matrix4data-matrix4, cmap='jet', aspect = 'auto')
        pcm5=ax[2,1].imshow(matrixBase4data-matrixBase4, cmap='jet', aspect = 'auto')
        for tot in ax :
            for axes in tot:
                print(axes)
                axes.set_xticks(numpy.arange(-0.5, 35.5, 6));
                axes.set_yticks(numpy.arange(-0.5, 35.5, 6));
                axes.set_xticklabels(numpy.arange(0, 36, 6));
                axes.set_yticklabels(numpy.arange(0, 36, 6));
                axes.grid(color='w', linestyle='-', linewidth=2)
        #ax[0].grid(color='w', linestyle='-', linewidth=2)

        for pp in [[pcm,ax[0,0]], [pcm1,ax[0,1]], [pcm2,ax[1,0]], [pcm3,ax[1,1]], [pcm4,ax[2,0]], [pcm5,ax[2,1]]]:
            plt.colorbar(pp[0] , ax=pp[1])

        #ax[1].grid(color='w', linestyle='-', linewidth=2)
        #ax[0].grid(color='w', linestyle='-', linewidth=2)
        #plt.colorbar(pcm3 , ax=ax[0])
        #plt.colorbar(pcm4, ax=ax[1])
        plt.savefig('../matrix/matrice_4_'+args.version+DOSSIER+'.png')  
        plt.clf()
        print('matrix fait')


    if args.distribution :
        
        for i in range(0,2):
            fig = plt.figure(figsize=(20,20))
            ax=fig.subplots(2,3)
            if i==0 :
                Histo=histo
                HISTOGMR=histobase
                nom="inactive"

            if i==1 :
                Histo=histobase
                HISTOGMR=histodatabase
                nom="base"
            #print(Histo)
            for j in range(0,6):
                if not Histo[j]:
                    print('c est vide' )

                else:
                    print(len(Histo[j]))
                
                    weights = numpy.ones_like(Histo[j])/float(len(Histo[j]))*20
                    m, bins, patches=ax[j%2,int(j/2)].hist(Histo[j],range=[0,5],bins=100,histtype='step',weights=weights,color=Color[j],linewidth=2.0)
                    
                    weightsGMR = numpy.ones_like(HISTOGMR[int(j)])/float(len(HISTOGMR[int(j)]))*20
                    m2, bins2, patches2=ax[j%2,int(j/2)].hist(HISTOGMR[int(j)],range=[0,5],bins=100,histtype='step',weights=weightsGMR, color='black',linewidth=2.0)
                    ax[j%2,int(j/2)].clear()
                    ax[j%2,int(j/2)].set_yscale('log')
                    ax[j%2,int(j/2)].plot(bins2[:-1]+0.5, m2, 'r--', linewidth=2,color=Color[j])
                    ax[j%2,int(j/2)].plot(bins[:-1]+0.5, m, linewidth=4,color=Color[j])
                    #ax[j%2,int(j/2)].set_xlabel('Temps')
                    ax[j%2,int(j/2)].legend()  
                    ax[j%2,int(j/2)].patch.set_facecolor(Color[j])
                    ax[j%2,int(j/2)].patch.set_alpha(0.05)    

            plt.savefig('../distrib/distribution_timet5'+args.version+DOSSIER+nom+'.png')  
            plt.clf()    
        print('distribution fait')


    if args.sequence :
        fig = plt.figure(figsize=(20,20))
        sequence.sequence_timecolor2(compbase,comp,Color)
        plt.savefig('../sequence/sequencet5'+args.version+DOSSIER+'.png')
        plt.clf()
    

    if args.Corre:
    
        fig = plt.figure(figsize=(20,20))
        # old=TAB[0,2]
        # Corre=[]
        # for i in range(0,len(TAB)) :
        #     if TAB[0,2]!=old :
        #         Corre.extend([list(TAB[a:i,2])])
        #         a=i
        #         old=TAB[0,2]
            
        # Fonction=numpy.zeros((end_time)*2)
        # for i in range(0, len(Corre)) :
        #     spatial=[]
        #     [spatial.append(spatial[-1]+Corre[j]) for j in range(1,len(Corre[i]))]
        #     for j in range(0,len(spatial)) :
        #         if int(spatial[j]*2)< (end_time)*2:
        #             Fonction[int(spatial[j]*2)]+=1/len(compbase)
        # plt.plot(numpy.linspace(0,end_time,end_time*2),Fonction/0.5)  
        #data_function.FonctionCorrelationEXP(LARVA,LARVAbase) 
        # print('EXP')
        # result.FonctionCorrelationMOD(comp,compbase)
        # plt.savefig('../Correlation'+args.version+name+'.png')
        # plt.clf()

        print('Correlation fait')

    if args.Cov:

        #Cov, Covbase =result.CovarianceTemporelle(LARVAdata,LARVAbasedata)
        
        fig = plt.figure(figsize=(40,20))
        ax=fig.subplots(1,2)
        pcm3=ax[0].imshow(Covdata, cmap='jet', aspect = 'auto')
        pcm4=ax[1].imshow(Covbasedata, cmap='jet', aspect = 'auto')
        fig.colorbar(pcm3, ax=ax[0])
        fig.colorbar(pcm4,ax=ax[1])
        plt.savefig('../Covariance/Covariancedata'+args.version+DOSSIER+'.png')
        plt.clf()
        
        timebase=[[]]
        time=[[]]
        
        for i in range(0,len(comp)):
            time.append([])
            timebase.append([])
            [time[-1].append(comp[i][j][1]) for j in range(0,len(comp[i]))]
            [timebase[-1].append(compbase[i][j][1]) for j in range(0,len(compbase[i]))]

        #Cov, Covbase=result.CovarianceTemporelle(time,timebase)         
        
        fig = plt.figure(figsize=(45,20))
        ax=fig.subplots(1,2)
        pcm3=ax[0].imshow(Cov, cmap='jet', aspect = 'auto')
        pcm4=ax[1].imshow(Covbase, cmap='jet', aspect = 'auto')
        fig.colorbar(pcm3, ax=ax[0])
        fig.colorbar(pcm4,ax=ax[1])
        plt.savefig('../Covariance/Covariance'+args.version+DOSSIER+'.png')
        plt.clf()

        # Cov2, Covbase2=result.CovarianceTemporelle2(comp,compbase)    
        # fig = plt.figure(figsize=(45,20))
        # ax=fig.subplots(1,2)
        # pcm3=ax[0].imshow(Cov2, cmap='jet', aspect = 'auto')
        # pcm4=ax[1].imshow(Covbase2, cmap='jet', aspect = 'auto')
        # fig.colorbar(pcm3, ax=ax[0])
        # fig.colorbar(pcm4,ax=ax[1])
        # plt.savefig('../Covariance/Covariance_2'+args.version+DOSSIER+'.png')
        # plt.clf()


        #compdata, compbasedata = data_function.EXPcomtime(TAB)
        #Covdata, Covbasedata=result.CovarianceTemporelle2(compdata,compbasedata)    
        # fig = plt.figure(figsize=(45,20))
        # ax=fig.subplots(1,2)
        # pcm3=ax[0].imshow(Covdata, cmap='jet', aspect = 'auto')
        # pcm4=ax[1].imshow(Covbasedata, cmap='jet', aspect = 'auto')
        # fig.colorbar(pcm3, ax=ax[0])
        # fig.colorbar(pcm4,ax=ax[1])
        # plt.savefig('../Covariance/Covariance_2_data_'+args.version+DOSSIER+'.png')
        # plt.clf()

        # fig = plt.figure(figsize=(20,40))
        # ax=fig.subplots(2,1)
        # for i in range(0,6): 
        #     ax[0].plot(Cov2[i,:],color=Color[i])
        #     ax[1].plot(Covbase2[i,:],color=Color[i])
        #     ax[0].plot(Covdata[i,:],color=Color[i])
        #     ax[1].plot(Covbasedata[i,:],color=Color[i])
        # plt.savefig('../Covariance/Covariance_line'+args.version+DOSSIER+'.png')
        # plt.clf()



    if args.CovRoll :
        Covdata=result.CovarianceTemporelleRoll(compdata)
        Cov=result.CovarianceTemporelleRoll(comp)

        fig = plt.figure(figsize=(20,40))
        ax=fig.subplots(6,1)
        for i in range(0,6):
            ax[i].plot(Cov[i,:],color=Color[i])
            ax[i].plot(Covdata[i,:],color='#04B486')

        plt.savefig('../Covariance/Covariance_Roll'+args.version+DOSSIER+'.png')
        plt.clf()














