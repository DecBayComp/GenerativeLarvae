if __name__ == "__main__":
    # setting the hyper parameters
    import numpy 
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors  
    import parameters 
    import gener 
    import result  
    import function
    import sequence 
    import variable as vble
    import glob
    import re 
    import argparse
    import data_function
    from variable import active_time
    from variable import activation_time
    from variable import end_time
    #from variable import files

    parser = argparse.ArgumentParser()
    parser.add_argument("version")
    parser.add_argument("-m","--matrix2", action="store_true")
    parser.add_argument("-mm","--matrix34", action="store_true")
    parser.add_argument("-s","--sequence", action="store_true")
    parser.add_argument("-d","--distribution", action="store_true")
    parser.add_argument("-c","--Corre", action="store_true")
    parser.add_argument("-v","--Cov", action="store_true")
    parser.add_argument("-vr","--CovRoll", action="store_true")
    args = parser.parse_args()
    print(args.version)

    numbers = re.compile(r'(\d+)')

    def numericalSort(value):
        parts = numbers.split(value)
        parts[1::2] = map(int, parts[1::2])
        return parts
    Color=['#17202a','#c0392b','#8bc34a','#2e86c1','#26c6da','#f1c40f']

    
    files=glob.glob('../donnee/t5/FCF_attP2_1500062@UAS_TNT_2_0003.txt')
    print(files)
    for name  in sorted(files,key=numericalSort): 

        TAB=numpy.loadtxt(name)
        TAB=TAB[~(TAB[:,0]==0.0)]
        TAB=TAB[~(TAB[:,2]==0.0)]
        name=name.split('/')[-1]
        larva=100000
        probabehavior=parameters.proba_ini(TAB)
        transbase, OMEGA_base= parameters.parameters_base(TAB)
        OMEGA, trans = parameters.parameters_active(TAB)
        OMEGA[0,:]=OMEGA_base[:]
        for i in range(0,6):
            plt.plot(numpy.linspace(0,active_time+1,active_time+1),OMEGA[:,i],Color[i])
        plt.savefig('../proba/OMEGA'+args.version+name+'.png') 
        plt.clf()


        print(files)
        gene, comp , compbase, temps_comp = gener.generatif(trans, transbase , OMEGA,probabehavior ,larva, args.version)


        probaS, suM=result.proba(gene)
        PROBA= data_function.probadata(TAB)
        
        fig = plt.figure(figsize=(20,20))
        ax=fig.subplots(1,1)
        for i in range(0,6):
            plt.plot(numpy.linspace(0,end_time-1,end_time),probaS[i,:]/suM,color=Color[i],linewidth=5.0)
            plt.plot(numpy.linspace(0,end_time,end_time*2+1),PROBA[:,i],color=Color[i],linewidth=1.5)
        plt.show()
        plt.ylim(0,1)
        plt.xlim(0, 45)
        plt.savefig('../proba/Proba'+args.version+name+'.png')  
        plt.clf()

        compdata, compbasedata, LARVAdata, LARVAbasedata=data_function.data(TAB)

        if args.matrix2 :

            transdata = result.matrix2(temps_comp)
            pcm3=[]
            [pcm3.append([]) for i in range(17)]
            pcm4=[]
            [pcm4.append([]) for i in range(17)]
            pcminter=[]
            [pcminter.append([]) for i in range(17)]

            fig = plt.figure(figsize=(110,12))
            ax=fig.subplots(3,16)
            pcm3base= ax[0,0].imshow(transbase, cmap='jet', aspect = 'auto')
            pcm4base=ax[1,0].imshow(transdata[0], cmap='jet', aspect = 'auto')
            pcmint=ax[2,0].imshow(transbase-transdata[0], cmap='jet', aspect = 'auto')
            plt.colorbar(pcm3base , ax=ax[0,0])
            plt.colorbar(pcm4base, ax=ax[1,0])
            plt.colorbar(pcmint, ax=ax[2,0])
            ax[0,0].set_title('Base, t='+'[0'+str(activation_time)+']')
            for i in range(1,16):
                ax[0,i].set_title('t='+str(activation_time+i-1))
                pcm3[i]=ax[0,i].imshow(trans[i], cmap='jet', aspect = 'auto')
                pcm4[i]=ax[1,i].imshow(transdata[i+1], cmap='jet', aspect = 'auto')
                pcminter[i]=ax[2,i].imshow(numpy.array(trans[i])-numpy.array(transdata[1+1]), cmap='jet', aspect = 'auto')
                plt.colorbar(pcm3[i] , ax=ax[0,i])
                plt.colorbar(pcm4[i], ax=ax[1,i])
                plt.colorbar(pcminter[i], ax=ax[2,i])
            plt.savefig('../proba/trans'+args.version+name+'.png') 
            plt.clf()

        if args.matrix34 :

            matrix3 , matrix4 = result.matrix34(comp)
            matrixBase3, matrixBase4=result.matrix34(compbase)
            matrix3data , matrix4data = result.matrix34(compdata)
            matrixBase3data, matrixBase4data=result.matrix34(compbasedata)

            fig = plt.figure(figsize=(20,20))
            ax=fig.subplots(2,2)
            pcm=ax[0,0].imshow(matrix3,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm1=ax[0,1].imshow(matrixBase3,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm2=ax[1,0].imshow(matrix3data,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm3=ax[1,1].imshow(matrixBase3data,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm4=ax[2,0].imshow(matrix3data-matrix3,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm5=ax[2,1].imshow(matrixBase3data-matrixBase3,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            for i in ax :
                i.set_xticks(numpy.arange(-0.5, 35.5, 6));
                i.set_yticks(numpy.arange(-0.5, 5.5, 1));
                i.set_xticklabels(numpy.arange(0, 36, 6));
                i.set_yticklabels(numpy.arange(0, 6, 1));

                i.grid(color='w', linestyle='-', linewidth=2)
            #ax[0].grid(color='w', linestyle='-', linewidth=2)

            for pp in [[pcm,ax[0,0]], [pcm1,ax[0,1]], [pcm2,ax[1,0]], [pcm3,ax[1,1]], [pcm4,ax[2,0]], [pcm5,ax[2,1]]]:
                plt.colorbar(pp[0] , ax=pp[1])
            #plt.colorbar(pcm4, ax=ax[1])
            plt.savefig('../matrix/matrice_3_'+args.version+name+'.png')  
            plt.clf()

            fig = plt.figure(figsize=(20,20))
            ax=fig.subplots(2,2)
            pcm=ax[0,0].imshow(matrix4,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm1=ax[0,1].imshow(matrixBase4,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm2=ax[1,0].imshow(matrix4data,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm3=ax[1,1].imshow(matrixBase4data,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm4=ax[2,0].imshow(matrix4data-matrix4,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            pcm5=ax[2,1].imshow(matrixBase4data-matrixBase4,norm=colors.PowerNorm(gamma=1./2.), cmap='jet', aspect = 'auto')
            for i in ax :
                i.set_xticks(numpy.arange(-0.5, 35.5, 6));
                i.set_yticks(numpy.arange(-0.5, 5.5, 1));
                i.set_xticklabels(numpy.arange(0, 36, 6));
                i.set_yticklabels(numpy.arange(0, 6, 1));

                i.grid(color='w', linestyle='-', linewidth=2)
            #ax[0].grid(color='w', linestyle='-', linewidth=2)

            for pp in [[pcm,ax[0,0]], [pcm1,ax[0,1]], [pcm2,ax[1,0]], [pcm3,ax[1,1]], [pcm4,ax[2,0]], [pcm5,ax[2,1]]]:
                plt.colorbar(pp[0] , ax=pp[1])

            #ax[1].grid(color='w', linestyle='-', linewidth=2)
            #ax[0].grid(color='w', linestyle='-', linewidth=2)
            #plt.colorbar(pcm3 , ax=ax[0])
            #plt.colorbar(pcm4, ax=ax[1])
            plt.savefig('../matrix/matrice_4_'+args.version+name+'.png')  
            plt.clf()
            print('matrix fait')


        if args.distribution :

            HISTOdata=[[],[]]
            for k in range(0,2):
                [HISTOdata[k].append([]) for i in range(0,6)]
                if k==0 :
                    TABB=TAB[(TAB[:,0]>activation_time)]
                if k==1 :
                    TABB=TAB[(TAB[:,0]>activation_time)]
                for i in range(0,6) :
                    GMR=TABB[~(TABB[:,1]!=i)]
                    if len(GMR)>0:
                        old=GMR[0,2]
                        for j in range(0,len(GMR)):
                            if(GMR[j,2]!=old):
                                HISTOdata[k][i].append(old)
                                old=GMR[j,2]
                        HISTOdata[k][i].append(old)        


            for i in range(0,2):
                fig = plt.figure(figsize=(20,20))
                ax=fig.subplots(2,3)
                
                Histo=result.distrib_time(comp,compbase)
                HISTOGMR=HISTOdata[i]
                for j in range(0,6):
                    
                    weights = numpy.ones_like(Histo[i][0][int(j)])/float(len(Histo[i][0][int(j)]))*20
                    m, bins, patches=ax[j%2,int(j/2)].hist(Histo[i][0][int(j)],range=[0,5],bins=100,histtype='step',weights=weights,color=Color[j],linewidth=2.0)
                    weightsGMR = numpy.ones_like(HISTOGMR[int(j)])/float(len(HISTOGMR[int(j)]))*20
                    m2, bins2, patches2=ax[j%2,int(j/2)].hist(HISTOGMR[int(j)],range=[0,5],bins=100,histtype='step',weights=weightsGMR, color='black',linewidth=2.0)
                    ax[j%2,int(j/2)].clear()
                    ax[j%2,int(j/2)].set_yscale('log')
                    ax[j%2,int(j/2)].plot(bins2[:-1]+0.5, m2, 'r--', linewidth=2,color=Color[j])
                    ax[j%2,int(j/2)].plot(bins[:-1]+0.5, m, linewidth=4,color=Color[j])
                    ax[j%2,int(j/2)].set_xlabel('Temps')
                    ax[j%2,int(j/2)].legend()  
                    ax[j%2,int(j/2)].patch.set_facecolor(Color[j])
                    ax[j%2,int(j/2)].patch.set_alpha(0.05)    

                plt.savefig('../distrib/distribution_timet5'+args.version+name+str(i)+'.png')  
                plt.clf()    
            print('distribution fait')


        if args.sequence :
            fig = plt.figure(figsize=(20,20))
            sequence.sequence_timecolor2(compbase,comp,Color)
            plt.savefig('../sequence/sequencet5'+args.version+name+'.png')
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

            Cov, Covbase =result.CovarianceTemporelle(LARVAdata,LARVAbasedata)
            
            fig = plt.figure(figsize=(40,20))
            ax=fig.subplots(1,2)
            pcm3=ax[0].imshow(Cov, cmap='jet', aspect = 'auto')
            pcm4=ax[1].imshow(Covbase, cmap='jet', aspect = 'auto')
            fig.colorbar(pcm3, ax=ax[0])
            fig.colorbar(pcm4,ax=ax[1])
            plt.savefig('../Covariance/Covariancedata'+args.version+name+'.png')
            plt.clf()
            
            timebase=[[]]
            time=[[]]
            
            for i in range(0,len(comp)):
                time.append([])
                timebase.append([])
                [time[-1].append(comp[i][j][1]) for j in range(0,len(comp[i]))]
                [timebase[-1].append(compbase[i][j][1]) for j in range(0,len(compbase[i]))]

            Cov, Covbase=result.CovarianceTemporelle(time,timebase)         
            
            fig = plt.figure(figsize=(45,20))
            ax=fig.subplots(1,2)
            pcm3=ax[0].imshow(Cov, cmap='jet', aspect = 'auto')
            pcm4=ax[1].imshow(Covbase, cmap='jet', aspect = 'auto')
            fig.colorbar(pcm3, ax=ax[0])
            fig.colorbar(pcm4,ax=ax[1])
            plt.savefig('../Covariance/Covariance'+args.version+name+'.png')
            plt.clf()

            Cov2, Covbase2=result.CovarianceTemporelle2(comp,compbase)    
            fig = plt.figure(figsize=(45,20))
            ax=fig.subplots(1,2)
            pcm3=ax[0].imshow(Cov2, cmap='jet', aspect = 'auto')
            pcm4=ax[1].imshow(Covbase2, cmap='jet', aspect = 'auto')
            fig.colorbar(pcm3, ax=ax[0])
            fig.colorbar(pcm4,ax=ax[1])
            plt.savefig('../Covariance/Covariance_2'+args.version+name+'.png')
            plt.clf()


            #compdata, compbasedata = data_function.EXPcomtime(TAB)
            Covdata, Covbasedata=result.CovarianceTemporelle2(compdata,compbasedata)    
            fig = plt.figure(figsize=(45,20))
            ax=fig.subplots(1,2)
            pcm3=ax[0].imshow(Covdata, cmap='jet', aspect = 'auto')
            pcm4=ax[1].imshow(Covbasedata, cmap='jet', aspect = 'auto')
            fig.colorbar(pcm3, ax=ax[0])
            fig.colorbar(pcm4,ax=ax[1])
            plt.savefig('../Covariance/Covariance_2_data_'+args.version+name+'.png')
            plt.clf()

            fig = plt.figure(figsize=(20,40))
            ax=fig.subplots(2,1)
            for i in range(0,6): 
                ax[0].plot(Cov2[i,:],color=Color[i])
                ax[1].plot(Covbase2[i,:],color=Color[i])
                ax[0].plot(Covdata[i,:],color=Color[i])
                ax[1].plot(Covbasedata[i,:],color=Color[i])
            plt.savefig('../Covariance/Covariance_line'+args.version+name+'.png')
            plt.clf()



        if args.CovRoll :
            Covdata=result.CovarianceTemporelleRoll(compdata)
            Cov=result.CovarianceTemporelleRoll(comp)

            fig = plt.figure(figsize=(20,40))
            ax=fig.subplots(6,1)
            for i in range(0,6):
                ax[i].plot(Cov[i,:],color=Color[i])
                ax[i].plot(Covdata[i,:],color='#04B486')

            plt.savefig('../Covariance/Covariance_Roll'+args.version+name+'.png')
            plt.clf()














