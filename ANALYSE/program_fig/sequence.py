import numpy 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors  
global  end_time, active_time, activation_time

def sequence_timecolor2(compbase,comp, Color): 
	for i in range(0,100) :
		x=0
		for j in range(len(compbase[i])):
			plt.plot([x,x+compbase[i][j][1]],[i,i],Color=Color[compbase[i][j][0]],linewidth=5.0)
			x+=compbase[i][j][1]
		x=30
		for j in range(len(comp[i])):
			plt.plot([x,x+comp[i][j][1]],[i,i],Color=Color[comp[i][j][0]],linewidth=5.0)
			x+=comp[i][j][1]
	return 


def sequence_timecolor(gene, Color): 
	for i in range(0,50) :
		for j in range(len(gene[i,:])):
			plt.plot([j,j+1],[i,i],Color=Color[int(gene[i,j])],linewidth=10.0)
	return			


def sequence_time(compbase,comp, Color): 
	for i in range(0,50) :
		x=0
		for j in range(len(compbase[i])):
			plt.plot([x,x+compbase[i][j][1]],[compbase[i][j][0],compbase[i][j][0]],Color=Color[compbase[i][j][0]],linewidth=5.0)
			x+=compbase[i][j][1]
		for j in range(len(comp[i])):
			plt.plot([x,x+comp[i][j][1]],[comp[i][j][0],comp[i][j][0]],Color=Color[comp[i][j][0]],linewidth=5.0)
			x+=comp[i][j][1]	
	return 


