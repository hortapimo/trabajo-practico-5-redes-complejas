from obtenerGraphNetworkx import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import funcionesAux as f


carpetaArchivosAlineados = "archivosFasta/Argentina/Alineados/"
sol = obtenerGrafo(carpetaArchivosAlineados, theshold=0.9, thresholdFrec = 0.05, nameGraph='Argentina_09_05')

##heat map
ax =sns.heatmap(sol[2][:,350:510])
ax.set_xticklabels(np.linspace(350, 510,27,dtype=int))
ax.set_yticklabels(np.linspace(1, 44,22,dtype=int))
ax.set_xlabel("# Residuo")
ax.set_ylabel("# Mes empezando en marzo 2020")

size = 18
ax.xaxis.label.set_size(size) 
ax.yaxis.label.set_size(size) 
ax.tick_params(axis='both', which='major', labelsize=size-4)

cbar = ax.collections[0].colorbar 
cbar.ax.tick_params(labelsize=size)


## columnas

fig,ax=f.plotColumnas([18,477,678,680,141,654,795,338,370,497,483,763,953,504,374,407,404,476,968,953,500,372,26,375], sol[2], pais="Arg")#le saque el mes de enero porque tenia muyy pocas frecuencias
size = 24
fig.set_size_inches(10, 7)
ax.xaxis.label.set_size(size) 
ax.yaxis.label.set_size(size) 
ax.tick_params(axis='both', which='major', labelsize=size-4)
ax.legend(prop={'size': size-2})

fig2,ax2=f.plotColumnas([477,338,370,497,483,504,374,407,404,476,500,372,375], sol[2], pais="Arg")#le saque el mes de enero porque tenia muyy pocas frecuencias
size = 24
fig2.set_size_inches(10, 7)
ax2.xaxis.label.set_size(size) 
ax2.yaxis.label.set_size(size) 
ax2.tick_params(axis='both', which='major', labelsize=size-4)

ax2.legend(prop={'size': size-2})

#f.plotColumnas([], sol[2], pais="Arg")
