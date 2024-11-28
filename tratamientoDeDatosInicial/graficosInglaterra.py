from obtenerGraphNetworkx import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import funcionesAux as f


carpetaArchivosAlineados = "archivosFasta/England/Alineados/"
sol = obtenerGrafo(carpetaArchivosAlineados, theshold=0.8, thresholdFrec = 0.05, nameGraph='England_08_0')

## columnas
fig,ax=f.plotColumnas([569,500,680,715,981,1117], sol[2][1:,:], pais="Eng")#le saque el mes de enero porque tenia muyy pocas frecuencias
size = 24
fig.set_size_inches(10, 7)
ax.xaxis.label.set_size(size) 
ax.yaxis.label.set_size(size) 
ax.tick_params(axis='both', which='major', labelsize=size-4)
ax.legend(prop={'size': size-2})

fig2,ax2=f.plotColumnas([17,221,261,271,582], sol[2][1:,:], pais="Eng")#le saque el mes de enero porque tenia muyy pocas frecuencias
size = 24
fig2.set_size_inches(10, 8)
ax2.xaxis.label.set_size(size) 
ax2.yaxis.label.set_size(size) 
ax2.tick_params(axis='both', which='major', labelsize=size-4)
ax2.legend(prop={'size': size-2})

fig3,ax3=f.plotColumnas([95,141,416,483,700,956], sol[2][1:,:], pais="Eng")#le saque el mes de enero porque tenia muyy pocas frecuencias
size = 24
fig3.set_size_inches(10, 8)
ax3.xaxis.label.set_size(size) 
ax3.yaxis.label.set_size(size) 
ax3.tick_params(axis='both', which='major', labelsize=size-4)
ax3.legend(prop={'size': size-2})