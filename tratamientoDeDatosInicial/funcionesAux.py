from Bio import SeqIO
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

def armarMatrizDeFasta(archivo):
    matriz = []
    i=0
    for record in SeqIO.parse(archivo, "fasta"):
        line =str(record.seq)
        matriz.append(line)

    
    return matriz

def calcularFrecuenciaColumna(matriz, columna, secuenciaRef):

    aminoacidoRef = secuenciaRef[columna]

    contador = 0
    i = 0
    for i in range(len(matriz)):
        elemento = matriz[i]
 
        if ((elemento[columna] != aminoacidoRef) and (elemento[columna] !='X')):
            contador = contador + 1

    return contador/len(matriz)


def calcularCorelacion(frecuenciasMes):
    rxy = np.zeros((len(frecuenciasMes[0]), len(frecuenciasMes[0])))
    
    
    for fila in range(len(rxy)):
        for columna in range(len(rxy)):
            if(fila == columna):
                rxy[fila,columna] =0
                pass
            else:
                correlacion,pvalue = pearsonr(x=frecuenciasMes[:,fila],y=frecuenciasMes[:,columna])
                rxy[fila,columna] = correlacion
                pass 
    
    return rxy

def calcularFrecuencias(archivo,matriz, secuenciaRef):
        
    numColumnas = len(secuenciaRef)
    columnasAltafrec = []
    frecuencias = np.zeros(numColumnas)
    for columna in range(numColumnas):
        frecuencias[columna] = calcularFrecuenciaColumna(matriz, columna, secuenciaRef)
        # if frecuecia>fCutOff:
        #     columnasAltafrec.append(columna)
    
   
    return frecuencias

def plotColumnas(columnas,FrecuenciaPorMes):
    columnas =[613,822,830,870,924,1050] 
    fig,ax = plt.subplots()

    for columna in columnas:
        plotColumna(ax,columna, FrecuenciaPorMes)
    
    return(fig,ax)

def plotColumna(ax,columna, matrizFrecuencias):
    frecuencias = [0]*20
    t=0
    for f in matrizFrecuencias:
        frecuencias[t] = f[columna]
        t=t+1
    meses =["5/20","6/20","7/20","8/20","9/20","10/20","11/20","12/20",
            "1/21","2/21","3/21","4/21","5/21","6/21","7/21","8/21","9/21","10/21","11/21","12/21",]
    ax.plot(meses,frecuencias, label=f"Posicion: {columna}")
    
           
    
#hay que analisar los caracteres x e - .

if __name__ == "__main__":
    
    carpetaArchivosAlineados = "archivosFasta/argentina/Alineados/"
    
    FrecuenciaPorMes = [[]]*20#20 meses vamos a tomar
    
    for file in os.listdir(carpetaArchivosAlineados): 
        if file.endswith('.fasta'):
            print(f"procesando archivo {file} ...")
            
            archivo = open( carpetaArchivosAlineados + file, "r")
            
            matriz = armarMatrizDeFasta(archivo)
            
            #secuencia de wohan
            secuenciaRef = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*"
           
            frecuenciasMes=calcularFrecuencias(archivo, matriz, secuenciaRef)
            
            mes = int(file.split('-')[1])
            if int(file.split('-')[2].split('_')[0]) == 2021 :
                FrecuenciaPorMes[mes+7] =frecuenciasMes
            else:
                FrecuenciaPorMes[mes-5] =frecuenciasMes
        
    FrecuenciaPorMes = np.array(FrecuenciaPorMes)
        
    matrizAdjacencia = calcularCorelacion(FrecuenciaPorMes)
        
    # columnas =[613,822,830,870,924,1050] 
    # fig,ax = plotColumnas(columnas,FrecuenciaPorMes)
  
    # ax.legend()
    # ax.set_ylabel("frecuencia de mutacion")
    # ax.set_xlabel("Mes")
    # ax.grid()
    theshold =0.7

    matrizAdjacencia[np.isnan(matrizAdjacencia)] = 0
    matrizAdjacencia[matrizAdjacencia<theshold]=0
    matrizAdjacencia[matrizAdjacencia>=theshold]=1
    
            
        
                    
            
            
    

    
    
    