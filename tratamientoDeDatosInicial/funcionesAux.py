from Bio import SeqIO
import numpy as np
import os
import matplotlib.pyplot as plt

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


def calcularCorelacion(columnasAltafrec, matriz, cutOff):
    rxy = np.zeros((len(columnasAltafrec), len(columnasAltafrec)))
    
    
    for columna in columnasAltafrec:
        for columna2 in columnasAltafrec:
            if(columna2== columna):
                pass
            else:
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

def plotColumna(ax,columna, matrizFrecuencias):
    frecuencias = [0]*20
    t=0
    for f in matrizFrecuencias:
        frecuencias[t] = f[columna]
        t=t+1
    
    ax.plot(range(1,21),frecuencias, label=f"Posicion: {columna}")
    
           
    
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
    
    columnas =[613,822,830,870,924,1050] 
    fig,ax = plt.subplots()

    for columna in columnas:
        plotColumna(ax,columna, FrecuenciaPorMes)
    
    ax.legend()
    ax.grid()
            
        
                    
            
            
    

    
    
    