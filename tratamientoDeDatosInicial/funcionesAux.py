from Bio import SeqIO
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
# 




import networkx as nx

def armarMatrizDeFasta(archivo):
    matriz = []
    i=0
    for record in SeqIO.parse(archivo, "fasta"):
        line =str(record.seq)
        matriz.append(line)

    
    return np.array(matriz)

def calcularFrecuenciaColumna(matriz, columna, secuenciaRef):

    aminoacidoRef = secuenciaRef[columna]

    contador = 0
    i = 0
    nFilas = len(matriz)
    for i in range(nFilas):
        elemento = matriz[i]
 
        if ((elemento[columna] != aminoacidoRef) and (elemento[columna] !='X') and (elemento[columna] !='-') ):
            contador = contador + 1

    return contador/nFilas


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

def plotColumnas(columnas,FrecuenciaPorMes, pais):
    fig,ax = plt.subplots()

    for columna in columnas:
        plotColumna(ax,columna, FrecuenciaPorMes, pais)
    
    ax.set_xlabel("Meses")
    ax.set_ylabel("Indice de mutacion [%]")
    return(fig,ax)

def plotColumna(ax,columna, matrizFrecuencias, pais=None):
    frecuencias = [0]*matrizFrecuencias.shape[0]
    t=0
    for f in matrizFrecuencias:
        frecuencias[t] = 100*f[columna]
        t=t+1
    if pais == "Arg":
        meses =["5/20","6/20","7/20","8/20","9/20","10/20","11/20","12/20",
            "1/21","2/21","3/21","4/21","5/21","6/21","7/21","8/21","9/21","10/21","11/21","12/21",
            "1/22","2/22","3/22","4/22","5/22","6/22","7/22","8/22","9/22","10/22","11/22","12/22",
            "1/23","2/23","3/23","4/23","5/23","6/23","7/23","8/23","9/23","10/23","11/23","12/23"]
    else:
        meses = ["1/20","2/20","3/20","4/20","5/20","6/20","7/20","8/20","9/20","10/20","11/20","12/20"]
        
    markers = ['o', 's', 'D', '^', 'v', '>', '<', 'p', '*', 'h']
    np.random.shuffle(markers)
    ax.plot(meses,frecuencias)
    ax.scatter(meses,frecuencias, marker=markers[0],label=f"Residuo: {columna+1}")
    ax.legend()

def filtrarMatriz(matriz,file):
    #no tenia la secuencia de wohan.
    if (file.split('-')[0] == "Argentina"):
        posicionWohan={"Argentina-1-2021_Alineado.fasta": 1,
                      "Argentina-10-2020_Alineado.fasta": 1,
                      "Argentina-10-2021_Alineado.fasta": 1,
                      "Argentina-11-2020_Alineado.fasta": 1,
                      "Argentina-11-2021_Alineado.fasta": 3,
                      "Argentina-12-2020_Alineado.fasta": 1,
                      "Argentina-12-2021_Alineado.fasta": 4,
                      "Argentina-2-2021_Alineado.fasta": 8,
                      "Argentina-3-2021_Alineado.fasta": 1,
                      "Argentina-4-2021_Alineado.fasta": 1,
                      "Argentina-5-2020_Alineado.fasta": 3,
                      "Argentina-5-2021_Alineado.fasta": 6,
                      "Argentina-6-2020_Alineado.fasta": 16,
                      "Argentina-6-2021_Alineado.fasta": 0,
                      "Argentina-7-2020_Alineado.fasta": 1,
                      "Argentina-7-2021_Alineado.fasta": 1,
                      "Argentina-8-2020_Alineado.fasta": 1,
                      "Argentina-8-2021_Alineado.fasta": 14,
                      "Argentina-9-2020_Alineado.fasta": 5,
                      "Argentina-9-2021_Alineado.fasta": 1,
                      "Argentina-1-2022_Alineado.fasta": 3,
                      "Argentina-2-2022_Alineado.fasta": 9,
                      "Argentina-3-2022_Alineado.fasta": 6,
                      "Argentina-4-2022_Alineado.fasta": 5,
                      "Argentina-5-2022_Alineado.fasta": 1,
                      "Argentina-6-2022_Alineado.fasta": 2,
                      "Argentina-7-2022_Alineado.fasta": 4,
                      "Argentina-8-2022_Alineado.fasta": 2,
                      "Argentina-9-2022_Alineado.fasta": 5,
                      "Argentina-10-2022_Alineado.fasta": 4,
                      "Argentina-11-2022_Alineado.fasta": 4,
                      "Argentina-12-2022_Alineado.fasta": 23,
                      "Argentina-1-2023_Alineado.fasta": 5,
                      "Argentina-2-2023_Alineado.fasta": 1,
                      "Argentina-3-2023_Alineado.fasta": 3,
                      "Argentina-4-2023_Alineado.fasta": 2,
                      "Argentina-5-2023_Alineado.fasta": 6,
                      "Argentina-6-2023_Alineado.fasta": 1,
                      "Argentina-7-2023_Alineado.fasta": 1,
                      "Argentina-8-2023_Alineado.fasta": 6,
                      "Argentina-9-2023_Alineado.fasta": 3,
                      "Argentina-10-2023_Alineado.fasta": 2,
                      "Argentina-11-2023_Alineado.fasta": 1,
                      "Argentina-12-2023_Alineado.fasta": 1 
                      }
    elif (file.split('-')[0] == "England"):
        posicionWohan={                      
                        "England-1-2020_Alineado.fasta":1,
                        "England-2-2020_Alineado.fasta":1,
                        "England-3-2020_Alineado.fasta": 1,
                        "England-4-2020_Alineado.fasta":1,
                        "England-5-2020_Alineado.fasta":1,
                        "England-6-2020_Alineado.fasta":1,
                        "England-7-2020_Alineado.fasta":213 ,
                        "England-8-2020_Alineado.fasta": 380,
                        "England-9-2020_Alineado.fasta": 1,
                        "England-10-2020_Alineado.fasta": 1,
                        "England-11-2020_Alineado.fasta": 2,
                        "England-12-2020_Alineado.fasta": 2 
                      }
        
    
    fila = posicionWohan[file]
    nCol = 1274 #largo de la secuencia de wohan
    i=0
    while i< nCol:
        if matriz[fila][i] == "-":
            for j in range(len(matriz)):
                secuenciaFiltrada = matriz[j][:i] + matriz[j][i+1:]
                matriz[j]=secuenciaFiltrada
        else:
            i =i+1
    
        
           
    
#hay que analisar los caracteres x e - .

if __name__ == "__main__":
    
    carpetaArchivosAlineados = "archivosFasta/England/Alineados/"
    
    FrecuenciaPorMes = [[]]* len([f for f in os.listdir(carpetaArchivosAlineados) if f.endswith(".fasta")])# meses vamos a tomar
    
    for file in os.listdir(carpetaArchivosAlineados): 
        if file.endswith('.fasta'):
            print(f"procesando archivo {file} ...")
            
            archivo = open( carpetaArchivosAlineados + file, "r")
            
            matriz = armarMatrizDeFasta(archivo)
            
            filtrarMatriz(matriz, file)
            
            #secuencia de wohan
            secuenciaRef = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*"
           
            frecuenciasMes=calcularFrecuencias(archivo, matriz, secuenciaRef)
            
            mes = int(file.split('-')[1])
            if file.split('-')[0] == "Argentina":
               if int(file.split('-')[2].split('_')[0]) == 2021 :
                   FrecuenciaPorMes[mes+7] =frecuenciasMes
               elif int(file.split('-')[2].split('_')[0]) == 2020:
                   FrecuenciaPorMes[mes-5] =frecuenciasMes
               elif int(file.split('-')[2].split('_')[0]) == 2022 :
                   FrecuenciaPorMes[mes+19] =frecuenciasMes
               elif int(file.split('-')[2].split('_')[0]) == 2023 :
                   FrecuenciaPorMes[mes+31] =frecuenciasMes
            elif file.split('-')[0] == "England":
               FrecuenciaPorMes[mes-1] =frecuenciasMes
                
        
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
    
    
    





    G = nx.from_numpy_array(matrizAdjacencia)
    
    # Save the graph to a SIF file
    with open('graph.sif', 'w') as f:
        for edge in G.edges():
            f.write(f"{edge[0]} pp {edge[1]}\n")

    print("Graph saved as graph.sif") #para ver el grapho con cytoscape
    
    local_clustering = nx.clustering(G)
    
    nodosConCoeficienteClusteringAlto = {}

    # connected_nodes = set()
    # for component in nx.connected_components(G):
    #     connected_nodes.update(component)
    
    # # Crear un subgrafo con solo los nodos conectados
    # G_connected = G.subgraph(connected_nodes).copy()

    # # Dibujar el subgrafo con nodos conectados
    # plt.figure(figsize=(8, 6))
    # nx.draw(G_connected, with_labels=True, node_color='lightblue', edge_color='gray')
    # plt.title('Grafo con Solo Nodos Conectados')
    # plt.show()

    # # Connect to Cytoscape
    # # client = CyRestClient()


            
        
                    
            
            
    

    
    
    