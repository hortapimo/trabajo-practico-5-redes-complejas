import numpy as np
import os
import funcionesAux as mf
import networkx as nx
# 


def obtenerGrafo(carpetaArchivosAlineados, theshold, thresholdFrec = 0, nameGraph = "default"):
    
   
    FrecuenciaPorMes = [[]]* len([f for f in os.listdir(carpetaArchivosAlineados) if f.endswith(".fasta")])# meses vamos a tomar
    cantSecuenciasProcesadas = 0
    
    for file in os.listdir(carpetaArchivosAlineados): 
        if file.endswith('.fasta'):
            print(f"procesando archivo {file} ...")
            
            archivo = open( carpetaArchivosAlineados + file, "r")
            
            matriz = mf.armarMatrizDeFasta(archivo)
            
            print(f"cantidad de secuencias procesadas : {len(matriz)}")
            cantSecuenciasProcesadas = cantSecuenciasProcesadas +len(matriz)
            mf.filtrarMatriz(matriz, file)
            
            #secuencia de wohan
            secuenciaRef = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*"
           
            frecuenciasMes=mf.calcularFrecuencias(archivo, matriz, secuenciaRef)
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
        
    print(f"cantidad total secuencias procesadas : {cantSecuenciasProcesadas}")
                
                
        
    FrecuenciaPorMes = np.array(FrecuenciaPorMes)
    FrecuenciaPorMes[FrecuenciaPorMes < thresholdFrec] = 0
    
    matrizAdjacencia = mf.calcularCorelacion(FrecuenciaPorMes)
        
    matrizAdjacencia[np.isnan(matrizAdjacencia)] = 0
    matrizAdjacencia[matrizAdjacencia<theshold]=0
    matrizAdjacencia[matrizAdjacencia>=theshold]=1
    
    G = nx.from_numpy_array(matrizAdjacencia)
    
    # Save the graph to a SIF file
    with open(f"{nameGraph}.sif", 'w') as f:
        for edge in G.edges():
            f.write(f"{edge[0]} pp {edge[1]}\n")

    print(f"Graph saved as {nameGraph}.sif") #para ver el grapho con cytoscape
    
    return G,matrizAdjacencia,FrecuenciaPorMes

if __name__ == "__main__":
    carpetaArchivosAlineados = "archivosFasta/England/Alineados/"
    sol = obtenerGrafo(carpetaArchivosAlineados, theshold=0.89, thresholdFrec = 0.05, nameGraph='England_09')
    
    # mf.plotColumnas([245, 74,75, 858,189,17,25,489,19,1026,416,500,137,1175,500,416,654,483], sol[2])
    
    
    
    
    
    