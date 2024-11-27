import numpy as np
import os
import funcionesAux as mf
import networkx as nx
# 


def obtenerGrafo(carpetaArchivosAlineados, theshold):
    
    FrecuenciaPorMes = [[]]*20#20 meses vamos a tomar
    
    for file in os.listdir(carpetaArchivosAlineados): 
        if file.endswith('.fasta'):
            print(f"procesando archivo {file} ...")
            
            archivo = open( carpetaArchivosAlineados + file, "r")
            
            matriz = mf.armarMatrizDeFasta(archivo)
            
            mf.filtrarMatriz(matriz, file)
            
            #secuencia de wohan
            secuenciaRef = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*"
           
            frecuenciasMes=mf.calcularFrecuencias(archivo, matriz, secuenciaRef)
            
            mes = int(file.split('-')[1])
            if int(file.split('-')[2].split('_')[0]) == 2021 :
                FrecuenciaPorMes[mes+7] =frecuenciasMes
            else:
                FrecuenciaPorMes[mes-5] =frecuenciasMes
        
    FrecuenciaPorMes = np.array(FrecuenciaPorMes)
        
    matrizAdjacencia = mf.calcularCorelacion(FrecuenciaPorMes)
        
    matrizAdjacencia[np.isnan(matrizAdjacencia)] = 0
    matrizAdjacencia[matrizAdjacencia<theshold]=0
    matrizAdjacencia[matrizAdjacencia>=theshold]=1
    
    G = nx.from_numpy_array(matrizAdjacencia)
    
    # Save the graph to a SIF file
    with open('graph.sif', 'w') as f:
        for edge in G.edges():
            f.write(f"{edge[0]} pp {edge[1]}\n")

    print("Graph saved as graph.sif") #para ver el grapho con cytoscape
    
    return G,matrizAdjacencia,FrecuenciaPorMes

if __name__ == "__main__":
    carpetaArchivosAlineados = "archivosFasta/argentina/Alineados/"
    sol = obtenerGrafo(carpetaArchivosAlineados, theshold=0.99)
    
    # mf.plotColumnas([245, 74,75, 858,189,17,25,489,19,1026,416,500,137,1175,500,416,654,483], sol[2])
    
    
    
    
    
    