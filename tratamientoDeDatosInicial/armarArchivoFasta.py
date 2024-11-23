from Bio import SeqIO

# Ruta al archivo FASTA
fasta_file = "spikeprot1103.fasta"

lista = []

i=0

pais="Brasil"
anio = 2021
meses = [1,2,3,4,5,6,7,8,9,10,11,12]
for mes in meses:    
    archivoAux = open(f"{pais}-{mes}-{anio}.fasta","x")
    minSeparador = 4
    identifiers = {}
    
    
    # Leer el archivo FASTA y las secuecuencias con ID similar les agrego un numero
    with open(fasta_file, "r",encoding="latin-1") as file:
        for record in SeqIO.parse(file, "fasta"):
            if(len(record.id.split('/'))<minSeparador):
                pass
            else:
                idRecord = record.id
                paisRecord = idRecord.split('/')[1]
                
                try:
                    try:
                        anioRecord = int(idRecord.split('/')[3].split('|')[0])
                        mesRecord= int(idRecord.split('/')[3].split('|')[1].split('-')[1])
                    except:
                        try:
                            anioRecord = int(idRecord.split('/')[4].split('|')[0])
                            mesRecord= int(idRecord.split('/')[4].split('|')[1].split('-')[1])
                        except:
                            anioRecord = int(idRecord.split('/')[5].split('|')[0])
                            mesRecord= int(idRecord.split('/')[5].split('|')[1].split('-')[1])
                except:
                    anioRecord = None
                    mesRecord= None
                    
     
                
                if ((paisRecord == pais) and(mesRecord == mes) and (anioRecord == anio)):
    
                        if record.id in identifiers:
                            # Incrementar el contador de duplicados
                            identifiers[record.id] += 1
        
                        else:
                             identifiers[record.id] = 0
                             archivoAux.write(">"+str(record.id)+"\n" +str(record.seq)+"\n")
                        i = i+1
    
    print(f"tiene {i} secuencias") 

#Separar paor pais
#Poner la secuencia de wohan
#Eliminar las secuencias repetidas

