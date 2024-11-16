from Bio import SeqIO

# Ruta al archivo FASTA
fasta_file = "spikeprot1103.fasta"

lista = []

i=0
limite = 100000
archivoAux = open("ArgentinaAhoraSi.fasta","x")
limiteCaracteresID = 0
identifiers = {}

pais="Argentina"

# Leer el archivo FASTA y las secuecuencias con ID similar les agrego un numero
with open(fasta_file, "r",encoding="latin-1") as file:
    for record in SeqIO.parse(file, "fasta"):
        if(len(str(record.id))<limiteCaracteresID):
            pass
        else:
            if (record.id.split('/')[1] == pais):
                if record.id in identifiers:
                    # Incrementar el contador de duplicados
                    identifiers[record.id] += 1

                else:
                     identifiers[record.id] = 0
                     archivoAux.write(">"+str(record.id)+"\n" +str(record.seq)+"\n")
                i = i+1
            
            if (i>limite):
                break
print(f"hubo {i} secuencias duplicadas")

#Separar paor pais
#Poner la secuencia de wohan
#Eliminar las secuencias repetidas

