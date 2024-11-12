from Bio import SeqIO

# Ruta al archivo FASTA
fasta_file = "spikeprot1103.fasta"

lista = []

i=0
limite = 400
archivo = open("archivo.fasta","x")
limiteCaracteresID = 0
# Leer el archivo FASTA
with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        if(len(str(record.id))<limiteCaracteresID):
            pass
        else:
            archivo.write(">"+str(record.id)+"\n" +str(record.seq)+"\n")
            i = i+1
        if (i>limite):
            break
