from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

# Ruta al ejecutable de Muscle y al archivo FASTA de entrada
muscle_exe = r"muscle.exe"
fasta_file = r"archivo4.fasta"
output_file = r"archivo4_aligned.fasta"

# Correr Muscle
cline = MuscleCommandline(muscle_exe, input=fasta_file, out=output_file)
cline()

# Leer y mostrar la alineación resultante
alignment = AlignIO.read(output_file, "fasta")
print(alignment)
