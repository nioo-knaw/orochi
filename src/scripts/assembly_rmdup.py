with open('scratch/assembly/megahit/minimus2/all.merged.contigs.fasta') as f:
    lines = f.readlines()

rango = list(range(0, int(round(len(lines)/2)), 2))
headers = list()
seqs = list()
for i in rango:
    headers.append(lines[i])
    seqs.append(lines[i+1])

filtered_headers = list()
for i in range(len(headers)):
    prueba = headers[i].split(" ")
    filtered_headers.append(prueba[0])

setOfElems = set()
filtered_seqs = list()
filtered_complete_headers = list()
repetidos_descartados = list()
for elem in range(len(filtered_headers)):
    if filtered_headers[elem] in setOfElems:
        repetidos_descartados.append(filtered_headers[elem])
    else:
        setOfElems.add(filtered_headers[elem]) 
        filtered_seqs.append(seqs[elem])
        filtered_complete_headers.append(headers[elem])

with open('scratch/assembly/megahit/minimus2/filtered_assembly.fasta', "a") as new_file:
    for i in range(len(setOfElems)):
        new_file.write(filtered_complete_headers[i] + '\n')
        new_file.write(filtered_seqs[i] + '\n')
