import os
from Bio import SeqIO

fileslist = os.listdir("results/binning/metabat/bins/")

with open("results/binning/metabat/bins/metabat.bins.txt", "a") as outfile:
    outfile.write("scaffold" + "\t" + "bin")
    for i in fileslist:
        bin_seqs = SeqIO.parse(open(os.path.join("bins/", i)),'fasta')
        for fasta in bin_seqs:
            name, sequence = fasta.id, str(fasta.seq)
            outfile.write('\n' + name + "\t" + i.replace('.', '')[:-2])
