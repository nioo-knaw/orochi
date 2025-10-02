# Functions to filter annotations so only contigs that are prokaryotic are kept


def filter_gff(gff, headerfile, outdir, size, sample_name):
    """
    Filter gff based on contig_headers.txt file.
    """

    header_file = open(headerfile, "r")
    contig_headers = header_file.read()
    contig_headers_list = contig_headers.split("\n")  # list of all headers
    header_file.close()

    gff_file = sample_name + "_" + "prokaryote" + "_" + str(size) + ".gff"

    with open(gff, "r") as f:
        with open(os.path.join(outdir, gff_file), "w") as g:
            for line in f:
                if line.startswith("#"):
                    pass

                else:
                    contig_name = line.split("\t")[0]
                    if len([x for x in contig_headers_list if contig_name in x]) > 0:
                        g.write(line)

if __name__ == "__main__":
    # Parse arguments from Snakemake
    gff = sys.argv[2]
    outdir = sys.argv[4]
    size = int(sys.argv[6])
    header_file = sys.argv[8]
    sample_name = sys.argv[10]

    size_filter(gff, outdir, size, header_file, sample_name)

