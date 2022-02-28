#!/usr/bin/python
### When chosen assembly is "per sample"

project = str(input("Name of your project (leave blank to use default): "))
if project == "":
    project_line = "\"project\" : " +"\"" + "My_project" +"\","
else:
    project_line = "\"project\" : " +"\"" + project +"\","

host_removal = str(input("Do you want to remove host DNA from the analysis? (yes/no): "))
if host_removal == "yes":
    host_removal_line = "\"host_removal\" : " + "True" + ","
elif host_removal == "no":
    host_removal_line = "\"host_removal\" : " + "False" + ","
else:
    print("Error: unvalid option, please type'yes' or 'no'")

reference = str(input("Do you have your own reference genome? (yes/no): "))
if reference == "no":
    reference_line = "\"reference\" : " +"\"" + "data/test-data/reference/beta-vulgaris-subset.fasta" +"\","
    reference_index_line = "\"reference_index\" : " +"\"" + "/data/shared/genomes/Eukaryotes/Phaseolus_vulgaris/GCA_000499845.1_PhaVulg1_0_genomic" +"\","
elif reference == "yes":
    own_ref = str(input("Path to your reference genome: "))
    own_ref_index = str(input("Path to the index file of your reference genome: "))
    reference_line = "\"reference\" : " +"\"" + own_ref +"\","
    reference_index_line = "\"reference_index\" : " +"\"" + own_ref_index +"\","

assembler = str(input("Choose an assembler (megahit/spades) (Default is megahit): "))
if assembler == "megahit" or assembler == "":
    assembler_line = "\"assembler\" : " +"\"" + "megahit" +"\","
elif assembler == "spades":
    assembler_line = "\"assembler\" : " +"\"" + assembler +"\","
else:
    print("Error: please type a valid option")

x_0 = "{{{0} " ## Open dict
x_1 = "}}{0}" ## Close dict
tmpdir_line = "\"tmpdir\" : " +"\"" + "/tmp" +"\","
adapters_line = "\"adapters\" : " +"\"" + "/data/shared/tools/Trimmomatic/0.36/adapters/NexteraPE-PE.fa" +"\","
min_qual_line = "\"min_qual\" : " +"\"" + "30" +"\","
min_length_line = "\"min_length\" : " +"\"" + "150" +"\","
emapper_database_line = "\"emapper_database\" : " +"\"" + "/data/shared/db/eggnogdb/5.0.0" +"\","
emapper_diamond_line = "\"emapper_diamond\" : " +"\"" + "/data/shared/db/eggnogdb/5.0.0/eggnog_proteins.dmnd" +"\","
CAT_database_line = "\"CAT_database\" : " +"\"" + "/data/shared/db/CAT/CAT_prepare_20200618/2020-06-18_CAT_database/" +"\","
CAT_taxonomy_line = "\"CAT_taxonomy\" : " +"\"" + "/data/shared/db/CAT/CAT_prepare_20200618/2020-06-18_taxonomy/" +"\","
filter_contigs_length_line = "\"filter_contigs_length\" : " +"\"" + "2000" +"\","
filter_contigs_antismash_line = "\"filter_contigs_antismash\" : " +"\"" + "5000" +"\","
binner_line = '\"binner\" : ' + '["metabat", "vamb"]' + ','
big_line = "\"big\" : " +"\"" + "bigscape" +"\","
kmers_line = "\"kmers\" : " +"\"" + "33,55,77,99,127" +"\","
assembly_klist_line = "\"assembly-klist\" : "
meta_large_line = "\"meta-large\" : " +"\"" + "27,37,47,57,67,77,87,97,107,117,127" +"\","
data_line = "\"data\" : "
treatment_line = "\"treatment\" : "

import os
data_path = str(input("Where are your input files? (Please follow the format /path/) (Leave blank to use default): "))
if data_path == "":
    folder = 'data/' 
else:
    folder = data_path


import os
list_with_treatments = list()
for i in os.listdir(folder):
    if os.path.isdir(os.path.join(folder, i)):
        list_with_treatments.append(i)

with open("config.json", "w") as file:
    file.write('{\n')
    file.write('    ' + project_line + '\n')
    file.write('    ' + tmpdir_line + '\n')
    file.write('    ' + adapters_line + '\n')
    file.write('    ' + min_qual_line + '\n')
    file.write('    ' + min_length_line + '\n')
    file.write('    ' + host_removal_line + '\n')
    file.write('    ' + reference_line + '\n')
    file.write('    ' + reference_index_line + '\n')
    file.write('    ' + emapper_database_line + '\n')
    file.write('    ' + emapper_diamond_line + '\n')
    file.write('    ' + CAT_database_line + '\n')
    file.write('    ' + CAT_taxonomy_line + '\n')
    file.write('    ' + filter_contigs_length_line + '\n')
    file.write('    ' + filter_contigs_antismash_line + '\n')
    file.write('    ' + assembler_line + '\n')
    file.write('    ' + binner_line + '\n')
    file.write('    ' + big_line + '\n')
    file.write('    ' + kmers_line + '\n')
    file.write('    ' + assembly_klist_line + x_0.format("") + '\n')
    file.write('    ' + '   ' + meta_large_line + '\n')
    file.write('     ' + x_1.format("") + "," '\n')
    file.write('     ' + treatment_line + x_0.format("") + '\n')
    for x in list_with_treatments:
        files = os.listdir(os.path.join(folder, x))
        files.sort()
        forward_list = list()
        reverse_list = list()
        for i in files:
            if "_R1_" in i:
                forward_list.append(i)
            elif "_R2_" in i:
                reverse_list.append(i)

        lista_para_nombresdelinea = files
        main_name_list = list()
        for i in lista_para_nombresdelinea:
            split_string = i.split(".")
            main_name = split_string[0]
            main_name_2 = main_name.split("_R")
            if len(main_name_2) == 2:
                main_name_list.append(main_name_2[0])
            else:
                print("An unvalid name has been found between your files. Please rename and execute again")
        main_name_list = list(dict.fromkeys(main_name_list))
        if list_with_treatments.index(x) == len(list_with_treatments)-1:
            file.write('    ' + '   "' + x + '" : ' + str(main_name_list) + '\n')
        else:
            file.write('    ' + '   "' + x + '" : ' + str(main_name_list) + ',' + '\n')
    file.write('     ' + x_1.format("") + "," + '\n')
    file.write('     ' + data_line + x_0.format("") + '\n')
    for x in list_with_treatments:
        files = os.listdir(os.path.join(folder, x,))
        files_path = os.path.join(folder, x)
        files_path = files_path + "/"
        files.sort()
        forward_list = list()
        reverse_list = list()
        for i in files:
            if "_R1_" in i:
                forward_list.append(i)
            elif "_R2_" in i:
                reverse_list.append(i)

        lista_para_nombresdelinea = files
        main_name_list = list()
        for i in lista_para_nombresdelinea:
            split_string = i.split(".")
            main_name = split_string[0]
            main_name_2 = main_name.split("_R")
            if len(main_name_2) == 2:
                main_name_list.append(main_name_2[0])
            else:
                print("An unvalid name has been found between your files. Please rename and execute again")
        main_name_list = list(dict.fromkeys(main_name_list))

        for i in range(len(main_name_list)):
            if i == len(main_name_list)-1 and list_with_treatments.index(x) == len(list_with_treatments)-1:
                file.write('    ' + '   ' + "\"" + main_name_list[i] + '": {"forward": ' + "\"" + str(files_path) + forward_list[i] + "\"" + ', "rev": ' + "\"" + str(files_path) + reverse_list[i] + "\"}" + '\n')
            else:
                file.write('    ' + '   ' + "\"" + main_name_list[i] + '": {"forward": ' + "\"" + str(files_path) + forward_list[i] + "\"" + ', "rev": ' + "\"" + str(files_path) + reverse_list[i] + "\"}," + '\n')
    file.write('     ' + x_1.format("") + '\n')
    file.write(x_1.format("") + '\n')
    print("Configuration file successfully created. Chosen assembly mode: pooled assembly")
