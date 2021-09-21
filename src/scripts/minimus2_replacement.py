import os
with open('minimus2_location.txt', 'r') as f:
    lineas = f.readlines()
    path_to_read = lineas[0]

with open('minimus2', 'w') as new:
    pass

with open('minimus2', 'a') as new:
    with open(path_to_read[0:-1], 'r') as f:
        lineas = f.readlines()
        i = 0
        while i < 45:
            new.write(lineas[i])
            i = i+1
        correct_line = lineas[43]
    correct_path = correct_line.split("=")
    path = correct_path[1]
    path = path[0:-1]
    new.write('DELTAFILTER = ' + path + '/delta-filter' + '\n')
    new.write('SHOWCOORDS = ' + path + '/show-coords' + '\n')

with open('minimus2', 'a') as new:
    with open(path_to_read[0:-1], 'r') as f:
        lineas = f.readlines()
        i = 47
        while i > 46 and i < 57:
            new.write(lineas[i])
            i = i+1
    new.write('20: $(NUCMER) --maxmatch -c $(OVERLAP) $(REFSEQ) $(QRYSEQ) -p $(PREFIX)' + '\n')

with open('minimus2', 'a') as new:
    with open(path_to_read[0:-1], 'r') as f:
        lineas = f.readlines()
        i = 58
        while i > 57 and i < 85:
            new.write(lineas[i])
            i = i+1

os.chmod("minimus2", 0o777)
os.replace("minimus2", path_to_read[0:-1])
