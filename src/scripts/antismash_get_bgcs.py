import json

def get_location(location):
    ## https://github.com/ohmeta/oral-assembly/blob/master/notebook/antismash.ipynb
    loc_ = location.split(":")
    return int(loc_[0].lstrip("[")), int(loc_[1].rstrip("]"))

json_dict = json.load(open(snakemake.input[0]))

outfile = open(snakemake.output[0],"w")

for record in json_dict['records']:
    for feat in record['features']:
            if feat['type'] == "region":
                start,stop = get_location(feat['location'])
                #print(record['seq']['data'][start:stop])
                #print(feat['qualifiers']['product'])
                #print(record['id'])
                outfile.write(">%s | %s\n%s\n" % (record['id'], feat['qualifiers']['product'][0], record['seq']['data'][start:stop]))
outfile.close()
