from Bio import GenBank
import csv

"""
Format for genbank files...

LOCUS       13C210_ACLSV      7541 bp    DNA     linear   UNA 
SOURCE      No definition (unknown)
FEATURES             Location/Qualifiers
     CDS             150..5801  
                     /cds_type="ORF1 replicase"
                     /note="Length: 5652"
                     /note="Found at strand: positive"
                     /note="Start codon: ATG"
     CDS             5713..7095
                     /cds_type="ORF2 movement protein"
                     /note="Length: 1383"
                     /note="Found at strand: positive"
                     /note="Start codon: ATG"
     CDS             6608..7360
                     /cds_type="ORF3 coat protein"
                     /note="Length: 753"
                     /note="Found at strand: positive"
                     /note="Start codon: ATG"
ORIGIN
        1 ATACTGATAC AGTGTACACT CACGTCGTGA GTAAACAGAT TGACGTAACG CCTCAATCGT....

WORKFLOW:
-store locus for each record
-store all locations + location data for each record
-output with cvs
"""

with open("ACLSV.gbk") as handle:

    Features = []

    for record in GenBank.parse(handle):
        #print(record.locus)
        for feature in record.features:
            Features.append(feature)
            #print(feature.qualifiers[0].value)
            #print("\n")



locations = []
for cds in Features:
    locations.append(cds.location)

location_tuples = []
for loc in locations:
    loc0 = ''
    loc1 = ''
    pos = 0
    for char in loc:
        if char != '.':
            loc0 += char
            pos += 1
        else:
            break

    for char in range(pos+2, len(loc)):
        loc1 += str(loc[char])

    temp_tuple = (loc0, loc1)
    location_tuples.append(temp_tuple)


print(location_tuples)


























