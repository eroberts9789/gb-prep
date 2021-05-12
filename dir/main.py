from Bio import GenBank
import csv

"""
FORMAT FOR INPUT/GENBANK FILE WHICH WE WANT TO EXTRACT DESIRED DATA FROM

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
			

HOW DATA IS STORED IN THIS FILE BEFORE OUTPUT:

we loop through all records in genbank file and create a list of record dictionaries named record_list that will conatain all info we want to output
record dictionaries will be formatted as:

record = { "locus" : "asdf1234",
           "features" : [     CDS             150..5801
                                             /cds_type="ORF1 replicase"
                                             /note="Length: 5652"
                                             /note="Found at strand: positive"
                                             /note="Start codon: ATG",
                                             
                              CDS             150..5801
                                             /cds_type="ORF1 replicase"
                                             /note="Length: 5652"
                                             /note="Found at strand: positive"
                                             /note="Start codon: ATG",
                                             
                              CDS             150..5801
                                             /cds_type="ORF1 replicase"
                                             /note="Length: 5652"
                                             /note="Found at strand: positive"
                                             /note="Start codon: ATG"
                                             ]
        }
"""
input_name = "ACLSV.gbk"
output_name = "ACLSVfeatures.txt"

record_list = []

with open(input_name) as handle:

    for record in GenBank.parse(handle):
        Feature = {
            "locus": "",
            "features": [],
        }

        Feature["locus"] = record.locus

        for feature in record.features:
            Feature["features"].append(feature)

        #print(Feature)
        record_list.append(Feature)

"""
NOW WE WRITE OUTPUT USING CSV IN THE FORMAT SHOWN BELOW:

>Feature	<LOCUS>
loc0	loc1	feature_key
                product	<first word in qualifier desc>
                product	<other words in qualidier desc>
loc0	loc1	feature_key
                product	<first word in qualifier desc>
                product	<other words in qualidier desc>
loc0	loc1	feature_key
                product	<first word in qualifier desc>
                product	<other words in qualidier desc>
                
lOOP BACK THROUGH, PRINTING AS WE GO!
"""
with open ('output_name', 'w', newline='') as csvfile:
    for record in record_list:
        pass




"""
#PROCEDURE FOR GETTING LOCATION TUPLES FROM 123..1234 FORMAT FOR LOCATIONS
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

"""
























