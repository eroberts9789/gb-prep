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

def format_locs(location):
    loc0 = ''
    loc1 = ''
    pos = 0
    for char in location:
        if char != '.':
            loc0 += char
            pos += 1
        else:
            break

    for char in range(pos + 2, len(location)):
        loc1 += str(location[char])

    loc_tuple = (loc0, loc1)
    return loc_tuple

def format_qualifiers(qualifier):
    quals0 = ''
    quals1 = ''
    pos = 0
    for char in qualifier:
        if char != ' ':
            quals0 += char
            pos += 1
        else:
            break

    for char in range(pos +1, len(qualifier)):
        quals1 += qualifier[char]
    if quals1 == '':
        quals1 = quals0
    quals0 = str.strip(quals0)
    quals1 = str.strip(quals1)
    print(quals0 + quals1)
    quals_tuple = (quals0, quals1)
    return quals_tuple


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
            #print(feature.qualifiers[0].value)
            format_qualifiers(feature.qualifiers[0].value)

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

with open(output_name, 'w', newline='') as file:
    for record in record_list:
        writer = csv.writer(file, delimiter='\t')

        header_vals = ['>Feature', str(record["locus"])]
        writer.writerow(header_vals)

        for feature in record["features"]:
            #print line: loc0	loc1	feature_key
            locs = format_locs(feature.location)

            loc_line_vals = [str(locs[0]), str(locs[1]), str(feature.key)]
            writer.writerow(loc_line_vals)

            #print lines: product	<first word in qualifier desc>
            qualifiers = format_qualifiers(feature.qualifiers[0].value)

            product_line_list = [['Product', qualifiers[0]], ['Product', qualifiers[1]]]
            writer.writerows(product_line_list)





































