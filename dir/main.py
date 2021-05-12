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

    for x in range(pos + 2, len(location)):
        loc1 += str(location[x])

    locs = [loc0, loc1]
    #print (loc0 + " " +loc1)
    return locs

def format_quals(qualifier):
    qualifier = str.strip(qualifier, '"')
    qual0 = ''
    qual1 = ''
    pos = 0
    for char in qualifier:
        if char != ' ':
            qual0 += char
            pos += 1
        else:
            break

    for x in range(pos +1, len(qualifier)):
        qual1 += qualifier[x]

    return [qual0, qual1]


input_name = "ACLSV.gbk"
output_name = "ACLSVfeatures.txt"

record_list = []

with open(input_name) as handle:

    for record in GenBank.parse(handle):
        record_info = {
            "locus": "",
            "features": []
        }

        record_info["locus"] = record.locus

        for feature in record.features:
            record_info["features"].append(feature)

        record_list.append(record_info)

"""
NOW WE WRITE OUTPUT USING CSV IN THE FORMAT SHOWN BELOW:

>Feature	<LOCUS>
loc0	loc1	feature_key
                product	<first word in qualifier desc>
                product	<other words in qualifier desc>
loc0	loc1	feature_key
                product	<first word in qualifier desc>
                product	<other words in qualifier desc>
loc0	loc1	feature_key
                product	<first word in qualifier desc>
                product	<other words in qualifier desc>
                
lOOP BACK THROUGH, PRINTING AS WE GO!
"""

with open(output_name, 'w', newline='') as output_file:
    writer = csv.writer(output_file, delimiter='\t', escapechar=' ', quoting=csv.QUOTE_NONE)

    for record in record_list:
        header_vals = ['>Feature', str(record["locus"])]
        writer.writerow(header_vals)

        for feature in record["features"]:
            """
            first change format '124..5675' given in feature.location to a list of locations stored in locs
            the format output so it prints line "loc0	loc1	feature_key"
            """
            locs = format_locs(feature.location)

            loc_line_vals = [locs[0], locs[1], str(feature.key)]
            writer.writerow(loc_line_vals)

            """
            first seperate the feature type into first word and following description (if exists, don't output if following description doens't exist)
            """
            quals = format_quals(feature.qualifiers[0].value)

            for qual in quals:
                if len(qual) != 0:
                    product_line_vals = ['\t','\t', 'Product', qual]
                    writer.writerow(product_line_vals)
        writer.writerow('\t')










































