from Bio import GenBank
import csv

"""
FORMAT FOR INPUT/GENBANK FILE WHICH WE WANT TO EXTRACT DESIRED DATA FROM

LOCUS       13C202_PNRSV_RNA2      2572 bp    DNA     linear   UNA 
FEATURES             Location/Qualifiers
     CDS             18..2417
                     /cds_type="ORF"
                     /product="polymerase p2"
                     /note="Length: 2400"
                     /note="Found at strand: positive"
                     /note="Start codon: ATG"
ORIGIN
        1 CTCGTGGTTG AGTTACAATG AATCCTTTGA ATCAGTTGTT GACACATGGG TGTACTGCTA
			

HOW DATA IS STORED IN THIS FILE BEFORE OUTPUT:

we loop through all records in genbank file and create a list of record dictionaries named record_list that will conatain all info we want to output
record dictionaries will be formatted as:

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


input_name = "PNRSV_RNA1.gbk"
output_name = "ACLSVfeatures.txt"

record_list = []

with open(input_name) as handle:

    first_record = True
    first_format = {
        'feature_keys' : [],
        'feature_type' : '',
        'product' : ''
    }

    for record in GenBank.parse(handle):

        record_info = {
            "locus": "",
            "features": []
        }

        record_info["locus"] = record.locus

        for feature in record.features:

            #only do this for first record, this is the record that is best formatted and the rest of the feature qualifiers should match this format in the output
            if first_record == True:
                first_format['feature_keys'].append(feature.key)
                for qualifier in feature.qualifiers:
                    if "type" in qualifier.key:
                        first_format['feature_type'] = qualifier.value
                        #print(first_format['feature_type'])
                    if "product" in qualifier.key:
                        first_format['product'] = qualifier.value
                        #print(first_format['product'])
            first_record = False

            #only save location output for rest of records because features output will match first record
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

        header_vals = ['>Feature', record["locus"]]
        writer.writerow(header_vals)

        for feature in record["features"]:
            if feature.key in first_format['feature_keys']:
                locs = format_locs(feature.location)

                loc_line_vals = [locs[0], locs[1], str(feature.key)]
                writer.writerow(loc_line_vals)

                feature_type = str.strip(first_format['feature_type'], '"')
                product = str.strip(first_format['product'], '"')

                if len(feature_type) > 0:
                    type_line_format = ['\t', '\t', 'Product', feature_type]
                    writer.writerow(type_line_format)
                if len(product) > 0:
                    product_line_format = ['\t', '\t', 'Product', product]
                    writer.writerow(product_line_format)

        writer.writerow('\t')














































