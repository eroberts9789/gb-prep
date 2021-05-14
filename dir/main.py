
from Bio import GenBank
import csv

"""
FORMAT FOR INPUT/GENBANK FILE WHICH WE WANT TO EXTRACT DESIRED DATA FROM

LOCUS       13C214_PNRSV_RNA3      1943 bp    RNA     linear   VRL 03-MAY-2012
FEATURES             Location/Qualifiers
     CDS             174..1025
                     /cds_type="ORF"
                     /product="movement protein"
                     /note="Length: 852"
                     /note="Found at strand: positive"
                     /note="Start codon: ATG"
     CDS             1100..1774
                     /cds_type="ORF"
                     /product="coat protein"
                     /note="Length: 675"
                     /note="Found at strand: positive"
                     /note="Start codon: ATG"
ORIGIN
        1 GTTTTTACAA TTGAAATCGT AACATCCAGC AATTGGTTGA TTTCACTTTT GCTTGACTAG


HOW DATA IS STORED IN THIS FILE BEFORE OUTPUT:

we loop through all records in genbank file and create a list of record dictionaries named record_list that will conatain all info we want to output
record dictionaries will be formatted as:

"""
global input_name
input_name = "PNRSV_RNA3.gbk"
global output_name
output_name = "features_output.txt"
formatted_right = False

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
    # print (loc0 + " " +loc1)
    return locs

def format_record_list(first_record, record_list):
    first_record_feature_keys = []
    for feature in first_record['features']:
        first_record_feature_keys.append(feature.key)

    must_fix_format = []
    result = []
    for record in record_list:
        """
        REMOVE ANY REPEATED FEATURES BY TURNING INTO DICT AND THEN BACK TO LIST
        """
        record['features'] = list(dict.fromkeys(record['features']))

        for feature in record['features']:
            """
            IF THERE IS A FEATURE KEY THAT DOESN'T BELONG IN THE OUTPUT, REMOVE IT FROM THE RECORD IN THE RECORD_LIST
            """
            if feature.key not in first_record_feature_keys:
                record['features'].remove(feature)

        if len(record['features']) != len(first_record_feature_keys):
            must_fix_format.append(record)

        if record not in must_fix_format:
            result.append(record)

    global formatted_right
    if len(must_fix_format) > 0:
        formatted_right = False
        print('SOME OF YOUR RECORDS ARE FORMATTED INCORRECTLY. THE FOLLOWING ARE NOT INCLUDED IN YOUR OUTPUT FILES.')
        for record in must_fix_format:
            print(record['locus'])
    else: formatted_right = True

    return (result)



"""
PARSE GENBANK FILE AND STORE ALL REQUIRED INFO IN RECORD DICTIONARIES IN RECORD_LIST LIST
"""

with open(input_name) as handle:
    record_list = []
    for gbk_record in GenBank.parse(handle):

        record = {
            "locus": "",
            "features": []
        }

        record["locus"] = gbk_record.locus

        for feature in gbk_record.features:
            record["features"].append(feature)

        record_list.append(record)

"""
update records in records list so that they all match the format of the first record
"""
first_record = record_list[0]
all_recs = format_record_list(first_record, record_list)


"""
NOW WRITE EVERYTHING TO FEATURES FILE
"""

with open(output_name, 'w', newline='') as output_file:
    writer = csv.writer(output_file, delimiter='\t', escapechar=' ', quoting=csv.QUOTE_NONE)
    if formatted_right == True:
        for record in record_list:
            """
            write header line in the format:
            >Feature	<LOCUS>
            """
            header_format = ['>Feature', record["locus"]]
            writer.writerow(header_format)

            for feature in record["features"]:
                """
                write location line in the format:
                <loc0>	<loc1>	<feature.key>
                """
                locs = format_locs(feature.location)
                loc_line_format = [locs[0], locs[1], str(feature.key)]
                writer.writerow(loc_line_format)

                for qualifier in feature.qualifiers:
                    """
                    looks for qualifier.key named "feature_type". If found, print line in the format:
                                    Product	<feature_type>
                    """
                    if 'type' in qualifier.key:
                        type = str.strip(qualifier.value, '"')
                        type_line_format = ['\t', '\t', 'Product', type]
                        writer.writerow(type_line_format)
                    """
                    looks for qualifier.key named "product". If found, print line in the format:
                                    Product	<product>
                    """
                    if 'product' in qualifier.key:
                        product = str.strip(qualifier.value, '"')
                        product_line_format = ['\t', '\t', 'Product', product]
                        writer.writerow(product_line_format)
                    """
                    IF OTHER QUALIFIERS ARE NEEDED THEN FORMAT LIKE ABOVE
                    """
            #newline
            writer.writerow('\t')

    else:
        for record in record_list:
            """
            write header line in the format:
            >Feature	<LOCUS>
            """
            header_format = ['>Feature', record["locus"]]
            writer.writerow(header_format)
            for feature, first_feature in zip(record["features"], first_record['features']):
                """
                write location line in the format:
                <loc0>	<loc1>	<feature.key>
                """
                locs = format_locs(feature.location)
                loc_line_format = [locs[0], locs[1], str(feature.key)]
                writer.writerow(loc_line_format)

                for qualifier in first_feature.qualifiers:
                    """
                    looks for qualifier.key named "feature_type". If found, print line in the format:
                                    Product	<feature_type>
                    """
                    if 'type' in qualifier.key:
                        type = str.strip(qualifier.value, '"')
                        type_line_format = ['\t', '\t', 'Product', type]
                        writer.writerow(type_line_format)
                    """
                    looks for qualifier.key named "product". If found, print line in the format:
                                    Product	<product>
                    """
                    if 'product' in qualifier.key:
                        product = str.strip(qualifier.value, '"')
                        product_line_format = ['\t', '\t', 'Product', product]
                        writer.writerow(product_line_format)
                    """
                    IF OTHER QUALIFIERS ARE NEEDED THEN FORMAT LIKE ABOVE
                    """
















































