
from Bio import GenBank
import csv


global input_name
input_name = "12G422__GRGV.gbk"
global output_name
output_name = "features_output.txt"

def format_locs(location):
    """
    Takes string describing loaction and returns nicely formatted tuple
    :param location: location string in the format 123..1234
    :return: location tuple (123, 1234)
    """
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

def format_record_list(first_record : list, record_list : list):
    """
    Takes full list of input records and removes any repeated features and any feature keys that don't match the feature keys in the first record. Prints list of incorrectly formatted features if contains extra features after these steps

    :param first_record: record formatted correctly, compare others to this format
    :param record_list: list of all records from gbk file
    :return: list of records only containing records with the correct number of features that match desired format
    """

    first_record_feature_keys = []
    for feature in first_record['features']:
        first_record_feature_keys.append(feature.key)

    must_fix_format = []
    result = []
    for record in record_list:
        #REMOVE ANY REPEATED FEATURES BY TURNING INTO DICT AND THEN BACK TO LIST
        record['features'] = list(dict.fromkeys(record['features']))

        for feature in record['features']:
            #IF THERE IS A FEATURE KEY THAT DOESN'T MATCH THE FEATURE KEYS IN THE FIRST FEATURE, REMOVE THAT FEATURE FROM RECORD
            if feature.key not in first_record_feature_keys:
                record['features'].remove(feature)

        #CREATE LIST OF POORLY FORMATTED RECORDS THAT NEED TO BE FIXED
        if len(record['features']) > len(first_record_feature_keys):
            must_fix_format.append(record)
        #CREATE LIST OF CORRECTLY FORMATTED RECORDS
        if record not in must_fix_format:
            result.append(record)

    #IF LIST CONTAINS ANY POORLY FORMATTED RECORDS, PRINT THEM OUT
    if len(must_fix_format) > 0:
        print('SOME OF YOUR RECORDS ARE FORMATTED INCORRECTLY. THE FOLLOWING ARE NOT INCLUDED IN YOUR OUTPUT FILES.')
        for record in must_fix_format:
            print(record['locus'])

    return (result)




#PARSE GENBANK FILE AND STORE ALL REQUIRED INFO IN RECORD DICTIONARIES IN RECORD_LIST LIST
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

first_record = record_list[0]
record_list = format_record_list(first_record, record_list)


#NOW WRITE EVERYTHING TO FEAUTURES FILE
with open(output_name, 'w', newline='') as output_file:
    writer = csv.writer(output_file, delimiter='\t', escapechar=' ', quoting=csv.QUOTE_NONE)

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

            #ITERATES THROUGH AND PRINTS QUALIFIERS IN FIRST FILE SO ALL OUTPUT MATCHES
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

            writer.writerow('\t')

















































