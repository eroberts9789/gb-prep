import os
import csv

from Bio import GenBank

INPUT_PATH = "GRVFV.gbk"
OUTPUT_PATH = "features_output.txt"


def format_locations(location):
    """
    Takes string describing location and returns nicely formatted tuple

    :param location: location string in the format 123..1234
    :return: location tuple (123, 1234)
    """
    return location.split("..")



def format_record_list(first_record : list, record_list : list):
    """
    Takes full list of input records and removes any repeated features and any feature keys that don't match the feature keys in the first record. Prints list of incorrectly formatted features if contains extra features after these steps

    :param first_record: record formatted correctly, compare others to this format
    :param record_list: list of all records from gbk file
    :return: list of records only containing records with the correct number of features that match desired format
    """

    first_record_feature_keys = []

    for feature in first_record["features"]:
        first_record_feature_keys.append(feature.key)

    must_fix_format = []
    result = []
    for record in record_list:
        record["features"] = list(dict.fromkeys(record['features']))

        for feature in record["features"]:
            if feature.key not in first_record_feature_keys:
                record["features"].remove(feature)

        #CREATE LIST OF POORLY FORMATTED RECORDS THAT NEED TO BE FIXED
        if len(record["features"]) != len(first_record_feature_keys):
            must_fix_format.append(record)
        #CREATE LIST OF CORRECTLY FORMATTED RECORDS
        if record not in must_fix_format:
            result.append(record)

    #IF LIST CONTAINS ANY POORLY FORMATTED RECORDS, PRINT THEM OUT
    if len(must_fix_format) > 0:
        print("SOME OF YOUR RECORDS ARE FORMATTED INCORRECTLY. THE FOLLOWING ARE NOT INCLUDED IN YOUR OUTPUT FILES.")
        for record in must_fix_format:
            print(record["locus"])

    return (result)



def run():
    with open(INPUT_PATH) as handle:

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


    with open(OUTPUT_PATH, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", escapechar=" ", quoting=csv.QUOTE_NONE)

        for record in record_list:
            writer.writerow([">Feature", record["locus"]])

            for feature, first_feature in zip(record["features"], first_record["features"]):

                positions = format_locations(feature.location)
                writer.writerow([positions[0], positions[1], str(feature.key)])

                for qualifier in first_feature.qualifiers:
                    if 'cds_type' in qualifier.key:
                        writer.writerow(['\t', '\t', '\t', 'product', str.strip(qualifier.value, '"')])

                    if 'product' in qualifier.key:
                        writer.writerow(['\t', '\t', '\t', 'product', str.strip(qualifier.value, '"')])

                writer.writerow('\t')


run()


