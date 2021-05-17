import os
import csv

from Bio import GenBank

INPUT_PATH = "GRVFV.gbk"
OUTPUT_PATH = "features_output.txt"


def format_locations(location : str):
    """
    Takes string describing location and returns nicely formatted tuple

    :param location: location string in the format 123..1234
    :return: location tuple (123, 1234)
    """
    return location.split("..")



def format_record_list(record_list : list):
    """
    Takes full list of input records and removes any repeated features and any feature keys that don't match the feature keys in the first record. Prints list of incorrectly formatted features if contains extra features after these steps

    :param record_list: list of all records from gbk file
    :return: list of records only containing records with the correct number of features that match desired format
    """
    first_record = record_list[0]
    first_record_feature_keys = [feature.key for feature in first_record["features"]]

    for record in record_list:

        record["features"] = list(dict.fromkeys(record['features']))
        for feature in record["features"]:
            if feature.key not in first_record_feature_keys:
                record["features"].remove(feature)

        records_to_fix = [record for record in record_list if len(record["features"]) != len(first_record_feature_keys)]

        records_with_correct_format = [record for record in record_list if record not in records_to_fix]


    if len(records_to_fix) > 0:
        print("SOME OF YOUR RECORDS ARE FORMATTED INCORRECTLY. THE FOLLOWING ARE NOT INCLUDED IN YOUR OUTPUT FILES.")
        for record in records_to_fix:
            print(record["locus"])

    return (records_with_correct_format)



def run():
    """
    Parses through and generates a list of dictionaries containing all required information about records. Makes calls to functions to format records and locations and writes necessary information to output file.
    """

    with open(INPUT_PATH) as handle:

        record_list = []
        for gbk_record in GenBank.parse(handle):

            record = {
                "locus": "",
                "features": []
            }

            record["locus"] = gbk_record.locus
            record["features"] = [feature for feature in gbk_record.features]
            record_list.append(record)

    record_list = format_record_list(record_list)

    first_record = record_list[0]
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


