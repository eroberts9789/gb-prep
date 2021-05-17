import os
import csv
from typing import Tuple

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



def format_record_list(record_list: list) -> list:
    """
    Takes full list of input records and removes any repeated features and any feature keys that don't match the feature keys in the first record. Prints list of incorrectly formatted features if contains extra features after these steps

    :param record_list: list of all records from gbk file
    :return: list of records only containing records with the correct number of features that match desired format
    """
    first_record = record_list[0]
    first_record_feature_keys = [feature.key for feature in first_record["features"]]

    records_with_correct_format = list()
    for record in record_list:
        record["features"] = list(dict.fromkeys(record['features']))

        for feature in record["features"]:
            if feature.key not in first_record_feature_keys:
                record["features"].remove(feature)

        if len(record["features"]) != len(first_record_feature_keys):
            print("Error: record " + record["locus"] + " is formatted incorrectly")
            sys.exit(1)
        else:
            records_with_correct_format.append(record)



def get_record_list():
    """
    Parses through and generates a list of dictionaries containing all required information about records.

    Makes calls to other functions to format location and record data
    """

    with open(INPUT_PATH) as handle:

        record_list = []

        for gbk_record in GenBank.parse(handle):
            record_list.append({
                "locus": gbk_record.locus,
                "features": [feature for feature in gbk_record.features]
            })

    return format_record_list(record_list)



def write_features_file(record_list : list):
    """
    Takes record list as input and writes information to features file

    :param record_list: list containing all properly formatted record data
    """
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



def run():
    """
    Makes calls to functions to collect necessary data from .gbk file and to write output in features file format
    """
    record_list = get_record_list()
    write_features_file(record_list)

