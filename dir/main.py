from Bio import GenBank
"""
GOAL:   SCAN THROUGH MULTIPLE SEQUENCE FILES IN GENBANK FILE, EXTRACT AND PRINT ALL CDS RANGES
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
    
"""

#open genbank file, look at each record
def get_features_list(filename):
    """
    CURRENTLY ONLY HAVE "ACLSV.gbk" AS WORKING FILENAME
    Iterate through all records in genbank file, add all Features in Feature block to Features list

    Feature format:
    Feature(key='CDS', location='150..5801')
    """
    with open(filename) as handle:
        Features = []
        for record in GenBank.parse(handle):
            for feature in record.features:
                Features.append(feature)
                print(feature)

    return Features




def format_locations(Features):
    """
    Takes list of CD features and outputs a nicely formatted list of location tuples (loc0, loc1)
    We first make a list of locations from Features list
    """
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

    return location_tuples


Features = get_features_list("ACLSV.gbk")
print(format_locations(Features))



















