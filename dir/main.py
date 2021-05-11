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
with open("ACLSV.gbk") as handle:

    #visit all feature blocks fro each record in genbank file and add each CDS feature to CDS list
    Feature_objects = []
    for record in GenBank.parse(handle):
        for feature in record.features:
            Feature_objects.append(feature)

#extract locations from Feature objects and add them to locations list
locations = []
for cds in Feature_objects:
    locations.append(cds.location)

print(locations)











