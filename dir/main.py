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
    """
    Iterate through all records in genbank file, add all Features in Feature block to Features list
    
    Feature format:
    Feature(key='CDS', location='150..5801')
    """
    Features = []
    for record in GenBank.parse(handle):
        for feature in record.features:
            Features.append(feature)

"""
We're interested in location field of Feature object, create list of locations from Features list
"""
locations = []
for cds in Features:
    locations.append(cds.location)

print(locations)

#output_file = open("ACLSVfeatures", 'w')
#output_file.close()













