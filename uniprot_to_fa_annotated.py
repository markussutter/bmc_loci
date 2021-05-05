#!/usr/bin/python3
import sys
import re
import Bio.SwissProt as sp 

def uniprot_to_geneid (name):     # some issues with uniprot files that contain entries from multiple genomes
       if "ORFNames=" in name:
          start_position = name.find("ORFNames=")+9
          gene_id = name[start_position:].split(' ',1)[0]
          return gene_id
       if "OrderedLocusNames=" in name:
          start_position = name.find("OrderedLocusNames=")+18
          gene_id = name[start_position:].split(' ',1)[0]
          return gene_id
       return "nogeneid"  

def geneid_to_genome_no (gene_name):
  for n in re.findall(r'\d+',gene_name):            # find last number in string
    lastn = n
  try:
     gene_no = lastn
     genome_id = gene_name[0:len(gene_name)-len(lastn)]
     return genome_id, lastn
  except NameError:
     return ""

def get_pfam(uniprot_cr_tag):
   pfam_list = ""
   for cr in uniprot_cr_tag:
        if cr[0] == "Pfam": 
            pfam_list = pfam_list+(cr[1].lower())     # extract pfams, concatenate 
   if pfam_list == "":
     pfam_list = "#"
   return pfam_list         

def main():
 try:
    filename = sys.argv[1]
 except IndexError:
    print ("usage:")
    print ("uniprot_to_fa_annotated.py uniprot_text_format_file.txt")
    exit()

 with open(filename) as handle: 
     records = sp.parse(handle) 
     for record in records: 
         gene_name = uniprot_to_geneid(record.gene_name)
         (genome_id, gene_id) = geneid_to_genome_no(gene_name)
         pfam = str(get_pfam(record.cross_references))
         taxid = str(record.taxonomy_id[0]).replace(" ","_")
         print(">"+genome_id+"|"+gene_id+"|#|#|"+pfam+"|"+taxid+"|"+str(record.accessions[0]))
         print(record.sequence)
         
main()
