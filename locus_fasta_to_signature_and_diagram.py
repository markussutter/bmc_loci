#!/usr/bin/python3

import sys
import os
from Bio import SeqIO
from os.path import normpath, basename
from graphics import *

def hmm_check_for_ep (filename, ep_ref_hmm):
  # run hmmscan and extract information from output files for encapsulation peptide HMM analysis
  ep_top_hit_dict = {}
  ep_score_cutoff = 5E-5
  os.system("hmmscan --tblout _ep_tbl.tmp "+ep_ref_hmm+" "+filename+" > _ep_log")
  with open('_ep_tbl.tmp') as file:
    for line in file:
      if not "#" in line[0]:
        target_id = line.split()[2]
        gene_id = target_id.split("|")[1]
        target_id = gene_id
        hit_id = line.split()[0]
        hit_score = float(line.split()[7])       
        if target_id in ep_top_hit_dict:
           current_top_hit = ep_top_hit_dict[target_id][1]
           if current_top_hit > hit_score:
              ep_top_hit_dict[target_id] = (hit_id, hit_score)
        else:
           if hit_score < ep_score_cutoff:
             ep_top_hit_dict[target_id] = (hit_id, hit_score) 
  with open('_ep_log') as file:
    current_id_top_score = "##nohit##"
    current_target_id = "##none##"
    location_information_below = False   
    done_with_query = False
    for line in file:
      if "Query:" in line:
        current_id = line.split()[1]
        gene_id = current_id.split("|")[1]
        try:
           current_id_top_score = str(ep_top_hit_dict[gene_id][1])
        except KeyError:
           current_id_top_score = "##nohit##" 
        try:
           current_target_id = str(ep_top_hit_dict[gene_id][0])
           done_with_query = False
           location_information_below = False           
        except KeyError:
           current_target_id = "##none##" 
      if  "envfrom" in line:
         location_information_below = True
      if not(done_with_query) and location_information_below and (current_id_top_score in line):
         ep_match_from = int(line[60:67])
         ep_match_to = int(line[69:74])         
         ep_top_hit_dict[gene_id] = (current_target_id, current_id_top_score, ep_match_from, ep_match_to)
         done_with_query = True
  return (ep_top_hit_dict)

def draw_operon (type_assignment, gene_list, organism, locusname, total_protein_length):
 # draw operon on screen using graphics library
 type_list = []
 # protein type vs color definitions
 type_list.append([("-", ),("grey")])  # if pfam field empty
 type_list.append([("alcdh",),("green")])  # alc_dh
 type_list.append([("HMMalddh",),("red")])   # alddh
 type_list.append([("HMMpropdeh", "HMMetly", "HMMgre","HMMaldol","HMMspu"), ("purple")])  #signature enzymes
 type_list.append([("HMMreg",), "orange"])
 type_list.append([("HMMptac","HMM01515PTA_B"),("magenta")])  # ptac
 type_list.append([("H_","Hp_"),("blue")])
 type_list.append([("T_","Tsp_","Ts_"),("cyan")])
 type_list.append([("Tdp_",),(color_rgb(0,160,160))])
 type_list.append([("P_", ),("yellow")])
 
# read in operon
 if total_protein_length < 3000:
     total_protein_length = 3000   # for scaling
 width = 1400
 height = 700
 win = GraphWin(organism, width+250, height, autoflush=False) # give title and dimensions
 win.setBackground('white')
 win.update() 
 current_pos = 0
 offset = 52
 previous_gene_number = int(gene_list[0][0])
 for field in gene_list:
    id_number = field[0]
    gene_length = int(width*(field[1]/total_protein_length))   # scaled to overall width
    assignment = field[2]
    hmm_bmc_type = field[3]
    direction = field[4]
    ep_string = field[5]
    id_string = id_number+" "+assignment+" "+hmm_bmc_type+" "+ep_string
    id_string_clean = id_string.replace("()","").replace("(-)","")
    difference = int(id_number) - previous_gene_number
    if abs(difference) > 50:             # if large break in gene id numbers, insert spacer
        offset = offset+10
        diff_label = Text(Point(current_pos, offset), "(.."+str(difference)+"..)")
        diff_label.draw(win)
        offset = offset+26
    # drawing of bar/arrow for single gene
    point1 = Point(current_pos,offset)
    point2 = Point(current_pos+gene_length,offset)  
    rect = Line(point1, point2)
    rect.setWidth(15)
    if direction == "-":          # directionality of arrows
      rect.setArrow('first')
    if direction == "+":
      rect.setArrow('last')
    if (direction != "-") and (direction != "+"):
      rect.setArrow('none')
    bar_color = 'black'
    for pftype in type_list:
      for pftypemember in pftype[0]:
         if pftypemember in assignment:
            bar_color = pftype[1]
    rect.setFill(bar_color)
    rect.draw(win)
    # drawing circle for encapsulation peptide position
    center_pos = Point(current_pos+int(gene_length/2), offset)
    left_pos = Point(current_pos, offset)
    right_pos = Point(current_pos+int(gene_length), offset)
    if ep_string != "":
       ep_pos = center_pos
       if "NTERM" in ep_string:
          if (direction == "-"):
            ep_pos = right_pos
          else:
            ep_pos = left_pos
       if "CTERM" in ep_string:
          if (direction == "-"):
             ep_pos = left_pos
          else:
             ep_pos = right_pos
       ep_circle = Circle(ep_pos, 6)
       ep_circle.setFill('brown')
       ep_circle.setOutline('yellow')
       ep_circle.draw(win) 
    label = Text(Point(current_pos+gene_length+len(id_string_clean)*5, offset), id_string_clean)
    # coloring text green if type assignment matches best type from HMM match
    if type_assignment in hmm_bmc_type:
       label.setFill('green')
    else:
       if "-" in field[3]:
         label.setFill('grey')
       else:
         label.setFill('red') 
    label.draw(win)
    # calculate difference between gene numbers and show that on right side 
    previous_gene_number = int(id_number)
    current_pos = current_pos + gene_length
    offset = offset+16
    if offset > height - 60:   # new column if at bottom of screen
        offset = 100
 # drawing legend
 legend_start = 250
 ls = legend_start
 message1 = Text(Point(89, ls), "BMC-P")
 message1.draw(win) 
 rect = Line(Point(40,ls-10), Point(40,ls+10))
 rect.setWidth(18)
 rect.setFill('yellow')
 rect.draw(win)

 ls = legend_start+20
 message1 = Text(Point(89, ls), "BMC-H")
 message1.draw(win) 
 rect = Line(Point(40,ls-10), Point(40,ls+10))
 rect.setWidth(18) 
 rect.setFill('blue')
 rect.draw(win)
 
 ls = legend_start+40
 message1 = Text(Point(89, ls), "BMC-T")
 message1.draw(win) 
 rect = Line(Point(40,ls-10), Point(40,ls+10))
 rect.setWidth(18)
 rect.setFill('cyan')
 rect.draw(win)

 ls = legend_start+60
 message1 = Text(Point(155, ls), "Aldehyde dehydrogenase")
 message1.draw(win) 
 rect = Line(Point(40,ls-10), Point(40,ls+10))
 rect.setWidth(18)
 rect.setFill('red')
 rect.draw(win)

 ls = legend_start+80
 message1 = Text(Point(87, ls), "PTAC")
 message1.draw(win) 
 rect = Line(Point(40,ls-10), Point(40,ls+10))
 rect.setWidth(18)
 rect.setFill('magenta')
 rect.draw(win)

 ls = legend_start+100
 message1 = Text(Point(129, ls), "Signature enzyme")
 message1.draw(win) 
 rect = Line(Point(40,ls-10), Point(40,ls+10))
 rect.setWidth(18)
 rect.setFill('purple')
 rect.draw(win)

 ls = legend_start+120
 message1 = Text(Point(149, ls), "Alcohol dehydrogenase")
 message1.draw(win) 
 rect = Line(Point(40,ls-10), Point(40,ls+10))
 rect.setWidth(18)
 rect.setFill('green')
 rect.draw(win)

 ls = legend_start+140
 message1 = Text(Point(100, ls), "Regulator")
 message1.draw(win) 
 rect = Line(Point(40,ls-10), Point(40,ls+10))
 rect.setWidth(18)
 rect.setFill('orange')
 rect.draw(win)
 
 message1 = Text(Point((width+250)/2, 20), organism)
 message2 = Text(Point((width+250)/2, 60), type_assignment)
 message3 = Text(Point((width+250)/2, 40), locusname)
 message1.draw(win)
 message2.setFill('green')
 message2.draw(win) 
 message3.draw(win)
# win.postscript(file=locusname+".eps")  output to eps file
 win.getKey()
 win.close()


def fa_header_to_fields(fasta_header):
   # read information from fasta file
   genome_id = fasta_header.split("|")[0]
   gene_id = fasta_header.split("|")[1]
   directionality = fasta_header.split("|")[2]
   try:
      cog_field = fasta_header.split("|")[3]
   except IndexError:
      cog_field = ""
   try:
      pfam_field = fasta_header.split("|")[4]
   except IndexError:
      pfam_field = "-"
   if "pfam" in cog_field:
      pfam_field = cog_field
      cog_field = ""
   pfam_list = pfam_field.replace(",","")
   return(genome_id, gene_id, directionality, cog_field, pfam_list)

def hmm_search_hit_list (filename, ref_hmm):
  # run hmmsearch with combined hmm library against full locus fasta file, store best hits in dictionary
  os.system("hmmsearch -T 40 --tblout _tbl.tmp "+ref_hmm+" "+filename+" > /dev/null")
  top_hit_dict = {}
  with open('_tbl.tmp') as file:
    for line in file:
      if not "#" in line[0]:
        target_id = line.split()[0]
        gene_id = target_id.split("|")[1]
        target_id = gene_id
        hit_id = line.split()[2]
        hit_score = float(line.split()[5])
        if target_id in top_hit_dict:
           current_top_hit = top_hit_dict[target_id][0]
           if current_top_hit < hit_score:
              top_hit_dict[target_id] = (hit_score, hit_id)
        else:
           top_hit_dict[target_id] = (hit_score, hit_id)
  return (top_hit_dict)  

def main():
 try:
   target_fa = sys.argv[1]
 except IndexError:
   print ("usage:")
   print ("locus_fasta_to_signature_and_diagram.py locus_fasta.fa [BMC type] > locus_signature.sig")
   exit()
 locus_name = str(basename(normpath(target_fa))).replace(".fa","")

 try:
   type_assignment = sys.argv[2]
 except IndexError:
   type_assignment = "unassigned"

# reference files
 ref_hmm = "reference_data/combined_hmms.hmm"
 ep_ref_hmm = "reference_data/ep_hmms.hmm"
 org_table = "reference_data/genome_vs_org_table.txt"
 org_dict = {}
 with open(org_table) as file:
    for org in file:
       line = org.strip()
       name = line.split()[0]
       org = " ".join(line.split()[1:])
       org_dict[name] = org

# initialize variables
 gene_list = []     # collection information for locus visualization
 sig_string = ""    # one line space separated locus information summary with all genes and HMM matches
 total_protein_length = 0

# run hmms against locus fasta file, store hits in dictionaries
 hmm_hits = hmm_search_hit_list (target_fa, ref_hmm)
 ep_hit_dict = hmm_check_for_ep (target_fa, ep_ref_hmm)

 with open(target_fa, "rU") as handle:
  for record in SeqIO.parse(handle, "fasta"):
    protein_length = len(str(record.seq))
    total_protein_length = total_protein_length + protein_length  
    (genome_id, gene_id, directionality, cog_field, pfam_list) = fa_header_to_fields(record.id)
    if (pfam_list == "pf00936") or (pfam_list == "pf03319"):   # don't show pfams for shell proteins
       pfam_list = ""
    if gene_id in hmm_hits:           # extract information from HMM hits
       protein_info = hmm_hits[gene_id][1].split(".")[0]+"("+pfam_list+")"
       bmctype = hmm_hits[gene_id][1].split(".")[1]             # extract bmctype from HMM name string
       sig_string = sig_string + str(gene_id)+"-"+hmm_hits[gene_id][1].split(".")[0]+"-"+bmctype+" "
    else:
       protein_info = "-"+pfam_list
       bmctype = "-"
       if pfam_list == "":
         pfam_list = "#"
       sig_string = sig_string + str(gene_id)+"-"+pfam_list+"-# "
    if gene_id in ep_hit_dict.keys():
       ep_match_from = ep_hit_dict[gene_id][2]
       ep_match_to = ep_hit_dict[gene_id][3]
       if ep_match_from < 25:
          ep_match_string = "EP_NTERM"
       elif (protein_length-ep_match_to) < 25:
          ep_match_string = "EP_CTERM"
       else:
          ep_match_string = ""                        
       ep_str = ep_match_string
    else:
       ep_str = ""
    gene_list.append((gene_id, protein_length, protein_info, bmctype, directionality, ep_str))
  try:
    organism = "".join(org_dict[genome_id])
  except KeyError:
    organism = "unknown"       

 print (locus_name, type_assignment, sig_string)
 draw_operon (type_assignment, gene_list, organism, locus_name, total_protein_length)

main()
