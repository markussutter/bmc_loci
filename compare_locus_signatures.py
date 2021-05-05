#!/usr/bin/python3

import sys

def is_shell_protein(protein_pfam):
    if ("P_" in protein_pfam) or ("H_" in protein_pfam) or ("T_" in protein_pfam) or ("Hp_" in protein_pfam) or ("Ts_" in protein_pfam) or ("Tsp_" in protein_pfam) or ("Tdp" in protein_pfam):
      return True
    else:
      return False

def calculate_score(locus1, locus2):
   match = 0
   shell_protein_multiplier = 1.5
   shell_protein_nearby_range = 3   
   for protein1 in locus1:
      protein1_pfam = protein1[1]
      protein1_position = locus1.index(protein1)      
      for protein2 in locus2:
         protein2_pfam = protein2[1]
         if (protein1_pfam == protein2_pfam) and (protein1_pfam != "#"):
           score = 1 	 
           # 1. increase weight of protein if shell proteins nearby
           if not is_shell_protein(protein1_pfam):    # only add weight from surrounding shell proteins to non-shell proteins
             for test_surrounding in range (protein1_position-shell_protein_nearby_range,protein1_position+shell_protein_nearby_range):
                #previous_protein_id = locus1.index(protein1)-1 #will wrap to end of list if negative
                #previous_protein = locus1[previous_protein_id]  
              if test_surrounding > 0 and test_surrounding < len(locus1)-1:      
                if is_shell_protein(locus1[test_surrounding][1]):
                  score = score*shell_protein_multiplier
                  # one shell protein nearby (1.5 multiplier): 1.5 weight, 2: 2.25, 3: 3.375, 4: 5.0625
           match = match + score                  
   return match

def compare_order(locus1, locus2): 
   locus1_list = []
   locus2_list = []   
   max_match = 0
   for protein1 in locus1:
      locus1_list.append(protein1[1])
   for protein2 in locus2:
      locus2_list.append(protein2[1])
   for protein1 in locus1_list:
     if protein1 != "#":
       if protein1 in locus2_list:
          locus1_match_position = locus1_list.index(protein1)
          locus2_match_position = locus2_list.index(protein1)
         # start search in positive and negative direction each
          match = 0
          counter1 = 1
          counter2 = 1
          try:
             while locus1_list[locus1_match_position+counter1] == locus2_list[locus2_match_position+counter2]:
               match += 1
               if match > max_match:
                 max_match = match
               counter1 += 1
               counter2 += 1	       
          except IndexError:
             pass
          match = 0
          counter1 = 1
          counter2 = -1	  
          try:
             while (locus1_list[locus1_match_position+counter1] == locus2_list[locus2_match_position+counter2]) and (locus2_match_position+counter2>-1):
               match += 1
               if match > max_match:
                 max_match = match
               counter1 += 1
               counter2 += -1	       
          except IndexError:
            pass     

          match = 0
          counter1 = -1
          counter2 = 1	  
          try:
             while (locus1_list[locus1_match_position+counter1] == locus2_list[locus2_match_position+counter2]) and (locus1_match_position+counter1>-1):
               match += 1
               if match > max_match:
                 max_match = match
               counter1 += -1
               counter2 += 1	       
          except IndexError:
            pass     

          match = 0
          counter1 = -1
          counter2 = -1	  
          try:
             while (locus1_list[locus1_match_position+counter1] == locus2_list[locus2_match_position+counter2]) and (locus1_match_position+counter1>-1) and (locus2_match_position+counter2>-1):
               match += 1
               if match > max_match:
                 max_match = match
               counter1 += -1
               counter2 += -1	       
          except IndexError:
            pass     
   return max_match


def main():
 order_match_multiplier = 10
 report_max_scores = 30
 scores_lower_cutoff = 25
 try:
   input_file = sys.argv[1]
 except IndexError:
   print ("usage: ")
   print ("compare_locus_signatures.py locus_signatures_file.txt > output.tsv")
   exit()  
 locus_collection = []
 with open(input_file) as file:
   for locus in file:
      locus = locus.strip() # preprocess line 
      locus_name = locus.split()[0]
      locus_type = locus.split()[1]      
      locus_attributes = locus.split()[2:]
      attributes_list = []
      for attributes in locus_attributes:
         attributes_list.append(attributes.split("-"))
      locus_collection.append([locus_name, locus_type, attributes_list])	 
 # print header for tab-delimited table
 print ('\t'.join(map(str, ('locus1','locus2','total_score','locus1_type','locus2_type'))))

 # all vs all comparison of loci
 for i in range(len(locus_collection)):
  top_scored = []
  for j in range(len(locus_collection)):
    if i != j:
      pfam_id_match = calculate_score(locus_collection[i][2], locus_collection[j][2])
      order_match = compare_order(locus_collection[i][2], locus_collection[j][2])
      total_score = pfam_id_match + order_match_multiplier*order_match
      if total_score > scores_lower_cutoff:
          top_scored.append([total_score, j, pfam_id_match, order_match]) 
  top_scored_sort = sorted(top_scored, reverse=True)  
  for topx in range (0,report_max_scores):
   try:
    total_score = top_scored_sort[topx][0]   
    top_id = top_scored_sort[topx][1]
    pfam_id_match = top_scored_sort[topx][2]      
    order_match = top_scored_sort[topx][3]      
    print ('\t'.join(map(str, (locus_collection[i][0], locus_collection[top_id][0], total_score, locus_collection[i][1], locus_collection[top_id][1]))))
   except IndexError:
    pass

main()

