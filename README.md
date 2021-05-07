# bmc_loci
A Catalog of the Diversity and Ubiquity of Bacterial Microcompartments

Scripts and example data

uniprot_to_fasta_annotated.py
convert Uniprot text file to fasta format, extracting gene id, pfams, taxid and accession and stores it in the fasta header
example:
./uniprot_to_fa_annotated.py example_data/Hoch_5814.txt > Hoch_5814.fa

locus_fasta_to_signature_and_diagram.py
runs hmms against locus in fasta format (with extra information in fasta header), creates a short signature and visualizes the locus
example:
./locus_fasta_to_signature_and_diagram.py example_data/ACI/A3B65__locus_1_main_sat.fa ACI > out.sig

compare_locus_signatures.py
uses the locus signatures from all BMC loci (or a subset) to compare them against each other, output is tab separated for import into Cytoscape
example:
./compare_locus_signatures.py example_data/SPU_loci_signatures.txt > out.tsv


Dependencies

biopython 1.78:
pip3 install biopython
(https://biopython.readthedocs.io/en/latest/index.html)

HMMER 3.1b2 / hmmer tools:
http://hmmer.org/
hmmscan and hmmsearch should be in path when running scripts that use HMM analysis 

Python zelle graphics library v5:
pip3 install --user http://bit.ly/csc161graphics
(https://www.pas.rochester.edu/~rsarkis/csc161/python/pip-graphics.html)

Cytoscape 3.7.2 for analysis of locus clustering:
https://cytoscape.org/
install clustermaker2 plugin:
-> Apps -> App Manager -> search for clustermaker2, then install
