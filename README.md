This repository contains code and GUI which simulates multiple ligand bindings on one-dimensional lattice. 
GUI version when runs asks for FASTA sequence (which is stored within files subfolder ) and 
then asks for number of ligands, ligand lengths, co-operativity parameters in the following 
dialogue boxes. It generates probability of binding of each ligands by averaging the 
'nrun' simulations and compares it with exact binding probability as the output.


In non-GUI version ligandnongui.m generates probability of binding of 3 ligands (size 4 lattice site each)
as discussed in the article.


The postionmap.m generates map of random postioning of ligands for 10 siulations and stack them vertically.
