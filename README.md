This repository contains code and GUI which simulates multiple ligand bindings on one-dimensional lattice. 
GUI version when runs asks for FASTA sequence (which is stored within files subfolder ) and 
then asks for number of ligands, ligand lengths, co-operativity parameters in the following 
dialogue boxes. It generates probability of binding of each ligands by averaging the 
'nrun' simulations and compares it with exact binding probability as the output.

The files within SRC are supporting files to generate the GUI.


In non-GUI version ligandnongui.m generates probability of binding of 3 ligands (size 4 lattice site each)
as discussed in the article.


The postionmap.m generates map of random postioning of ligands for 10 siulations and stack them vertically.

We have introduced an application section which describes how our random  sequential sampling method can be used to predict optimal netropsin (a non-fluorescent dye) concentration in a particular a competitive binding assay with YOYO-1 used in optical DNA mapping. To find the optimal concentration we use a self information score based on the intensity profile (DNA barcode) of the molecule. We have utilized article 'https://doi.org/10.1371/journal.pone.0121905' eqn 4 for calculating self information score. info_oct6.m calculates the score for each netropsin concentration (500 samples for each concentration and average them). 
The code info_oct6.m uses two functions robustextrema2.m and info_score.m 
