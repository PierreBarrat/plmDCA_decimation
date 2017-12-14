# plmDCA_decimation
Research script -- 

--------------------------------------
pseudolikelihood maximization DCA (modified) code used here was taken from 

 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013)
  
  and

	M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
	maximization for direct-coupling analysis of protein structure
	from many homologous amino-acid sequences, J. Comput. Phys. 276, 341-356 (2014)
  
--------------------------------------

$./main.sh to start

main.sh simply starts matlab and launches main\_script.m

main\_script.m does three things
- Reads a line from a file called list.txt. This should be the absolute path to a fasta file. An out directory is created based on this name.
- Preprocesses this data for my decimation script to use. Translates the fasta file to a number mapping, and computes weights for every sequences. Preprocessed data are stored in the out directory.
- the decimation.m function is started in this out directory. It performs the decimation and stores results in a sub directory. 

Parameters of the decimation can be set inside the "decimation.m" file, in the "sources" directory (not practical, but it could be changed). These include
- Number of decimation rounds
- Fraction of couplings set to 0 at each round (I recommend 0.01 since this what I've worked with so far).
- Regularization strength
- Possibility to have plmDCA use multiple cores. I am not sure it works in this case, because matlab is started from command line with no java virtual machine. I can check this if necessary

For now, I have tried with a list.txt containing two alignments. One of them is not particularly small, so it's not ideal to test things. For testing purposes, I recommend choosing a low number of decimation rounds in "decimation.m"


