This package contains 3 sub-folders containing of source code for two programs involved in tumor-specific causal 
inference and a set of example files.

	1. TCI_GenerateGlobalDriver: Source code for a program that search a population wide driver for a phenotype using Bayesian approach.  

	2. TCI:  Source code for a program that perform tumor-specific causal inference

	3. Data: A folder contains a set of training data.


Compile program:
	Run Makefile within TCI and TCI_GenerateGlobalDriver to compile programs, which will yield two executable programs:
		1. TCI_GD.  This program searches for a "global driver" for a phenotype at the population level using  
		    Bayesian causal framework 

		2. TCI.  This program perfrom TCI analysis
		
Preparing data:
   	Both TCI and TCI_GD take 3 input matrices:
   		1. An N-by-G matrix, referred to as A matrix, where N is the number of cases, and G is the number of
	   		genes.  Each row represent the somatic genome alteration (SGA) data of a tumor, where a "1" indicates
	   		that the corresponding gene is altered in the current tumor, and "0" otherwise.

	   	2. An N-by-G matrix, referred to as P matrix, where N and G should match those of the A matrix.  Again 
	   		each row represents a tumor, where each element is the prior probability that an SGA, corresponding
	   		to a gene with "1" in A matrix, is a driver of a phenotype (e.g., a differentially expressed gene) in
	   		a tumor.  The sum of the prior probability of a tumor (a row) should be 1.  One can prepare this 
	   		matrix according to available prior knowledge, e.g., the probability that TP53 is a driver gene of a 
	   		tumor and normalize among SGAs in a given tumor, or using a uniform prior if no prior knowledge is 
	   		available. 

   		2. An N-by-D matrix, referred to as E matrix, where N is the number of cases, and D is the number of 
   			phenotypes of interes, e.g., differentially expressed genes (DEGs).  Each row represents phenotypes
   			observed in a tumor, and rows should be matched with those in A & P matrices.  An element in a row 
   			represents whether the phenotype (indexed by D) is present ("1") or not ("0") in a given tumor.  

Performing TCI analysis:
	First, search population-wide driver of each phenotype using TCI_GD
	TCI_GD:
		./TCI_GD -p PmatrixFilePathname -f AmatrixFilePathname -d EmatrixFilePathname -o populationDriverFilePathname


		TCI_GD takes 3 input matrices as described above.  It outputs a D-by-2 common-separated-values (CSV) file,
		in which each each row corresponding to a phenotype and its most probable driver at the population 
		level.  

	Then, perform TCI by utilizing the results from TCI_GD	
	TCI:
		./TCI  -p PmatrixFilePathname -f AmatrixFilePathname -d EmatrixFilePathname -g populationDriverFilePathname
		 -o outputFileDirectoryName [-s startingRow -e endingRow]

		TCI takes 3 input matrices as described above plus the population-wide driver information.  If no addition 
		optional argument (-s and -e) provided, it iterates through each tumor (rows in 3 matrices) as a test case 
		and use the rest of the matrix as training cases to perform a tumor-specific causal inference.  
		For each phenotype that present in a tumor, e.g., a DEG event indicated by a "1" in E matrix, TCI tries to 
		examine the posterior probabilities of each SGA that is present in the same tumor as the candidate cause of 
		the phenotype.  Thus, for each tumor, TCI output a D-by-G matrix in the director provided by "-o" option, in which
		each row contains assigned posterior proability of SGAs in the tumor as the cause of a phenotype 
		(indicated by the row).  

		If one only wants to run a block of cases as test cases, use "-s" and "-e" to indicate the beginning row
		and end row of the block.  TCI will iterate through the cases within the block to perform TCI analysis
		on these cases.


