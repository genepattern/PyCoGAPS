# PyCoGAPS (v1)

**Description**: 
See [https://github.com/FertigLab/pycogaps](https://github.com/FertigLab/pycogaps) for more details of PyCoGAPS.

Coordinated Gene Activity in Pattern Sets (CoGAPS) implements a Bayesian MCMC matrix factorization algorithm, GAPS, and links it to gene set statistic methods to infer biological process activity. It can be used to perform sparse matrix factorization on any data, and when this data represents biomolecules, to do gene set analysis.

Gene Association in Pattern Sets (GAPS) infers underlying patterns in a matrix of measurements that can
be interpreted as arising from the multiplication of two lower dimensional matrices. The first development of
this code in R/Biocondcutor was focused on gene expression analysis, however the original concept was used
in spectral imaging. The approach is a general form of matrix factorization using a stochastic algorithm.
 While in this doc we will focus on gene expression analysis for concreteness, but the factorization is applicable more broadly.

The Markov chain Monte Carlo (MCMC) matrix factorization that infers patterns also infers the extent
to which individual genes belong to these patterns. The CoGAPS algorithm extends GAPS to infer the
coordinated activity in sets of genes for each of the inferred patterns based upon and to refine gene set
membership based upon.


**Authors**: Fertig Lab, Johns Hopkins University; wrapped as a genePattern module by Ted Liefeld - Mesirov Lab, UCSD

**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help)


## Summary

GAPS seeks a pattern matrix (P) and the corresponding distribution matrix of weights (A) whose product forms a mock data matrix (M) that represents the expression data D within noise limits (ε). That is,

D = M + ε = AP + ε.

The number of rows in P (columns in A) defines the number of biological patterns that GAPS will infer from the measured microarray data or equivalently the number of nonorthogonal basis vectors required to span the data space. As in the Bayesian Decomposition algorithm, the matrices A and P in GAPS are assumed to have the atomic prior. In the GAPS implementation, αA and αP are corresponding parameters for the expected number of atoms which map to each matrix element in A and P, respectively. The corresponding matrices A and P are found by MCMC sampling.

 

CoGAPS infers coordinated activity in gene sets active in each row of the pattern matrix P found by GAPS in a single step, by running both GAPS and then performing the statistical analysis of calcCoGAPSStat. Specifically, CoGAPS computes a Z-score based statistic on each column of the A matrix. The resulting Z-score for pattern p and gene set i, Gi , with G elements is given by Zi,p = 1 / G SUM g∈Gi (Agp/ Asdgp ) where g indexes the genes in the set and Asdgp is the standard deviation of Agp obtained from the MCMC sampling in GAPS. CoGAPS then uses random sample tests to convert the Z-scores from eq. (3.2) to p values for each gene set.

## References

If you use the CoGAPS package for your analysis please cite:  EJ Fertig, J Ding, AV Favorov, G
Parmigiani, and MF Ochs (2010) CoGAPS: an R/C++ package to identify patterns and biological process
activity in transcriptomic data. Bioinformatics 26: 2792-2793.

To cite the CoGAPS algorithm use:  MF Ochs (2003) Bayesian Decomposition in The Analysis of
Gene Expression Data: Methods and Software G Parmigiani, E Garrett, R Irizarry, and S Zeger, ed. New
York: Springer Verlag.


To cite the gene set statistic use:  MF Ochs, L Rink, C Tarn, S Mburu, T Taguchi, B Eisenberg, and
AK Godwin (2009) Detection of treatment-induced changes in signaling pathways in gastrointestinal stromal
tumors using transcriptomic data. Cancer Research 69: 9125-9132.


To site the set-membership refinement statistic use:  EJ Fertig, AV Favorov, and MF Ochs (2012)
Identifying context-specific transcription factor targets from prior knowledge and gene expression data. 2012
IEEE International Conference on Bioinformatics and Biomedicine, B310, in press.
Please contact Elana J. Fertig ejfertig@jhmi.edu or Michael F. Ochs ochsm@tcnj.edu for assistance.


## Source Links
* [The Fertig Lab PyCoGAPS source repository](https://github.com/FertigLab/pycogaps)
* [Genepattern PyCoGAPS Module source repository](https://github.com/genepattern/PyCoGAPS/)

## Parameters

### Standard Parameters
| Name | Description  | Default Value |
|---------|--------------|----------------|
| input file * | 	Input data file in csv or gct format. |  |
| output filename * | The result output file name (output is saved as a .pkl file).  |  |
| num patterns * | The number of patterns PyCoGAPS will learn. | 3 |
| num iterations * | The number of iterations for each phase of the algorithm. | 1000 |
| seed * | Random number generator seed. | 0 |
| use sparse optimization * | When true, speeds up performance with sparse data (roughly >80% of data is zero), note this can only be used with the default uncertainty. | False |
| transpose data * | 	Transpose the dataset before processing. | False |

### Run Parameters
| Name | Description  | Default Value |
|---------|--------------|----------------|
| num threads * | 		The maximum number of threads to run on. |  1 |
| messages * | 		When True, display additional outputs to stdout.txt. | False |
| output frequency * | 		The number of iterations between each output (set to 0 to disable status updates). | 500 |
| uncertainty * | 	The uncertainty matrix - in csv format. |
| checkpoint out filename * | 	The name of the checkpoint file to create. |
| checkpoint in file | 		If this is provided, CoGAPS runs from the checkpoint contained in this file. |
| worker ID * | 		If calling CoGAPS in parallel the worker ID can be specified. |
| asynchronous updates * | Enable asynchronous updating which allows for multi-threaded runs. | False |
| n snapshots * | 	Sets how many snapshots to take in each phase, setting this to 0 disables snapshots. | 0 |
| snapshot phase * | 		During which phase to take snapsjots in e.g. "equilibration", "sampling", "all". | sampling |


### Sparsity Parameters
| Name | Description  | Default Value |
|---------|--------------|----------------|
| alpha A * | 		The sparsity parameter for the feature matrix. |  .01 |
| alpha P * | 		The sparsity parameter for the sample matrix. | .01 |
| max Gibbs Mass A * | 	The atomic mass restriction for the feature matrix. | 100 |
| max Gibbs Mass P * | 	The atomic mass restriction for the sample matrix. | 100 |

### Distributed Parameters
| Name | Description  | Default Value |
|---------|--------------|----------------|
| distributed | 		Either null (None) or genome-wide. | None  |
| num sets * | 		The number of sets to break data into. | 4 |
| cut * | 	The number of branches at which to cut dendrogram used in pattern matching.  | <num patterns> |
| min NS * | 	The minimum of individual set contributions a cluster must contain. | math.ceil(cut / 2) |
| explicit sets | 		Whether to specify subsets by index or name. | None  |
| sampling annotation | 		Specify categories along the rows (cols) to use for weighted sampling. |   |
| sampling weight | 		The weights associated with sampling annotation. |   |

### Additional Parameters
| Name | Description  | Default Value |
|---------|--------------|----------------|
| subset indices * | 		The set of indices to use from the data. |   |
| subset dimension * | 		Which dimension (rows, columns) to subset. |  |
| gene names * | 	The vector of names of genes in the data. |  |
| sample names * | 		The vector of names of samples in the data. |  |
| fixed patterns * | 		The provided fixedPatterns matrix (either ‘A’ or ‘P’ matrix) allows for manual pattern matching. In distributed CoGAPS (distributed=’genome-wide’), the first phase is skipped and fixedPatterns is used for all sets. In standard CoGAPS, fixedPatterns allows for fixed runs. |   |
| which matrix fixed * | 			Either 'A' or 'P', indicating which matrix is fixed. |  A |
| take pump samples * | 		Whether or not to take PUMP samples. |    |
| hdf key * | 		Hdf key for reading .h5 files. |    |
| hdf row key * | 		Hdf row key for reading .h5 files. |   |
| hdf col key * | 			Hdf column key for reading .h5 files. |    |
 
\*  required

## Input Files

1. input file  
    A file containing the matrix to test. Formatted as gct or csv.
    
2. uncertainty 
    A file containing an uncertainty matrix. Formatted as csv.

3. checkpoint in file  
    A file a previously generated checkpoint to start from. Formatted as a pkl file.

4. explicit sets  
    Whether to specify subsets by index or name.  csv format that will convert to a python list.

5. sampling annotation 
    Specify categories along the rows (cols) to use for weighted sampling. csv format that will convert to a python list.
 
6. sampling weight  
   The weights associated with sampling annotation.. txt format that will converted to a python dictionary.
 
7. subset indices   
   	The set of indices to use from the data.  csv format that will convert to a python list.

8. gene names  
    	The vector of names of genes in the data. csv format that will convert to a python list.

9. sample names  
    The vector of names of samples in the data. csv format that will convert to a python list.

10. fixed patterns
  
   	The provided fixedPatterns matrix (either ‘A’ or ‘P’ matrix). Formatted as csv.


 
## Output Files
<!-- list and describe any files output by the module -->

1. \<output_filename\>.pkl  
    Pickle (.pkl) formatted dictionary of the result as two representations is stored: an anndata object. CoGAPS stores the lower dimensional representation of the samples (P matrix) in the .var slot and the weight of the features (A matrix) in the .obs slot. The standard deviation across sample points for each matrix are stored in the .uns slots.
2. stdout.txt
    This is standard output from the Python script. Sometimes helpful for debugging.

## Example Data
<!-- provide links to example data so that users can see what input & output should look like and so that they and we can use it to test -->

Input file:  
[GIST.csv](https://github.com/FertigLab/pycogaps/blob/master/data/GIST.csv)


## License

`PyCoGAPS` is distributed under a modified BSD license available at [https://github.com/genepattern/PyCoGAPS/blob/v2/LICENSE.](https://github.com/genepattern/PyCoGAPS/blob/v2/LICENSE)

## Version Comments

| Version | Release Date | Description                                 |
----------|--------------|---------------------------------------------|
| 1 | August 21, 2022 | Initial version. |
