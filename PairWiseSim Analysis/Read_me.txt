Pair-wise Similarity Analysis

Purpose: Calculating pair-wise similarity values (Jaccard coefficients) for all pairs of neurons within a data set and 
	 examining if there is a difference between Jaccard coefficients of neurons that are members of the rich club,
	 and that of other neurons.

Sequence of running the scripts:

			1) fr_sim_analysis.m --> Calculates the Jaccard coefficients for all possible pairs of neurons and
						 100 shuffled versions (ISIs are shuffled in asdf files) and saves corresponding files.

			2) jittered pairwise sim?

			2) PairSimFigs.m --> Plots histograms and probability distributions of Jaccard coefficients based on 
					     neurons' memberships to the RC. Also plots similarity matrices for original,
					     shuffled and jittered data sets sorted in descending order of out-degrees of neurons.