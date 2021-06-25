Avalanche Rich Profile:

Purpose: Detect where along avalanches members of the rich club are more active and if that changes for avalanches of different durations.

Sequence of running the scripts:

			1) AvalRichProfileRev_DoItAll.m --> Using "asdf," "asdf_shuf," "PDF" and "wgt" files, this script detects avalanches,
							    detects where along those avalanches the RC members are more active (proportion 
							    of active members of the RC) and plots the proportions for the assigned avalanche 
							    durations.
						
							    The script could also load the previously generated "aval_list" flies
							    containing lists of avalanches detected in "asdf" files. This makes 
							    the process faster (especially when doing the part for shuffled data sets).