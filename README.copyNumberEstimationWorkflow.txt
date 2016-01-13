== Description

Some description

== Run flags / switches

Switch                      Default Description
runWithDelly          		true    Run with delly output (from EMBLs delly workflow) or with crest
runQualityCheckOnly      	false   Stop the workflow after the gc correction.
runMetaCNVGeneration        false   Run a meta step for the cnv generation. This submits a single,
									slightly optimized step instead of several (per chromosome) jobs.

== Changelist

* Version update to 1.0.189


pscbs_all.R: bugfix that screwed up chrlength

* Version update to 1.0.187


-stabilizing addition to pscbs_all.R

* Version update to 1.0.185

-bugfixes for pscbs_all.R

* Version update to 1.0.183






-pscbs_all.R and PSCBSall.sh:
	*scientific format of start and end here already, replace +Inf/-Inf with chromosome length/0
-adjustsAlleleFreqs and functions.R and purity_ploidy.R:
	*don't consider dbSNP position with 0 reads mapped in tumor
-manual_pruning.R
	*don't consider dbSNP position with 0 reads mapped in tumor and remove bug in case no main cluster is found


* Version update to 1.0.181


-convertTabTovcf.sh:
	*removed usage of id option (for pcawg output)
-convertToVcf.py:
	*libraries removed
	* default for id argument set (NA)
	*changed version number (for pcawg output)
	* added missing "\n" for header
	* added SAMPLE_ID line for header (pancan)
	*changed sample_$pid column name to TUMOR
-correctGCBias_functions.R:
	*corrected coordinates in coverage plots
-datatablePSCBSgaps.sh:
	*tabix for position not window (col5 instead of 1)
-haplotypes.sh:
	*tabix corrected (added -e)
-manual_pruning.R:
	*adjusted for no cluster within cluster limits
-merge_and_filter_cnv.py/merge_and_filter_snp.py:
	*removed coverage filter for tumor
-pscbs_all.R:
	*round coordinates with .5 (ceiling for start, floor for end)
-PSCBSgabs_plus_delly_points.py:
	*removed library
-pscbs_plots_functions.R:
	*fullPloidy calculation by largest fraction of genome instead of most segments in genome
	*annotation of LOH (cn,gain,loss) corrected
	*LOH definition adjusted instead of dhMean </> 0.8 using round(c1) and round(c2) (==0 for LOHs, != 0 for everything else)
-pscbs_plots.R:
	*added gc()
-purity_ploidy.R:
	*average coverage on autosomes instead of all chromosomes (better for males)
	*allow missing autosomes
-vcfAnno.sh:
	*rearranged order of script (first gender estimation)

-addition of fake control option:
	added scripts:
		*replaceControl.sh
		*replaceControlACEseq.R
 

* Version update to 1.0.158

extracted into single plugin

* Version update to 1.0.131

- plots					tab seperateed file with cnv parameters written
- correctGCBias:			extra file with qc parameters printed, conversion to json, parameter for scale adjustion added
- Change workflow class to override another execute method. This makes the workflow a bit cleaner.

* Version update to 1.0.114

-PSCBSgabs_plus_CRESTpoints.py :	add identifier column to sv_points file (used with DELLY calls)
-PSCBSgabs_plus_delly_points.py:	add identifier column to sv_points file, convert start coordinates to 1-based
-homozygous_deletion.pl:			added id column for sv file, adjust length calculation to inclusion of start AND End coordinate
-correctGCBias.R:					density estimation without restrictions for final corrected covRatio
-manual_pruning.R:					check for haploblock files to prevent script from failing after 2 hours due to missing files, plot Names with PID,
									new feature: merge clusters again after outlier removal and choose random cluster if 2 or more "mainClusters are found"
-purity_ploidy.R:					estimate and report average coverage, add minCoverage as optional parameter
-purity_ploidy_estimation_final.R:	bug fixes for balanced segments using wrong matrix,
									new feature: disallow copy number states of 0, don't punish balanced segments with copy number down to -0.3
-functions.R:						soft limits for control SNPs used for peak calling (coverage dependant), no limits for tumor SNPs, single threat for running
-purityPloidity_EstimateFinal.sh:	add PID as parameter
-analysisCopyNumberEstimation.xml:	adjust settings for peak calling to actually used resources
-pscbs_plots_functions.R:			convert 0.5 coordinates to integer values, don't allow scientific format for output files

* Version update to 1.0.109

* Version update to 1.0.105

* Version update to 1.0.104

-bug fix: read all lines of breakpoint.txt to avoid missing information on chr 1 
* Version update to 1.0.103