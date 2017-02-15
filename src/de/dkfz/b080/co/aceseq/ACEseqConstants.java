package de.dkfz.b080.co.aceseq;

/**
 * A static class for the storage of computational oncology based constants.
 */
@groovy.transform.CompileStatic
public final class ACEseqConstants {
    /**
     * Tool entries
     */
    public static final String TOOL_CNV_SNP_GENERATION = "cnvSnpGeneration";
    public static final String TOOL_REPLACE_BAD_CONTROL = "replaceBadControl";
    public static final String TOOL_MERGE_AND_FILTER_SNP_FILES = "mergeAndFilterSnpFiles";
    public static final String TOOL_MERGE_AND_FILTER_CNV_FILES = "mergeAndFilterCnvFiles";
    public static final String TOOL_ANNOTATE_COV_WIN = "annotateCnvFiles";
    public static final String TOOL_GET_GENOTYPES = "getGenotypes";
    public static final String TOOL_CREATE_UNPHASED_GENOTYPE = "createUnphasedGenotype";
    public static final String TOOL_IMPUTE_GENOTYPES = "imputeGenotypes";
    public static final String TOOL_IMPUTE_GENOTYPES_NOMPILEUP = "imputeGenotypes_noMpileup";
    public static final String TOOL_IMPUTE_GENOTYPEX = "imputeGenotypes_X";
    public static final String TOOL_IMPUTE_GENOTYPEX_NOMPILEUP = "imputeGenotypes_X_noMpileup";
    public static final String TOOL_ADD_HAPLOTYPES_TO_SNP_FILE = "addHaplotypesToSnpFile";
    public static final String TOOL_CORRECT_GC_BIAS = "correctGcBias";
    public static final String TOOL_GET_BREAKPOINTS = "getBreakpoints";
    public static final String TOOL_MERGE_BREAKPOINTS_AND_SV_DELLY = "mergeBreakpointsAndSvDelly";
    public static final String TOOL_MERGE_BREAKPOINTS_AND_SV_CREST = "mergeBreakpointsAndSvCrest";
    public static final String TOOL_GET_SEGMENTS_AND_SNPS = "getSegmentsAndSnps";
    public static final String TOOL_MARK_HOMOZYGOUS_DELETIONS = "markHomozygousDeletions";
    public static final String TOOL_SEGMENTS_TO_SNP_DATA_HOMODEL = "segmentsToSnpDataHomodel";
    public static final String TOOL_CLUSTER_AND_PRUNE_SEGMENTS = "clusterAndPruneSegments";
    public static final String TOOL_SEGMENTS_TO_SNP_DATA_PRUNED = "segmentsToSnpDataPruned";
    public static final String TOOL_ESTIMATE_PEAKS_FOR_PURITY = "estimatePeaksForPurity";
    public static final String TOOL_ESTIMATE_PURITY_AND_PLOIDY = "estimatePurityAndPloidy";
    public static final String TOOL_GENERATE_RESULTS_AND_PLOTS = "generateResultsAndPlots";
    public static final String TOOL_GENERATE_VCF_FROM_TAB = "generateVcfFromTab";
    
    private ACEseqConstants() {
    }
}
