package de.dkfz.b080.co.aceseq;

import de.dkfz.b080.co.common.ParallelizationHelper;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.Roddy;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.core.ExecutionContextError;
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider;
import de.dkfz.roddy.execution.jobs.*;
import de.dkfz.roddy.knowledge.files.*;
import de.dkfz.roddy.knowledge.methods.GenericMethod;


import javax.xml.soap.Text;
import java.io.File;
import java.nio.file.spi.FileSystemProvider;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider.*;

/**
 * Created by kleinhei on 6/16/14.
 */
@StaticScriptProviderClass
public final class ACESeqMethods {

    @ScriptCallingMethod
    public static CnvSnpGeneratorResultByType generateCNVSNPs(BamFile controlBam, BamFile tumorBam) {
        IndexedFileObjects indexedFileObjects = ParallelizationHelper.runParallel(COConstants.CVALUE_CHROMOSOME_INDICES, COConstants.TOOL_CNV_SNP_GENERATION, tumorBam, controlBam, "PARM_CHR_INDEX=");
        return new CnvSnpGeneratorResultByType(indexedFileObjects, controlBam.getExecutionContext());
    }

    public static CnvSnpGeneratorResultByType generateCNVSNPsMeta(BamFile controlBam, BamFile tumorBam) {
        return GenericMethod.callGenericTool("cnvSnpGenerationMeta", tumorBam, controlBam);
    }


    @ScriptCallingMethod
    public static ImputeGenotypeByChromosome imputeGenotypes(BamFile controlBam) {
        IndexedFileObjects indexedFileObjects = ParallelizationHelper.runParallel(COConstants.CVALUE_AUTOSOME_INDICES, COConstants.TOOL_IMPUTE_GENOTYPES, controlBam, null, "PARM_CHR_INDEX=");
        return new ImputeGenotypeByChromosome(indexedFileObjects, controlBam.getExecutionContext());
    }

    @ScriptCallingMethod
    public static Tuple2<PhasedGenotypeFile, HaploblockGroupFile> imputeGenotypeX(TextFile sexFile, BamFile controlBam) {
        return (Tuple2<PhasedGenotypeFile, HaploblockGroupFile>) GenericMethod.callGenericTool(COConstants.TOOL_IMPUTE_GENOTYPEX, controlBam, sexFile);
    }


    @ScriptCallingMethod
    public static TextFile addHaploTypes(TextFile mergedAndFilteredSNPFile, PhasedGenotypeFileGroupByChromosome phasedGenotypes, PhasedGenotypeFile phasedGenotypeX) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_ADD_HAPLOTYPES_TO_SNP_FILE, mergedAndFilteredSNPFile, phasedGenotypes, phasedGenotypeX);
    }

    public static TextFile replaceControl(TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool("replaceBadControl", genderFile);
    }

    @ScriptCallingMethod
    public static Tuple3<TextFile, TextFile, TextFile> correctGC(TextFile mergedAndFilteredCovWinFile) {
        return (Tuple3<TextFile, TextFile, TextFile>) GenericMethod.callGenericTool(COConstants.TOOL_CORRECT_GC_BIAS, mergedAndFilteredCovWinFile);
    }

    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> pscbsGaps(TextFile haplotypedSNPFile, TextFile correctedCovWinFile, TextFile genderFile) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(COConstants.TOOL_GET_BREAKPOINTS, haplotypedSNPFile, correctedCovWinFile, genderFile);
    }

    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> mergeDelly(TextFile knownSegmentsFile) {
        TextFile svFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "dellyFileTag", null, null);
        svFile.setAsSourceFile();
        JobResult result = new JobResult(knownSegmentsFile.getExecutionContext(), null, JobDependencyID.getFileExistedFakeJob(knownSegmentsFile.getExecutionContext()), false, null, null, null);
        svFile.setCreatingJobsResult(result);
        boolean b = FileSystemAccessProvider.getInstance().checkBaseFiles(svFile);
        if (b)
            return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(COConstants.TOOL_MERGE_BREAKPOINTS_AND_SV_DELLY, knownSegmentsFile, svFile);

        knownSegmentsFile.getExecutionContext().addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("Delly files were not found in input path."));
        return null;
    }

    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> mergeCrest(TextFile knownSegmentsFile) {
        TextFile svFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "crestDelDupInvFileTag", null, null);
        TextFile translocFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "crestTranslocFileTag", null, null);

        boolean b = FileSystemAccessProvider.getInstance().checkBaseFiles(svFile, translocFile);
        if (!b) {
            knownSegmentsFile.getExecutionContext().addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("Crest files were not found in input path."));
            return null;
        }

        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(COConstants.TOOL_MERGE_BREAKPOINTS_AND_SV_CREST, knownSegmentsFile, svFile, translocFile);
    }

    @ScriptCallingMethod
    public static TextFile getSegmentAndGetSnps(TextFile breaks, TextFile pscbsSnps) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_GET_SEGMENTS_AND_SNPS, breaks, pscbsSnps);
    }

    @ScriptCallingMethod
    public static TextFile markSegsWithHomozygDel(TextFile segments, TextFile svPoints) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_MARK_HOMOZYGOUS_DELETIONS, segments, svPoints);
    }

    @ScriptCallingMethod
    public static TextFile segsToSnpDataHomodel(TextFile segments, TextFile pscbsSnpsFile) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_SEGMENTS_TO_SNP_DATA_HOMODEL, segments, pscbsSnpsFile);
    }

    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> clusterPruneSegments(TextFile segmentsFile, TextFile snpsFile, TextFile genderFile, TextFile correctParams, HaploblockFileGroupByChromosome haplogroups, HaploblockGroupFile haplogroupsX) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(COConstants.TOOL_CLUSTER_AND_PRUNE_SEGMENTS, segmentsFile, snpsFile, genderFile, correctParams, haplogroups, haplogroupsX);
    }

    @ScriptCallingMethod
    public static TextFile segsToSnpDataPruned(TextFile segments, TextFile snpsFile) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_SEGMENTS_TO_SNP_DATA_PRUNED, segments, snpsFile);
    }

    @ScriptCallingMethod
    public static TextFile estimatePeaks(TextFile segments, TextFile snpsFile, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_ESTIMATE_PEAKS_FOR_PURITY, segments, snpsFile, genderFile);
    }

    @ScriptCallingMethod
    public static TextFile estimatePurityPloidy(TextFile segments, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_ESTIMATE_PURITY_AND_PLOIDY, segments, genderFile);
    }

    @ScriptCallingMethod
    public static TextFile generatePlots(TextFile segments, TextFile snpsFile, TextFile svPoints, TextFile purityPloidyFile, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_GENERATE_RESULTS_AND_PLOTS, segments, snpsFile, svPoints, purityPloidyFile, genderFile);
    }

    @ScriptCallingMethod
    public static TextFile convertToVcf(TextFile checkpointFile) {
        return (TextFile) GenericMethod.callGenericTool(COConstants.TOOL_GENERATE_VCF_FROM_TAB, checkpointFile);
    }

}
