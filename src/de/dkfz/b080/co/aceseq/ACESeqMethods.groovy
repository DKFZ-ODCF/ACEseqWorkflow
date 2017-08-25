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
@groovy.transform.CompileStatic
@StaticScriptProviderClass
public final class ACESeqMethods {

    public static CnvSnpGeneratorResultByType generateCNVSNPs(BamFile controlBam, BamFile tumorBam) {
        IndexedFileObjects indexedFileObjects = ParallelizationHelper.runParallel(COConstants.CVALUE_CHROMOSOME_INDICES, ACEseqConstants.TOOL_CNV_SNP_GENERATION, tumorBam, controlBam, "PARM_CHR_INDEX=");
        return new CnvSnpGeneratorResultByType(indexedFileObjects, controlBam.getExecutionContext());
    }

    public static ImputeGenotypeByChromosome imputeGenotypes(BamFile controlBam) {
        IndexedFileObjects indexedFileObjects = ParallelizationHelper.runParallel(COConstants.CVALUE_AUTOSOME_INDICES, ACEseqConstants.TOOL_IMPUTE_GENOTYPES, controlBam, null, "PARM_CHR_INDEX=");
        return new ImputeGenotypeByChromosome(indexedFileObjects, controlBam.getExecutionContext());
    }

    public static ImputeGenotypeByChromosome imputeGenotypes( UnphasedGenotypeFileGroupByChromosome unphasedGenotypeFiles) {
	Map<String, UnphasedGenotypeFile> mapOfFiles = [:]
	mapOfFiles += unphasedGenotypeFiles.getFiles();
	mapOfFiles.remove("X")
	IndexedFileObjects indexedFileObjects = runParallel(ACEseqConstants.TOOL_IMPUTE_GENOTYPES_NOMPILEUP, new UnphasedGenotypeFileGroupByChromosome(mapOfFiles.keySet() as List<String>, mapOfFiles, unphasedGenotypeFiles.getExecutionContext()), null, "PARM_CHR_INDEX=");
        return new ImputeGenotypeByChromosome(indexedFileObjects, unphasedGenotypeFiles.getExecutionContext());
    }

    public static Tuple2<PhasedGenotypeFile, HaploblockGroupFile> imputeGenotypeX(TextFile sexFile, BamFile controlBam) {
        return (Tuple2<PhasedGenotypeFile, HaploblockGroupFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_IMPUTE_GENOTYPEX, controlBam, sexFile);
    }

    public static Tuple2<PhasedGenotypeFile, HaploblockGroupFile> imputeGenotypeX(TextFile sexFile, UnphasedGenotypeFileGroupByChromosome unphasedGenotypeFiles ) {
        return (Tuple2<PhasedGenotypeFile, HaploblockGroupFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_IMPUTE_GENOTYPEX_NOMPILEUP, unphasedGenotypeFiles.getFiles().get("X"), sexFile);
    }

    public static TextFile getGenotypes( TextFile mergedAndFilteredSNPFile ) {
	return (TextFile) GenericMethod.callGenericTool( ACEseqConstants.TOOL_GET_GENOTYPES, mergedAndFilteredSNPFile );
    }

    public static UnphasedGenotypeFileGroupByChromosome createUnphased( TextFile genotypeSNPFile ) {
        Map<String, UnphasedGenotypeFile> listOfFiles = new LinkedHashMap<>();
        List<BaseFile> filesToCheck = new LinkedList<>();
	List<String> keyset =Arrays.asList("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X");
        for (String chrIndex : keyset) {
            UnphasedGenotypeFile unphasedGenotypeFile = (UnphasedGenotypeFile)BaseFile.constructManual(UnphasedGenotypeFile.class, genotypeSNPFile);
            String path = unphasedGenotypeFile.getAbsolutePath();
            unphasedGenotypeFile.setPath(new File(path.replace("#CHROMOSOME_INDEX#", chrIndex)));
            listOfFiles.put(chrIndex, unphasedGenotypeFile);
            filesToCheck.add(unphasedGenotypeFile);
        }

//        TextFile unphasedGenotypeChrXFile = (TextFile)BaseFile.constructManual(TextFile.class, files.get("1"));
//        UnphasedGenotypeChrXFile.overrideFilenameUsingSelectionTag("unphasedGenotypeChrXFile");
//        filesToCheck.add(unphasedGenotypeFile);

        ExecutionContext run = filesToCheck[0].getExecutionContext();
        Map<String, Object> parameters = run.getDefaultJobParameters(ACEseqConstants.TOOL_CREATE_UNPHASED_GENOTYPE);
	parameters["FILENAME_SNP_POSITIONS_WG_FAKE"]=genotypeSNPFile.getAbsolutePath();
        parameters.put("FILENAME_UNPHASED_GENOTYPE", "( " + filesToCheck.collect { BaseFile file -> file.getAbsolutePath() }.join(" ") + ' )');

	Job job = new Job(run, run.createJobName((UnphasedGenotypeFile)listOfFiles.get("1"), ACEseqConstants.TOOL_CREATE_UNPHASED_GENOTYPE, true), ACEseqConstants.TOOL_CREATE_UNPHASED_GENOTYPE, null, parameters, [ genotypeSNPFile ] as List<BaseFile>, filesToCheck);
        JobResult jobResult = job.run();
        for (BaseFile baseFile : filesToCheck) {
            baseFile.setCreatingJobsResult(jobResult);
        }

        return new UnphasedGenotypeFileGroupByChromosome(keyset, listOfFiles, genotypeSNPFile.getExecutionContext());
    }


    public static TextFile addHaploTypes(TextFile mergedAndFilteredSNPFile, PhasedGenotypeFileGroupByChromosome phasedGenotypes, PhasedGenotypeFile phasedGenotypeX) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ADD_HAPLOTYPES_TO_SNP_FILE, mergedAndFilteredSNPFile, phasedGenotypes, phasedGenotypeX);
    }

    public static TextFile createControlBafPlot(TextFile haplotypedSNPFile , TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool( ACEseqConstants.TOOL_CREATE_CONTROL_BAF_PLOTS, haplotypedSNPFile, genderFile );
    }

    public static TextFile replaceControl(TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool("replaceBadControl", genderFile);
    }

    public static Tuple3<TextFile, TextFile, TextFile> correctGC(TextFile mergedAndFilteredCovWinFile) {
        return (Tuple3<TextFile, TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_CORRECT_GC_BIAS, mergedAndFilteredCovWinFile);
    }

    public static Tuple2<TextFile, TextFile> pscbsGaps(TextFile haplotypedSNPFile, TextFile correctedCovWinFile, TextFile genderFile) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GET_BREAKPOINTS, haplotypedSNPFile, correctedCovWinFile, genderFile);
    }

    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> mergeSv(TextFile knownSegmentsFile) {
        TextFile svFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "svFileTag", null, null);
        svFile.setAsSourceFile();
        JobResult result = new JobResult(knownSegmentsFile.getExecutionContext(), null, JobDependencyID.getFileExistedFakeJob(knownSegmentsFile.getExecutionContext()), false, null, null, null);
        svFile.setCreatingJobsResult(result);
        boolean b = FileSystemAccessProvider.getInstance().checkBaseFiles(svFile);
        if (b)
            return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MERGE_BREAKPOINTS_AND_SV, knownSegmentsFile, svFile);

        knownSegmentsFile.getExecutionContext().addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("SV files were not found in input path."));
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

        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MERGE_BREAKPOINTS_AND_SV_CREST, knownSegmentsFile, svFile, translocFile);
    }

    @ScriptCallingMethod
    public static TextFile getSegmentAndGetSnps(TextFile breaks, TextFile pscbsSnps) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GET_SEGMENTS_AND_SNPS, breaks, pscbsSnps);
    }

    @ScriptCallingMethod
    public static TextFile markSegsWithHomozygDel(TextFile segments, TextFile svPoints) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MARK_HOMOZYGOUS_DELETIONS, segments, svPoints);
    }

    @ScriptCallingMethod
    public static TextFile segsToSnpDataHomodel(TextFile segments, TextFile pscbsSnpsFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_SEGMENTS_TO_SNP_DATA_HOMODEL, segments, pscbsSnpsFile);
    }

    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> clusterPruneSegments(TextFile segmentsFile, TextFile snpsFile, TextFile genderFile, TextFile correctParams) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_CLUSTER_AND_PRUNE_SEGMENTS, segmentsFile, snpsFile, genderFile, correctParams);
    }

    @ScriptCallingMethod
    public static TextFile segsToSnpDataPruned(TextFile segments, TextFile snpsFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_SEGMENTS_TO_SNP_DATA_PRUNED, segments, snpsFile);
    }

    @ScriptCallingMethod
    public static TextFile estimatePeaks(TextFile segments, TextFile snpsFile, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ESTIMATE_PEAKS_FOR_PURITY, segments, snpsFile, genderFile);
    }

    @ScriptCallingMethod
    public static TextFile estimatePurityPloidy(TextFile segments, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ESTIMATE_PURITY_AND_PLOIDY, segments, genderFile);
    }

    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> generatePlots(TextFile segments, TextFile snpsFile, TextFile svPoints, TextFile purityPloidyFile, TextFile genderFile) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GENERATE_RESULTS_AND_PLOTS, segments, snpsFile, svPoints, purityPloidyFile, genderFile);
    }

    @ScriptCallingMethod
    public static TextFile estimateHRD(TextFile genderFile, TextFile cnvParameterFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ESTIMATE_HRD_SCORE, genderFile, cnvParameterFile);
    }

    public static IndexedFileObjects runParallel( String toolID, IndexedFileObjects fileGroup, BaseFile otherFile, String indexParameterName) {
        List<String> indices = fileGroup.getIndices();
        Map<String, FileObject> map = new LinkedHashMap<>();

        //First one executes locally or via ssh but without a cluster system.
        def stream = JobManager.getInstance().executesWithoutJobSystem() ? indices.parallelStream() : indices.stream();
        stream.each{String index -> ParallelizationHelper.callWithIndex(toolID, index, indexParameterName, map, (BaseFile) fileGroup.getIndexedFileObjects().get(index), otherFile)};

        return new IndexedFileObjects(indices, map, fileGroup.getExecutionContext());
    }

}
