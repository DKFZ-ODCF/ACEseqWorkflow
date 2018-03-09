/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

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


import static de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider.*;

/**
 * Created by kleinhei on 6/16/14.
 */
@groovy.transform.CompileStatic
@StaticScriptProviderClass
public final class ACESeqMethods {

    static LinkedHashMap<String, String> getGlobalJobSpecificParameters(Configuration config) {
        return new LinkedHashMap<String, String>(
                (ACEseqConstants.CHR_NAME): config.configurationValues.getString(ACEseqConstants.CHR_NAME),
                (ACEseqConstants.CHR_NR): config.configurationValues.getString(ACEseqConstants.CHR_NR),
                (ACEseqConstants.GENETIC_MAP_FILE): config.configurationValues.getString(ACEseqConstants.GENETIC_MAP_FILE),
                (ACEseqConstants.KNOWN_HAPLOTYPES_FILE): config.configurationValues.getString(ACEseqConstants.KNOWN_HAPLOTYPES_FILE),
                (ACEseqConstants.KNOWN_HAPLOTYPES_LEGEND_FILE): config.configurationValues.getString(ACEseqConstants.KNOWN_HAPLOTYPES_LEGEND_FILE))
    }

    static CnvSnpGeneratorResultByType generateCNVSNPs(BamFile controlBam, BamFile tumorBam) {
        IndexedFileObjects indexedFileObjects = ParallelizationHelper.runParallel(
                COConstants.CVALUE_CHROMOSOME_INDICES,
                ACEseqConstants.TOOL_CNV_SNP_GENERATION,
                tumorBam,
                controlBam,
                ACEseqConstants.PARM_CHR_INDEX,
                getGlobalJobSpecificParameters(controlBam.executionContext.configuration))
        return new CnvSnpGeneratorResultByType(indexedFileObjects, controlBam.getExecutionContext())
    }

    static ImputeGenotypeByChromosome imputeGenotypes(BamFile controlBam) {
        IndexedFileObjects indexedFileObjects = ParallelizationHelper.runParallel(
                COConstants.CVALUE_AUTOSOME_INDICES,
                ACEseqConstants.TOOL_IMPUTE_GENOTYPES,
                controlBam,
                null,
                ACEseqConstants.PARM_CHR_INDEX,
                getGlobalJobSpecificParameters(controlBam.executionContext.configuration))
        return new ImputeGenotypeByChromosome(indexedFileObjects, controlBam.getExecutionContext())
    }

    static ImputeGenotypeByChromosome imputeGenotypes(UnphasedGenotypeFileGroupByChromosome unphasedGenotypeFiles) {
        Map<String, UnphasedGenotypeFile> mapOfFiles = [:]
        mapOfFiles += unphasedGenotypeFiles.getFiles()
        mapOfFiles.remove("X")
        IndexedFileObjects indexedFileObjects = runParallel(
                ACEseqConstants.TOOL_IMPUTE_GENOTYPES_NOMPILEUP,
                new UnphasedGenotypeFileGroupByChromosome(mapOfFiles.keySet() as List<String>, mapOfFiles, unphasedGenotypeFiles.getExecutionContext()),
                null,
                ACEseqConstants.PARM_CHR_INDEX,
                getGlobalJobSpecificParameters(unphasedGenotypeFiles.executionContext.configuration))
        return new ImputeGenotypeByChromosome(indexedFileObjects, unphasedGenotypeFiles.getExecutionContext())
    }

    static Tuple2<PhasedGenotypeFile, HaploblockGroupFile> imputeGenotypeX(TextFile sexFile, BamFile controlBam) {
        return (Tuple2<PhasedGenotypeFile, HaploblockGroupFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_IMPUTE_GENOTYPEX, controlBam, sexFile);
    }

    static Tuple2<PhasedGenotypeFile, HaploblockGroupFile> imputeGenotypeX(TextFile sexFile, UnphasedGenotypeFileGroupByChromosome unphasedGenotypeFiles) {
        return (Tuple2<PhasedGenotypeFile, HaploblockGroupFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_IMPUTE_GENOTYPEX_NOMPILEUP, unphasedGenotypeFiles.getFiles().get("X"), sexFile);
    }

    static TextFile getGenotypes(TextFile mergedAndFilteredSNPFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GET_GENOTYPES, mergedAndFilteredSNPFile);
    }

    static UnphasedGenotypeFileGroupByChromosome createUnphased(TextFile genotypeSNPFile) {
        Map<String, UnphasedGenotypeFile> listOfFiles = new LinkedHashMap<>();
        List<BaseFile> filesToCheck = new LinkedList<>();
        List<String> keyset = Arrays.asList("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X");
        for (String chrIndex : keyset) {
            UnphasedGenotypeFile unphasedGenotypeFile = (UnphasedGenotypeFile) BaseFile.constructManual(UnphasedGenotypeFile.class, genotypeSNPFile);
            String path = unphasedGenotypeFile.getAbsolutePath();
            unphasedGenotypeFile.setPath(new File(path.replace("#CHROMOSOME_INDEX#", chrIndex)));
            listOfFiles.put(chrIndex, unphasedGenotypeFile);
            filesToCheck.add(unphasedGenotypeFile);
        }

//        TextFile unphasedGenotypeChrXFile = (TextFile)BaseFile.constructManual(TextFile.class, files.get("1"));
//        UnphasedGenotypeChrXFile.overrideFilenameUsingSelectionTag("unphasedGenotypeChrXFile");
//        filesToCheck.add(unphasedGenotypeFile);

        ExecutionContext context = filesToCheck[0].getExecutionContext();
        Map<String, Object> parameters = context.getDefaultJobParameters(ACEseqConstants.TOOL_CREATE_UNPHASED_GENOTYPE);
        parameters["FILENAME_SNP_POSITIONS_WG_FAKE"] = genotypeSNPFile.getAbsolutePath();
        parameters.put("FILENAME_UNPHASED_GENOTYPE", "( " + filesToCheck.collect { BaseFile file -> file.getAbsolutePath() }.join(" ") + ' )');

        Job job = new Job(context, context.createJobName((UnphasedGenotypeFile) listOfFiles.get("1"), ACEseqConstants.TOOL_CREATE_UNPHASED_GENOTYPE,
                true), ACEseqConstants.TOOL_CREATE_UNPHASED_GENOTYPE, null, parameters,
                [genotypeSNPFile] as List<BaseFile>, filesToCheck)
        BEJobResult jobResult = job.run();
        for (BaseFile baseFile : filesToCheck) {
            baseFile.setCreatingJobsResult(jobResult);
        }

        return new UnphasedGenotypeFileGroupByChromosome(keyset, listOfFiles, genotypeSNPFile.getExecutionContext());
    }


    static TextFile addHaploTypes(TextFile mergedAndFilteredSNPFile, PhasedGenotypeFileGroupByChromosome phasedGenotypes, PhasedGenotypeFile phasedGenotypeX) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ADD_HAPLOTYPES_TO_SNP_FILE, mergedAndFilteredSNPFile, phasedGenotypes, phasedGenotypeX);
    }

    static TextFile createControlBafPlot(TextFile haplotypedSNPFile, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_CREATE_CONTROL_BAF_PLOTS, haplotypedSNPFile, genderFile);
    }

    static TextFile replaceControl(TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool("replaceBadControl", genderFile);
    }

    static Tuple3<TextFile, TextFile, TextFile> correctGC(TextFile mergedAndFilteredCovWinFile) {
        return (Tuple3<TextFile, TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_CORRECT_GC_BIAS, mergedAndFilteredCovWinFile);
    }

    static Tuple2<TextFile, TextFile> pscbsGaps(TextFile haplotypedSNPFile, TextFile correctedCovWinFile, TextFile genderFile) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GET_BREAKPOINTS, haplotypedSNPFile, correctedCovWinFile, genderFile);
    }

    static BEJobResult getFileExistedFakeJobResult() {
        return new BEJobResult(null, new FakeBEJob(new BEFakeJobID(BEFakeJobID.FakeJobReason.FILE_EXISTED)), null, null, null, null)
    }

    @ScriptCallingMethod
    static Tuple2<TextFile, TextFile> mergeSv(TextFile knownSegmentsFile) {
        TextFile svFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "svFileTag", null, null);
        svFile.setAsSourceFile();
        BEJobResult result = getFileExistedFakeJobResult()
        svFile.setCreatingJobsResult(result);

        // TODO: The following checks should be done by ACEseqWorkflow.checkExecutability(). Unfortunately currently it is not trivial to get the name of the file whose existence to check.
        boolean b = FileSystemAccessProvider.getInstance().checkBaseFiles(svFile);
        if (b)
            return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MERGE_BREAKPOINTS_AND_SV, knownSegmentsFile, svFile);

        knownSegmentsFile.getExecutionContext().addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("SV files were not found in input path."));
        return null;
    }

    @ScriptCallingMethod
    static Tuple2<TextFile, TextFile> mergeNoSv(TextFile knownSegmentsFile) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MERGE_BREAKPOINTS_WITHOUT_SV, knownSegmentsFile)
    }

    @ScriptCallingMethod
    static Tuple2<TextFile, TextFile> mergeCrest(TextFile knownSegmentsFile) {
        TextFile svFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "crestDelDupInvFileTag", null, null)
        TextFile translocFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "crestTranslocFileTag", null, null)

        boolean b = FileSystemAccessProvider.getInstance().checkBaseFiles(svFile, translocFile)
        if (!b) {
            knownSegmentsFile.getExecutionContext().addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("Crest files were not found in input path."))
            return null
        }

        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MERGE_BREAKPOINTS_AND_SV_CREST, knownSegmentsFile, svFile, translocFile)
    }

    @ScriptCallingMethod
    static TextFile getSegmentAndGetSnps(TextFile breaks, TextFile pscbsSnps) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GET_SEGMENTS_AND_SNPS, breaks, pscbsSnps);
    }

    @ScriptCallingMethod
    static TextFile markSegsWithHomozygDel(TextFile segments, TextFile svPoints) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MARK_HOMOZYGOUS_DELETIONS, segments, svPoints);
    }

    @ScriptCallingMethod
    static TextFile segsToSnpDataHomodel(TextFile segments, TextFile pscbsSnpsFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_SEGMENTS_TO_SNP_DATA_HOMODEL, segments, pscbsSnpsFile);
    }

    @ScriptCallingMethod
    static Tuple2<TextFile, TextFile> clusterPruneSegments(TextFile segmentsFile, TextFile snpsFile, TextFile genderFile, TextFile correctParams) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_CLUSTER_AND_PRUNE_SEGMENTS, segmentsFile, snpsFile, genderFile, correctParams);
    }

    @ScriptCallingMethod
    static TextFile segsToSnpDataPruned(TextFile segments, TextFile snpsFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_SEGMENTS_TO_SNP_DATA_PRUNED, segments, snpsFile);
    }

    @ScriptCallingMethod
    static TextFile estimatePeaks(TextFile segments, TextFile snpsFile, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ESTIMATE_PEAKS_FOR_PURITY, segments, snpsFile, genderFile);
    }

    @ScriptCallingMethod
    static TextFile estimatePurityPloidy(TextFile segments, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ESTIMATE_PURITY_AND_PLOIDY, segments, genderFile);
    }

    @ScriptCallingMethod
    static Tuple2<TextFile, TextFile> generatePlots(TextFile segments, TextFile snpsFile, TextFile svPoints, TextFile purityPloidyFile, TextFile genderFile) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GENERATE_RESULTS_AND_PLOTS, segments, snpsFile, svPoints, purityPloidyFile, genderFile);
    }

    @ScriptCallingMethod
    static TextFile convertToVcf(TextFile purityPloidyFile, TextFile checkpointFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GENERATE_VCF_FROM_TAB, purityPloidyFile, checkpointFile);
    }

    @ScriptCallingMethod
    static TextFile estimateHRD(TextFile genderFile, TextFile cnvParameterFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ESTIMATE_HRD_SCORE, genderFile, cnvParameterFile);
    }

    static IndexedFileObjects runParallel(String toolID, IndexedFileObjects fileGroup, BaseFile otherFile, String indexParameterName, LinkedHashMap<String, String> parameters = [:]) {
        List<String> indices = fileGroup.getIndices();

        //First one executes locally or via ssh but without a cluster system.
        def stream = Roddy.jobManager.executesWithoutJobSystem() ? indices.parallelStream() : indices.stream();
        Map<String, FileObject> map = stream.collect { String index ->
            LinkedHashMap<String, String> indexMap = new LinkedHashMap((indexParameterName): index)
            indexMap.putAll(parameters)
            new MapEntry(index, ParallelizationHelper.callWithOptionalSecondaryBam(
                    toolID,
                    (BaseFile) fileGroup.getIndexedFileObjects().get(index),
                    otherFile,
                    indexMap)
            )
        } as Map<String, FileObject>

        return new IndexedFileObjects(indices, map, fileGroup.getExecutionContext());
    }

}
