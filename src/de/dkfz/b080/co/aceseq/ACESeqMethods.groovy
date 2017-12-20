/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */
package de.dkfz.b080.co.aceseq

import de.dkfz.b080.co.common.ParallelizationHelper
import de.dkfz.b080.co.files.*
import de.dkfz.roddy.Roddy
import de.dkfz.roddy.config.Configuration
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.core.ExecutionContextError
import de.dkfz.roddy.execution.io.fs.FileSystemAccessProvider
import de.dkfz.roddy.execution.jobs.*
import de.dkfz.roddy.knowledge.files.*
import de.dkfz.roddy.knowledge.methods.GenericMethod
import de.dkfz.roddy.plugins.LibrariesFactory

import java.util.logging.Level

/**
 * Created by kleinhei on 6/16/14.
 */
@groovy.transform.CompileStatic
@StaticScriptProviderClass
final class ACESeqMethods {

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
                ACEseqConstants.PARM_CHR_INDEX)
//        ,
//                getGlobalJobSpecificParameters(controlBam.executionContext.configuration))
        return new CnvSnpGeneratorResultByType(indexedFileObjects, controlBam.getExecutionContext())
    }

    static ImputeGenotypeByChromosome imputeGenotypes(BamFile controlBam) {
        IndexedFileObjects indexedFileObjects = ParallelizationHelper.runParallel(
                COConstants.CVALUE_AUTOSOME_INDICES,
                ACEseqConstants.TOOL_IMPUTE_GENOTYPES,
                controlBam,
                null,
                ACEseqConstants.PARM_CHR_INDEX)
//        ,
//                getGlobalJobSpecificParameters(controlBam.executionContext.configuration))
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

    public static Tuple2<PhasedGenotypeFile, HaploblockGroupFile> imputeGenotypeX(TextFile sexFile, BamFile controlBam) {
        return (Tuple2<PhasedGenotypeFile, HaploblockGroupFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_IMPUTE_GENOTYPEX, controlBam, sexFile);
    }

    public static Tuple2<PhasedGenotypeFile, HaploblockGroupFile> imputeGenotypeX(TextFile sexFile, UnphasedGenotypeFileGroupByChromosome unphasedGenotypeFiles) {
        return (Tuple2<PhasedGenotypeFile, HaploblockGroupFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_IMPUTE_GENOTYPEX_NOMPILEUP, unphasedGenotypeFiles.getFiles().get("X"), sexFile);
    }

    public static TextFile getGenotypes(TextFile mergedAndFilteredSNPFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_GET_GENOTYPES, mergedAndFilteredSNPFile);
    }

    public static UnphasedGenotypeFileGroupByChromosome createUnphased(TextFile genotypeSNPFile) {
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


    public static TextFile addHaploTypes(TextFile mergedAndFilteredSNPFile, PhasedGenotypeFileGroupByChromosome phasedGenotypes, PhasedGenotypeFile phasedGenotypeX) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_ADD_HAPLOTYPES_TO_SNP_FILE, mergedAndFilteredSNPFile, phasedGenotypes, phasedGenotypeX);
    }

    public static TextFile createControlBafPlot(TextFile haplotypedSNPFile, TextFile genderFile) {
        return (TextFile) GenericMethod.callGenericTool(ACEseqConstants.TOOL_CREATE_CONTROL_BAF_PLOTS, haplotypedSNPFile, genderFile);
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
    static Tuple2<SVFile, TextFile> mergeSv(TextFile knownSegmentsFile, boolean runWithSv) {
        if (runWithSv) {
            BaseFile svFile = getSVFile(knownSegmentsFile)
            boolean b = FileSystemAccessProvider.getInstance().checkBaseFiles(svFile);
            if (b)
                return (Tuple2<SVFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MERGE_BREAKPOINTS_AND_SV, knownSegmentsFile, svFile);

            knownSegmentsFile.getExecutionContext().addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("SV files were not found in input path.", Level.WARNING));
            return null;
        } else {
            return (Tuple2<SVFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MERGE_BREAKPOINTS_WITHOUT_SV, knownSegmentsFile);
        }
    }


    static BaseFile getSVFile(BaseFile anyFile) {
        BaseFile svFile = BaseFile.deriveFrom(anyFile, SVFile.class.name)
        svFile.setAsSourceFile()
        return svFile
    }

    @ScriptCallingMethod
    public static Tuple2<SVFile, TextFile> mergeCrest(TextFile knownSegmentsFile) {
        TextFile svFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "crestDelDupInvFileTag", null, null);
        TextFile translocFile = (TextFile) BaseFile.constructManual(TextFile.class, knownSegmentsFile, null, null, null, null, "crestTranslocFileTag", null, null);

        boolean b = FileSystemAccessProvider.getInstance().checkBaseFiles(svFile, translocFile);
        if (!b) {
            knownSegmentsFile.getExecutionContext().addErrorEntry(ExecutionContextError.EXECUTION_NOINPUTDATA.expand("Crest files were not found in input path."));
            return null;
        }

        return (Tuple2<SVFile, TextFile>) GenericMethod.callGenericTool(ACEseqConstants.TOOL_MERGE_BREAKPOINTS_AND_SV_CREST, knownSegmentsFile, svFile, translocFile);
    }

    @ScriptCallingMethod
    public static TextFile getSegmentAndGetSnps(BaseFile breaks, TextFile pscbsSnps) {
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
