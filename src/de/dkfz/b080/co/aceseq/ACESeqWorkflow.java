/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.aceseq;

import de.dkfz.b080.co.common.*;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.config.*;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.files.Tuple3;
import de.dkfz.roddy.knowledge.files.GenericFileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

/**
 */
public class ACESeqWorkflow extends WorkflowUsingMergedBams {


    private boolean allowMissingSVFile(ExecutionContext context) {
        return getflag(context, "allowMissingSVFile", false);
    }

    private boolean runWithSV(ExecutionContext context) {
        return getflag(context, "runWithSv", true);
    }

    @Override
    public boolean execute(ExecutionContext context, BasicBamFile _bamControlMerged, BasicBamFile _bamTumorMerged) {


        boolean runWithCrest = getFlag("runWithCrest", false);
        boolean runQualityCheckOnly = getFlag("runQualityCheckOnly", false);
        boolean runWithFakeControl = getFlag("runWithFakeControl", false);
        boolean runWithoutControl = getFlag("runWithoutControl", false);

        BamFile bamControlMerged = new BamFile(_bamControlMerged);
        BamFile bamTumorMerged = new BamFile(_bamTumorMerged);
        context.getConfigurationValues().add(new ConfigurationValue("tumorSample", ((COFileStageSettings) _bamTumorMerged.getFileStage()).getSample().getName()));
        context.getConfigurationValues().add(new ConfigurationValue("controlSample", ((COFileStageSettings) _bamControlMerged.getFileStage()).getSample().getName()));

        CnvSnpGeneratorResultByType resultByType;
        resultByType = ACESeqMethods.generateCNVSNPs(bamControlMerged, bamTumorMerged);

        //TODO The annotate job tool id is not visible in the jobstate logfile.
        CoverageWindowsFileAnnotationResult annotationResult = resultByType.getCoverageWindowsFiles().annotate();
        TextFile replaceControlFile;
        TextFile mergedAndFilteredCoverageWindowFiles;
        if (runWithFakeControl || runWithoutControl) {
            replaceControlFile = ACESeqMethods.replaceControl(annotationResult.getGenderFile());
            mergedAndFilteredCoverageWindowFiles = GenericMethod.callGenericTool("mergeAndFilterCnvFiles_withReplaceBadControl", replaceControlFile, new GenericFileGroup(annotationResult.getListOfFiles()));
        } else {
            mergedAndFilteredCoverageWindowFiles = annotationResult.mergeAndFilterCoverageWindowFiles();
        }
        Tuple3<TextFile, TextFile, TextFile> correctedWindowFile = ACESeqMethods.correctGC(mergedAndFilteredCoverageWindowFiles);

        if (runQualityCheckOnly)
            return true;

        ImputeGenotypeByChromosome imputedGenotypeByChromosome;
        Tuple2<PhasedGenotypeFile, HaploblockGroupFile> phasedGenotypeX;
        TextFile haplotypedSNPFile;

        if (runWithoutControl) {
            TextFile mergedAndFilteredSNPFile = resultByType.getPositionFiles().mergeAndFilter();
            TextFile genotypeSNPFile = ACESeqMethods.getGenotypes(mergedAndFilteredSNPFile);
            UnphasedGenotypeFileGroupByChromosome unphasedGenotypeFile = ACESeqMethods.createUnphased(genotypeSNPFile);
            imputedGenotypeByChromosome = ACESeqMethods.imputeGenotypes(unphasedGenotypeFile);
            phasedGenotypeX = ACESeqMethods.imputeGenotypeX(annotationResult.getGenderFile(), unphasedGenotypeFile);
            haplotypedSNPFile = ACESeqMethods.addHaploTypes(genotypeSNPFile, imputedGenotypeByChromosome.getPhasedSnpFiles(), phasedGenotypeX.value0);

        } else {

            TextFile mergedAndFilteredSNPFile = resultByType.getPositionFiles().mergeAndFilter();
            imputedGenotypeByChromosome = ACESeqMethods.imputeGenotypes(bamControlMerged);
            phasedGenotypeX = ACESeqMethods.imputeGenotypeX(annotationResult.getGenderFile(), bamControlMerged);
            haplotypedSNPFile = ACESeqMethods.addHaploTypes(mergedAndFilteredSNPFile, imputedGenotypeByChromosome.getPhasedSnpFiles(), phasedGenotypeX.value0);
            TextFile baffile = ACESeqMethods.createControlBafPlot(haplotypedSNPFile, annotationResult.getGenderFile());
        }


        Tuple2<TextFile, TextFile> breakpoints = ACESeqMethods.pscbsGaps(haplotypedSNPFile, correctedWindowFile.value0, annotationResult.getGenderFile());
        Tuple2<BreakpointsFile, TextFile> mergedSvs = null;


        if (runWithSV(context)) {
            mergedSvs = ACESeqMethods.mergeSv(breakpoints.value0, _bamTumorMerged);
            if (mergedSvs == null) {
                return allowMissingSVFile(context); // Here, exit with error (false) is possible
            }
        } else if (runWithCrest) {
            mergedSvs = ACESeqMethods.mergeCrest(breakpoints.value0);
            if (mergedSvs == null)
                return false;  // Getting no merged SVs from the Crest step is always wrong.
        } else {
            ACESeqMethods.mergeNoSv(breakpoints.value0);
            return true;
        }

        TextFile pscbsSegments = ACESeqMethods.getSegmentAndGetSnps(mergedSvs.value0, breakpoints.value1);
        TextFile homoDelSegments = ACESeqMethods.markSegsWithHomozygDel(pscbsSegments, mergedSvs.value1);
        TextFile homoDelSnps = ACESeqMethods.segsToSnpDataHomodel(homoDelSegments, breakpoints.value1);
        Tuple2<TextFile, TextFile> clusteredSegments = ACESeqMethods.clusterPruneSegments(homoDelSegments, homoDelSnps, annotationResult.getGenderFile(), correctedWindowFile.value1);
        TextFile clusteredSnps = ACESeqMethods.segsToSnpDataPruned(clusteredSegments.value0, clusteredSegments.value1);
        TextFile peakSegments = ACESeqMethods.estimatePeaks(clusteredSegments.value0, clusteredSnps, annotationResult.getGenderFile());
        TextFile purityPloidy = ACESeqMethods.estimatePurityPloidy(peakSegments, annotationResult.getGenderFile());
        Tuple2<TextFile, TextFile> results = ACESeqMethods.generatePlots(peakSegments, clusteredSnps, mergedSvs.value1, purityPloidy, annotationResult.getGenderFile());
        TextFile finalVcf = ACESeqMethods.estimateHRD(annotationResult.getGenderFile(), results.value1);

        return true;
    }

    private boolean checkSvFileIfNecessary(ExecutionContext context) {
        // Just return true, if the file is not used.
        if (!runWithSV(context))
            return true;


        BasicBamFile bamControlMerged = new BasicBamFile(loadInitialBamFilesForDataset(context)[0]);
        BasicBamFile bamTumorMerged = new BasicBamFile(loadInitialBamFilesForDataset(context)[1]);
        context.getConfigurationValues().add(new ConfigurationValue("tumorSample", ((COFileStageSettings) bamTumorMerged.getFileStage()).getSample().getName()));
        context.getConfigurationValues().add(new ConfigurationValue("controlSample", ((COFileStageSettings) bamControlMerged.getFileStage()).getSample().getName()));


        // Check, if the file exists
        boolean fileIsAccessible = context.fileIsAccessible(ACESeqMethods.getSVFile(bamTumorMerged).getPath());

        // It is either allowed to run with the file or we see, if the file really exists.
        return allowMissingSVFile(context) || fileIsAccessible;
    }

//    @Override
//    public boolean setupExecution(ExecutionContext context) {
//        return super.setupExecution(context);
//    }

    @Override
    public boolean checkExecutability(ExecutionContext context) {
        boolean executable = super.checkExecutability(context);
        executable &= checkSvFileIfNecessary(context);
        return executable;
    }
}
