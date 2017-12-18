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
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.files.Tuple3;
import de.dkfz.roddy.knowledge.files.GenericFileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

/**
 */
public class ACESeqWorkflow extends WorkflowUsingMergedBams {

    @Override
    public boolean execute(ExecutionContext context, BasicBamFile _bamControlMerged, BasicBamFile _bamTumorMerged) {

        BamFile bamControlMerged = new BamFile(_bamControlMerged);
        BamFile bamTumorMerged = new BamFile(_bamTumorMerged);

        boolean runWithSv = context.getConfiguration().getConfigurationValues().getBoolean("runWithSv", true);
        boolean runWithCrest = context.getConfiguration().getConfigurationValues().getBoolean("runWithCrest", false);
        boolean runQualityCheckOnly = context.getConfiguration().getConfigurationValues().getBoolean("runQualityCheckOnly", false);
        boolean runWithFakeControl = context.getConfiguration().getConfigurationValues().getBoolean("runWithFakeControl", false);
        boolean runWithoutControl = context.getConfiguration().getConfigurationValues().getBoolean("runWithoutControl", false);

        context.getConfiguration().getConfigurationValues().add(new ConfigurationValue("tumorSample", ((COFileStageSettings) _bamTumorMerged.getFileStage()).getSample().getName()));
        context.getConfiguration().getConfigurationValues().add(new ConfigurationValue("controlSample", ((COFileStageSettings) _bamControlMerged.getFileStage()).getSample().getName()));

        CnvSnpGeneratorResultByType resultByType = null;
        resultByType = ACESeqMethods.generateCNVSNPs(bamControlMerged, bamTumorMerged);

        //TODO The annotate job tool id is not visible in the jobstate logfile.
        CoverageWindowsFileAnnotationResult annotationResult = resultByType.getCoverageWindowsFiles().annotate();
        TextFile replaceControlFile = null;
        TextFile mergedAndFilteredCoverageWindowFiles = null;
        if (runWithFakeControl || runWithoutControl) {
            replaceControlFile = ACESeqMethods.replaceControl(annotationResult.getGenderFile());
            mergedAndFilteredCoverageWindowFiles = GenericMethod.callGenericTool("mergeAndFilterCnvFiles_withReplaceBadControl", replaceControlFile, new GenericFileGroup(annotationResult.getListOfFiles()));
        } else {
            mergedAndFilteredCoverageWindowFiles = annotationResult.mergeAndFilterCoverageWindowFiles();
        }
        Tuple3<TextFile, TextFile, TextFile> correctedWindowFile = ACESeqMethods.correctGC(mergedAndFilteredCoverageWindowFiles);

        if (runQualityCheckOnly)
            return true;

        ImputeGenotypeByChromosome imputedGenotypeByChromosome = null;
        Tuple2<PhasedGenotypeFile, HaploblockGroupFile> phasedGenotypeX = null;
        TextFile haplotypedSNPFile = null;

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
        Tuple2<TextFile, TextFile> mergedSvs = null;

        if (runWithSv) {
            mergedSvs = ACESeqMethods.mergeSv(breakpoints.value0, runWithSv); // true is passed
        } else if (runWithCrest) {
            mergedSvs = ACESeqMethods.mergeCrest(breakpoints.value0);
        } else {
            ACESeqMethods.mergeSv(breakpoints.value0, runWithSv); // false is passed
            return true;
        }

        if (mergedSvs == null)
            return false;

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

    private boolean checkSvFile(ExecutionContext executionContext) {
        // TODO getInitialBamFiles(executionContext)[0]
        return true;
    }

    @Override
    public boolean checkExecutability(ExecutionContext context) {
        boolean executable = super.checkExecutability(context);
        executable &= checkSvFile(context);
        return executable;
    }
}
