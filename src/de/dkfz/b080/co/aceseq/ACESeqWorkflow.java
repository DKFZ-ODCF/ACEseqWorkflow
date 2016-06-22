package de.dkfz.b080.co.aceseq;

import de.dkfz.b080.co.common.WorkflowUsingMergedBams;
import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.core.ExecutionContext;
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

        boolean runWithDelly = context.getConfiguration().getConfigurationValues().getBoolean("runWithDelly", true);
        boolean runQualityCheckOnly = context.getConfiguration().getConfigurationValues().getBoolean("runQualityCheckOnly", false);
        boolean runMetaCNVGeneration = context.getConfiguration().getConfigurationValues().getBoolean("runMetaCNVGeneration", false);
        boolean runWithFakeControl = context.getConfiguration().getConfigurationValues().getBoolean("runWithFakeControl", false);

        CnvSnpGeneratorResultByType resultByType = null;
        if (runMetaCNVGeneration) {
            resultByType = ACESeqMethods.generateCNVSNPsMeta(bamControlMerged, bamTumorMerged);
            return true;
        } else {
            resultByType = ACESeqMethods.generateCNVSNPs(bamControlMerged, bamTumorMerged);
        }

        //TODO The annotate job tool id is not visible in the jobstate logfile.
        CoverageWindowsFileAnnotationResult annotationResult = resultByType.getCoverageWindowsFiles().annotate();
        TextFile replaceControlFile = null;
        TextFile mergedAndFilteredCoverageWindowFiles = null;
        if (runWithFakeControl) {
            replaceControlFile = ACESeqMethods.replaceControl(annotationResult.getGenderFile());
            mergedAndFilteredCoverageWindowFiles = GenericMethod.callGenericTool("mergeAndFilterCnvFiles_withReplaceBadControl", replaceControlFile, new GenericFileGroup(annotationResult.getListOfFiles()));
        } else {
            mergedAndFilteredCoverageWindowFiles = annotationResult.mergeAndFilterCoverageWindowFiles();
        }
        Tuple3<TextFile, TextFile, TextFile> correctedWindowFile = ACESeqMethods.correctGC(mergedAndFilteredCoverageWindowFiles);

        if (runQualityCheckOnly)
            return true;

        ImputeGenotypeByChromosome imputedGenotypeByChromosome = ACESeqMethods.imputeGenotypes(bamControlMerged);
        TextFile mergedAndFilteredSNPFile = resultByType.getPositionFiles().mergeAndFilter();
        Tuple2<PhasedGenotypeFile, HaploblockGroupFile> phasedGenotypeX = ACESeqMethods.imputeGenotypeX(annotationResult.getGenderFile(), bamControlMerged);
        TextFile haplotypedSNPFile = ACESeqMethods.addHaploTypes(mergedAndFilteredSNPFile, imputedGenotypeByChromosome.getPhasedSnpFiles(), phasedGenotypeX.value0);
        Tuple2<TextFile, TextFile> breakpoints = ACESeqMethods.pscbsGaps(haplotypedSNPFile, correctedWindowFile.value0, annotationResult.getGenderFile());
        Tuple2<TextFile, TextFile> mergedSvs = null;

        if (runWithDelly)
            mergedSvs = ACESeqMethods.mergeDelly(breakpoints.value0);
        else
            mergedSvs = ACESeqMethods.mergeCrest(breakpoints.value0);
        if (mergedSvs == null)
            return true;

        TextFile pscbsSegments = ACESeqMethods.getSegmentAndGetSnps(mergedSvs.value0, breakpoints.value1);
        TextFile homoDelSegments = ACESeqMethods.markSegsWithHomozygDel(pscbsSegments, mergedSvs.value1);
        TextFile homoDelSnps = ACESeqMethods.segsToSnpDataHomodel(homoDelSegments, breakpoints.value1);
        Tuple2<TextFile, TextFile> clusteredSegments = ACESeqMethods.clusterPruneSegments(homoDelSegments, homoDelSnps, annotationResult.getGenderFile(), correctedWindowFile.value1, imputedGenotypeByChromosome.getHaploblockFiles(), phasedGenotypeX.value1);
        TextFile clusteredSnps = ACESeqMethods.segsToSnpDataPruned(clusteredSegments.value0, clusteredSegments.value1);
        TextFile peakSegments = ACESeqMethods.estimatePeaks(clusteredSegments.value0, clusteredSnps, annotationResult.getGenderFile());
        TextFile purityPloidy = ACESeqMethods.estimatePurityPloidy(peakSegments, annotationResult.getGenderFile());
        TextFile results = ACESeqMethods.generatePlots(peakSegments, clusteredSnps, mergedSvs.value1, purityPloidy, annotationResult.getGenderFile());
        TextFile finalVcf = ACESeqMethods.convertToVcf(results);

        return true;
    }
}
