/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.files;

import de.dkfz.b080.co.aceseq.ACEseqConstants;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.execution.jobs.BEJobResult;
import de.dkfz.roddy.execution.jobs.Job;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.files.GenericFileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by michael on 11.06.14.
 */
public class CoverageWindowsFileGroupByChromosome extends FileGroup {

    private Map<String, CoverageWindowsFile> files;

    public CoverageWindowsFileGroupByChromosome(Map<String, CoverageWindowsFile> files) {
        super(new LinkedList<>(files.values()));
        this.files = files;
    }

    public Map<String, CoverageWindowsFile> getFiles() {
        return files;
    }

    public CoverageWindowsFileAnnotationResult annotate() {
        GenderFile genderFile = (GenderFile)BaseFile.constructManual(GenderFile.class, files.get("1"), null, null, null, null, "genderFile", null, null);
        GenericFileGroup outFiles = GenericMethod.callGenericToolWithFileGroupOutput(ACEseqConstants.TOOL_ANNOTATE_COV_WIN,  genderFile, new LinkedList<>(files.keySet()), this);

        return new CoverageWindowsFileAnnotationResult(outFiles.getFilesInGroup(), genderFile);
    }
}
