/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.files.GenericFileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.List;

import static de.dkfz.b080.co.aceseq.ACEseqConstants.TOOL_MERGE_AND_FILTER_CNV_FILES;

/**
 * Created by kleinhei on 6/16/14.
 */
public class CoverageWindowsFileAnnotationResult extends FileGroup {
    private final List<TextFile> listOfFiles;
    private final GenderFile genderFile;

    public CoverageWindowsFileAnnotationResult(List<TextFile> listOfFiles, GenderFile genderFile) {
        super(listOfFiles);
        super.addFile(genderFile);
        this.listOfFiles = listOfFiles;
        this.genderFile = genderFile;
    }

    public List<TextFile> getListOfFiles() {
        return listOfFiles;
    }

    public TextFile mergeAndFilterCoverageWindowFiles() {
        TextFile file = GenericMethod.callGenericTool(TOOL_MERGE_AND_FILTER_CNV_FILES, genderFile, new GenericFileGroup(listOfFiles));
        return file;
     }

    public GenderFile getGenderFile() {
        return genderFile;
    }
}
