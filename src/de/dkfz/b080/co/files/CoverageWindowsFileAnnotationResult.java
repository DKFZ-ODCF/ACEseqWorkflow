package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.files.GenericFileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.LinkedList;
import java.util.List;

import static de.dkfz.b080.co.files.COConstants.TOOL_MERGE_AND_FILTER_CNV_FILES;

/**
 * Created by kleinhei on 6/16/14.
 */
public class CoverageWindowsFileAnnotationResult extends FileGroup {
    private final List<TextFile> listOfFiles;
    private final TextFile genderFile;

    public CoverageWindowsFileAnnotationResult(List<TextFile> listOfFiles, TextFile genderFile) {
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

    public TextFile getGenderFile() {
        return genderFile;
    }
}
