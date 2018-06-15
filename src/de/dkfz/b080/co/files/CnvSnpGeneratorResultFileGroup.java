/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.List;

/**
 * Created by michael on 11.06.14.
 */
public class CnvSnpGeneratorResultFileGroup extends FileGroup {
    private SNPPositionFile positionFile;
    private CoverageWindowsFile windowsFile;

    public CnvSnpGeneratorResultFileGroup(List<BaseFile> files) {
        super(files);
        positionFile = (SNPPositionFile) files.get(0);
        windowsFile = (CoverageWindowsFile) files.get(1);
    }

    public SNPPositionFile getPositionFile() {
        return positionFile;
    }

    public CoverageWindowsFile getWindowsFile() {
        return windowsFile;
    }
}
