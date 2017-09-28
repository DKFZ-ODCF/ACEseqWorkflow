/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 */
public class CoverageWindowsFile extends COBaseFile {
    public CoverageWindowsFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public CoverageWindowsFile(BamFile parentFile) {
        super(parentFile);
    }
}
