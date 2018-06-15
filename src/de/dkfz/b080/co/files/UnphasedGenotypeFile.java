/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 */
public class UnphasedGenotypeFile extends COBaseFile {
    public UnphasedGenotypeFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public UnphasedGenotypeFile(BamFile parentFile) {
        super(parentFile);
    }
}
