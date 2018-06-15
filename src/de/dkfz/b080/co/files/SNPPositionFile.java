/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.files;

/**
 */
public class SNPPositionFile extends COBaseFile {
    public SNPPositionFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public SNPPositionFile(BamFile parentFile) {
        super(parentFile);
    }
}
