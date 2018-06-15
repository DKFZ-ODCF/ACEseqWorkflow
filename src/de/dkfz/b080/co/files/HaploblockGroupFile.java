/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;

/**
 * Created by kleinhei on 6/13/14.
 */
public class HaploblockGroupFile  extends COBaseFile {
    public HaploblockGroupFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }

    public HaploblockGroupFile(BamFile parentFile) {
        super(parentFile);
    }
}
