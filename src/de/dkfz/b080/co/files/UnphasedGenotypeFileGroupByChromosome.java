/*
 * Copyright (c) 2017 The ACEseq workflow developers.
 *
 * Distributed under the MIT License (license terms are at https://www.github.com/eilslabs/ACEseqWorkflow/LICENSE.txt).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.IndexedFileObjects;

import java.util.List;
import java.util.Map;

/**
 * Created by michael on 11.06.14.
 */
public class UnphasedGenotypeFileGroupByChromosome extends IndexedFileObjects {

    private Map<String, UnphasedGenotypeFile> files;

    public UnphasedGenotypeFileGroupByChromosome(List<String> keyset, Map<String, UnphasedGenotypeFile> files, ExecutionContext context) {
        super(keyset, files, context);
        this.files = files;
    }

    public Map<String, UnphasedGenotypeFile> getFiles() {
        return files;
    }
}
