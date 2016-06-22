package de.dkfz.b080.co.files;

import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.files.FileObject;
import de.dkfz.roddy.knowledge.files.IndexedFileObjects;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by michael on 11.06.14.
 */
public class CnvSnpGeneratorResultByType extends IndexedFileObjects {
    public CnvSnpGeneratorResultByType(IndexedFileObjects indexedFileObjects, ExecutionContext executionContext) {
        super(indexedFileObjects.getIndices(), indexedFileObjects.getIndexedFileObjects(), executionContext);
    }

    public SNPPositionFileGroupByChromosome getPositionFiles() {
        Map<String, SNPPositionFile> files = new LinkedHashMap<>();
        Map<String, FileObject> indexedFileObjects = super.getIndexedFileObjects();
        for (String chr : (List<String>)indices) {
            CnvSnpGeneratorResultFileGroup fg = (CnvSnpGeneratorResultFileGroup) indexedFileObjects.get(chr);
            files.put(chr, fg.getPositionFile());
        }
        return new SNPPositionFileGroupByChromosome(files);
    }

    public CoverageWindowsFileGroupByChromosome getCoverageWindowsFiles() {
        Map<String, CoverageWindowsFile> files = new LinkedHashMap<>();
        Map<String, FileObject> indexedFileObjects = super.getIndexedFileObjects();
        for (String chr : (List<String>)indices) {
            CnvSnpGeneratorResultFileGroup fg = (CnvSnpGeneratorResultFileGroup) indexedFileObjects.get(chr);
            files.put(chr, fg.getWindowsFile());
        }
        return new CoverageWindowsFileGroupByChromosome(files);
    }
}
