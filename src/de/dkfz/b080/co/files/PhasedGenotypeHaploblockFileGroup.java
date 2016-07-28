package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.List;

/**
 * Created by kleinhei on 6/13/14.
 */

public class PhasedGenotypeHaploblockFileGroup extends FileGroup {

    private PhasedGenotypeFile phasedSnpFile;
    private HaploblockGroupFile haploblockFile;

    public PhasedGenotypeHaploblockFileGroup(List<BaseFile> files) {
        super(files);
        phasedSnpFile = (PhasedGenotypeFile) files.get(0);
        haploblockFile = (HaploblockGroupFile) files.get(1);
    }

    public PhasedGenotypeFile getPhasedSnpFile() {
        return phasedSnpFile;
    }

    public HaploblockGroupFile getHaploblockFile() {
        return haploblockFile;
    }
}