package de.dkfz.b080.co.files;

import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.execution.jobs.Job;
import de.dkfz.roddy.execution.jobs.JobResult;
import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;

import java.io.File;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by michael on 11.06.14.
 */
public class HaploblockFileGroupByChromosome extends FileGroup {

    private Map<String, HaploblockGroupFile> files;

    public HaploblockFileGroupByChromosome(Map<String, HaploblockGroupFile> files) {
        super(new LinkedList<>(files.values()));
        this.files = files;
    }

    public Map<String, HaploblockGroupFile> getFiles() {
        return files;
    }

}
