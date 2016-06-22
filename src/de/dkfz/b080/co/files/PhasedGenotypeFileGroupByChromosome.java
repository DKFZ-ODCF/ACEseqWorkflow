package de.dkfz.b080.co.files;

import de.dkfz.roddy.execution.jobs.ScriptCallingMethod;
import de.dkfz.roddy.knowledge.files.FileGroup;
import de.dkfz.roddy.knowledge.methods.GenericMethod;

import java.util.LinkedList;
import java.util.Map;
/**
 * Created by kleinhei on 6/16/14.
 */
public class PhasedGenotypeFileGroupByChromosome extends FileGroup<PhasedGenotypeFile> {

        private Map<String, PhasedGenotypeFile> files;

        public PhasedGenotypeFileGroupByChromosome(Map<String, PhasedGenotypeFile> files) {
            super(new LinkedList<>(files.values()));
            this.files = files;
        }

        public Map<String, PhasedGenotypeFile> getFiles() {
            return files;
        }

//       @ScriptCallingMethod
//       public TextFile mergeAndFilter() {
//            BamFile bf = (BamFile)(files.get("1")).getParentFiles().get(0); //Should be merged tumor bam file//
//            TextFile file = GenericMethod.callGenericTool(COConstants.TOOL_MERGE_AND_FILTER_SNP_FILES, bf, this);
//            return file;
//        }
//    }
}
