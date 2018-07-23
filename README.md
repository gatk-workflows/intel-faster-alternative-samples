# intel-faster-alternative-samples
Miscellaneous workflows optimized by Intel to be fast-running.

### WORKFLOWS AND JSONS
This repository {currently} contains a faster version of the Germline Single Sample Calling variant workflow. 

├── faster\_alternative\_workflow\_with\_bwa.wdl *&rarr;* Workflow enabled with Bwa, Samblaster, Samtools, HaplotypeCaller \
├── faster\_alternative\_workflow\_with\_minimap2.wdl *&rarr;* Workflow enabled with Minimap2, Samblaster, Samtools, HaplotypeCaller \
├── faster\_alternative\_workflow.inputs.json *&rarr;* JSON file to support above WDL workflows

In the JSON file, please follow the comments to make modifications to certail requirements. Also modify the paths to the datasets and tools where they reside in your cluster. \

#### FPGA CHANGES
Assuming the environemnt has been setup to offload the pairhmm kernel of HaplotypeCaller to FPGA - the below changes must be enabled in the WDL/JSON files (based on the comments) to make use of the FPGA. 

a. In the WDL files, for task Haplotype Caller runtime section, uncomment the line:
require\_fpga: "yes"

b. In the JSON file, change the "PairedEndSingleSampleWorkflow.gatk\_gkl\_pairhmm\_implementation" to “EXPERIMENTAL\_FPGA\_LOGLESS\_CACHING” from “AVX\_LOGLESS\_CACHING”.

### DATASETS
The sample JSON file uses the HG02769 fastq datasets which can be obtained from the Broad Institute. 
This workflow can also be run with other fastq files. 

### TOOLS
For on-prem, the workflow uses non-dockerized tools. The latest version of the tools used in the WDL workflow can be downloaded from their respective websites. 
Once the tools are compiled and ready to use, please ensure https://github.com/gatk-workflows/intel-faster-alternative-samples/blob/master/faster\_alternative\_workflow.inputs.json#L47 is modified to reflect the same.
