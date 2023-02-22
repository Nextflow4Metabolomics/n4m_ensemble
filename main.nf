#!/usr/bin/env nextflow

/**
    A Nextflow-based reproducible pipeline for untargeted metabolomics data analysis
    Description  : Post processing for results obtained from other repo in the Nextflow4Metabolomics Suite
    Author       : Xinsong Du
    Creation Date: 02/20/2022
    License      : MIT License
          
    This script is free software: you can redistribute it and/or modify
    it under the terms of the MIT License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    MIT License for more details.
    
    You should have received a copy of the MIT License
    along with this script. If not, see <https://opensource.org/licenses/MIT>.
    
    For any bugs or problems found, please contact us at
    - xinsongdu@ufl.edu (or xinsongdu@gmail.com), djlemas@ufl.edu
    - https://github.com/Nextflow4Metabolomics/n4m_results_summary
*/

// Those variable names which are all uppercase are channel names

version='1.0dev'
timestamp='20230220'

DATA_DIR = Channel.fromPath(params.input_dir, type: 'dir') // Location of folder storing positive data

// Config files
MSDIAL_CONFIG = Channel.fromPath(params.msdial_config)
MSFLO_CONFIG = Channel.fromPath(params.msflo_config)

// Library
MS1_LIBRARY = Channel.fromPath(params.ms1_library)
MS2_LIBRARY = Channel.fromPath(params.ms2_library)

// Reference files
REF = Channel.fromPath(params.ref)

/**
    Basic running information
*/

println "Project : $workflow.projectDir"
println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Manifest's pipeline version: $workflow.manifest.version"

/**
    Prints help when asked for
*/

if (params.help) {
    System.out.println("")
    System.out.println("A Nextflow-based reproducible pipeline for untargeted metabolomics data analysis - Version: $version ($timestamp)")
    System.out.println("This pipeline is distributed in the hope that it will be useful")
    System.out.println("but WITHOUT ANY WARRANTY. See the MIT License for more details.")
    System.out.println("")
    System.out.println("Please report comments and bugs to the issue tracking on GitHub")
    System.out.println("at https://github.com/Nextflow4Metabolomics/n4m_results_summaryissues.")
    System.out.println("Check https://github.com/Nextflow4Metabolomics/n4m_results_summary for updates, and refer to")
    System.out.println("https://github.com/Nextflow4Metabolomics/n4m_results_summary/wiki")
    System.out.println("")
    System.out.println("Usage:  ")
    System.out.println("   nextflow main.nf -profile [options: functional_test; docker; singularity]")
    System.out.println("")
    System.out.println("Arguments:")
    System.out.println("----------------------------- common parameters ----------------------------------")
    System.out.println("    -profile                                docker (run pipeline locally), singularity (run pipeline on supercomputer), or test (run test data locally with docker)")
    System.out.println("    --help                                  whether to show help information or not, default is null")
    System.out.println("Please refer to nextflow.config for more options.")
    System.out.println("")
    System.out.println("The workflow supports .csv format files.")
    System.out.println("The workflow supports MacOS and Linux operating systems.")
    System.out.println("")
    exit 1
}

custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Header log info
def summary = [:]
summary['Pipeline Name']  = 'N4M Summary'
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Input']            = params.input
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/** 
    Process description: Process for running MS-DIAL with negative mode batchfile and data to generate peak table of negative mode.
    Inputs: MS-DIAL config file; raw .mzML data; MS1 library; MS2 library; reference file.
    Outputs: The peak table produced by MS-DIAL after processing all raw data.
*/
process peak_merge {

    echo true

    publishDir './results/', mode: 'copy'

    input:
    file table_1 from TABLE_1_P // Location of the first peak table
    file table_2 from TABLE_2_P // location of the second peak table
    file peak_merge_code from PEAK_MERGE_CODE // Locatio of the python code for merging peaks
    val m from COLUMN_NAME_FOR_MZ_VALUES // 
    val r from COLUMN_NAME_FOR_RT_VALUES //
    val p from PPM_THRESHOLD_FOR_MERGING //
    val t from RT_THRESHOLD_FOR_MERGING // retention tiem threshold for merging
    val o from MERGING_RESULT // output file path and name

    output:
    file "*.csv" into TABLE_MERGED // merged peak table.

    shell:
    """
    echo "peak merging" &&
    python ${peak_merge_code} -i ${table_1} -j ${table_2} -m ${m} -r ${r} -p ${p} -t ${t} -o ${o}
    """
}

/** 
    Process description: MS-FLO download.
    inputs: N/A.
    outputs: The MS-FLO repo.
*/
process quantification {
    
    publishDir './results/', mode: 'copy'

    input:
    file table_1 from TABLE_1_Q // Location of the first peak table
    file table_2 from TABLE_2_Q // location of the second peak table
    file table_merged from TABLE_MERGED //
    file quantification_code from QUANTIFICATION_CODE
    val a from COLUMN_NAME_FOR_ANNOTATION_STATUS // column_name_for_annotation_status
    val d from DESIRED_IDENTIFICATION_STATUS // desired_identification_status
    val m from COLUMN_NAME_FOR_ANNOTATED_METABOLITES // column_name_for_annotated_metabolites
    val o from QUANTIFICATION_RESULT // output file location

    output:
    file "*.csv" into QUANTIFICATION

    shell:
    """
    echo "quantification" &&
    python ${quantification_code} -i ${table_1} -j ${table_2} -k ${table_merged} -a ${a} -d ${d} -m ${m} -o ${o}
    """

}