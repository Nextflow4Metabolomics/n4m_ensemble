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

TABLE_1 = Channel.fromPath(params.table_1)
TABLE_2 = Channel.fromPath(params.table_2)
TABLE_1.into{TABLE_1_P; TABLE_1_Q}
TABLE_2.into{TABLE_2_P; TABLE_2_Q}

PEAK_MERGE_CODE = Channel.fromPath(params.peak_merge_code)
QUANTIFICATION_CODE = Channel.fromPath(params.quantification_code)

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

println "After help"

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
summary['Config Files'] = workflow.configFiles.join(', ')
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/** 
    Process description: Merging two peak tables.
    Inputs: Two peak tables that will be merges. Note that "table_2" should have less rows than "table_1".
    Outputs: The merged peak table.
*/
process peak_merge {

    echo true

    println "enter peak merge process"

    publishDir './results/', mode: 'copy'

    input:
    file table_1 from TABLE_1_P // Location of the first peak table
    file table_2 from TABLE_2_P // location of the second peak table
    file peak_merge_code from PEAK_MERGE_CODE // Locatio of the python code for merging peaks

    output:
    file "*.csv" into TABLE_MERGED // merged peak table

    shell:
    """
    echo "peak merging" &&
    python3 ${peak_merge_code} -i ${table_1} -j ${table_2} -m ${params.column_name_for_mz_values} -r ${params.column_name_for_rt_values} -p ${params.ppm_threshold_for_merging} -t ${params.rt_threshold_for_merging} -o ${params.merging_result}
    """
}

/** 
    Process description: Quantify peak detection and merging results.
    inputs: Two peak tables and the merged peak table.
    outputs: The summary table.
*/
process quantification {
    
    publishDir './results/', mode: 'copy'

    println "enter quantification process"

    input:
    file table_1 from TABLE_1_Q // Location of the first peak table
    file table_2 from TABLE_2_Q // location of the second peak table
    file table_merged from TABLE_MERGED //
    file quantification_code from QUANTIFICATION_CODE

    output:
    file "*.csv" into QUANTIFICATION

    shell:
    """
    echo "quantification" &&
    python3 ${quantification_code} -i ${table_1} -j ${table_2} -k ${table_merged} -a ${params.column_name_for_annotation_status} -d ${params.desired_identification_status} -m ${params.column_name_for_annotated_metabolites} -o ${params.quantification_result}
    """

}