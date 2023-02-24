#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description : This code quantifies obtained peak tables.
Author      : Xinsong Du
License     : MIT License
Maintainer  : xinsongdu@ufl.edu
Usage       : python quantification.py -i $first_input_peak_table
                                       -j $second_input_peak_table
                                       -k $merged_peak_table
                                       -a $column_name_for_annotation_status
                                       -d $desired_identification_status
                                       -m $column_name_for_annotated_metabolites
                                       -o $output_quantification_table
"""
import warnings
import logging
import logging.handlers
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s]: %(levelname)s: %(message)s')

warnings.filterwarnings('ignore')

def quantification(t_1, t_2, t_merged, col_name_annotation_status, desired_identification_status, col_name_metabolite_name):
    """Add threshold for blank subtraction algorithm.

    # Arguments:
        t_1                           : the first peak table.
        t_2                           : the second peak table.
        t_merged                      : the merged peak table.
        col_name_identification_status: the name of the column including information about how each metabolite was annotated.
        desired_identification_status : the value indicating the desired annotation method.
        col_name_metabolite_name      : the column name of annotated metabolites (true positive peaks).

    # Returns:
        A dataframe showing the quantification result table.
    """

    n_rows_t_1 = t_1.shape[0]
    n_rows_t_2 = t_2.shape[0]
    n_rows_t_merged = t_merged.shape[0]
    
    n_annotated_rows_t_1 = t_1[t_1[col_name_annotation_status] == desired_identification_status].shape[0]
    n_annotated_rows_t_2 = t_2[t_2[col_name_annotation_status] == desired_identification_status].shape[0]
    n_annotated_rows_t_merged = t_merged[t_merged[col_name_annotation_status] == desired_identification_status].shape[0]
    
    tpr_t_1 = round(n_annotated_rows_t_1/n_rows_t_1, 4)
    tpr_t_2 = round(n_annotated_rows_t_2/n_rows_t_2, 4)
    tpr_t_merged = round(n_annotated_rows_t_merged/n_rows_t_merged, 4)
    
    n_deduplicated_metabolites_t_1 = t_1[t_1[col_name_annotation_status] == desired_identification_status].drop_duplicates(subset = col_name_metabolite_name).shape[0]
    n_deduplicated_metabolites_t_2 = t_2[t_2[col_name_annotation_status] == desired_identification_status].drop_duplicates(subset = col_name_metabolite_name).shape[0]
    n_deduplicated_metabolites_t_merged = t_merged[t_merged[col_name_annotation_status] == desired_identification_status].drop_duplicates(subset = col_name_metabolite_name).shape[0]
    
    d = {
        'Peak Tables': ['Peak Table 1', 'Peak Table 2', 'Merged Peak Table'],
        'Total Number of Peaks': [n_rows_t_1, n_rows_t_2, n_rows_t_merged],
        'Number of True Peaks': [n_annotated_rows_t_1, n_annotated_rows_t_2, n_annotated_rows_t_merged],
        'True Positive Rate': [tpr_t_1, tpr_t_2, tpr_t_merged],
        'Number of Distinct Metabolites': [n_deduplicated_metabolites_t_1, n_deduplicated_metabolites_t_2, n_deduplicated_metabolites_t_merged]
    }
    
    df = pd.DataFrame(data=d)
    
    return df

if __name__ == '__main__':

    logger.info('generating quantification results')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input_1', help="the location of the first peak table;", \
        default="table_1.csv", dest="input_1", required=False)
    parser.add_argument(
        '-j', '--input_2', help="the location of the second peak table;", \
        default="pos_design.csv", dest="input_2", required=False)
    parser.add_argument(
        '-k', '--input_merged', help="the location of the merged peak table;", \
        default="pos_withstats.csv", dest="input_merged", required=False)
    parser.add_argument(
        '-a', '--col_name_annotation_status', help="the column name for the metabolite annotation method;", \
        default="Annotation tag (VS1.0)", dest="col_name_annotation_status", required=False)
    parser.add_argument(
        '-d', '--desired_identification_status', help="the desired identification method;", \
        default="annotated by user-defined text library", dest="desired_identification_status", required=False)
    parser.add_argument(
        '-m', '--col_name_metabolite_name', help="the column name for annotated metabolites;", \
        default="Metabolite name", dest="col_name_metabolite_name", required=False)
    parser.add_argument(
        '-o', '--output', help="the location of the output summary table;", \
        default="quantification.csv", required=False)

    args = parser.parse_args()

    table_1 = pd.read_csv(args.input_1, delimiter="\t")
    table_2 = pd.read_csv(args.input_2, delimiter="\t")
    table_merged = pd.read_csv(args.input_merged)
    summary_df = quantification(table_1, table_2, table_merged, args.col_name_annotation_status, \
                                args.desired_identification_status, args.col_name_metabolite_name)

    summary_df.to_csv(args.output, index=False)
