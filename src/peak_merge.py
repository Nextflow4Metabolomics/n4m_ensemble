#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description : This code merged peak tables produced by different platforms.
Author      : Xinsong Du
License     : MIT License
Maintainer  : xinsongdu@ufl.edu
Usage       : python peak_merge.py -i $first_input_peak_table
                                   -j $second_input_peak_table
                                   -m $column_name_for_mz_values
                                   -r $column_name_for_rt_values
                                   -p $ppm_threshold_for_merging
                                   -t $rt_threshold_for_merging
                                   -o $output_file_location
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

def ppm_calculation(x, y):
    """calculating ppm value based on two mz values.

    # Arguments:
        x            : the measured mz value.
        y            : the exact mz value.

    # Returns:
        The merged peak table.
    """

    return abs(x-y)*10e6/y

def nearest_value(lst, K):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

def MS1_match(mz, rt, mz_list, rt_list, col_name_mz, col_name_rt, ppm_threshold, rt_tolerance):
#    print("enter match")
    df = pd.DataFrame(list(zip(mz_list, rt_list)),
               columns =[col_name_mz, col_name_rt])
    if len(df[(ppm_calculation(mz, df[col_name_mz])<=ppm_threshold) & (abs(df[col_name_rt] - rt) < rt_tolerance)])>0:
        return True
    else:
        return False

def overlapping_peaks(t_1, t_2, col_name_mz, col_name_rt, ppm_threshold, rt_tolerance):
    """Merging peaks from two peak tables.

    # Arguments:
        t_1            : the first peak table.
        t_2            : the second peak table, which should have less peaks than the first one.
        col_name_mz    : the column name of the column including the mass to charge ratio.
        col_name_rt    : the column name of the column including the retention time.
        ppm_threshold  : the threshold of part per million used for merging peaks.
        rt_tolerance   : the threshold of retention time tolerance for merging peaks.

    # Returns:
        The merged peak table.
    """

    flag = 0
    flagged_rows = []
    total = t_2.shape[0]
    if type(t_1[col_name_mz].iloc[0]) == np.float64:
        no_string = 1
    else:
        no_string = 0      
    mz_list = []
    rt_list = []
    if no_string == 1:
        mz_list = list(t_1[col_name_mz])
        rt_list = list(t_1[col_name_rt])
    else:
        for mz_element in t_1[col_name_mz].str.split("_"):
            mz_list += [float(x) for x in mz_element]
        for rt_element in t_1[col_name_rt].str.split("_"):
            rt_list += [float(x) for x in rt_element]
    for i in range(total):
        print(str(i*100/total) + "%", end="\r")
        second_flag = 0
        if (type(t_2[col_name_mz].iloc[i]) == np.float64) or ("_" not in t_2[col_name_mz].iloc[i]):
            mz = float(t_2[col_name_mz].iloc[i])
            rt = float(t_2[col_name_rt].iloc[i])
            if MS1_match(mz, rt, mz_list, rt_list, col_name_mz, col_name_rt, ppm_threshold, rt_tolerance):
                continue
            else:
                flag = 1
        else:
            for j in range(len(t_2[col_name_mz].iloc[i].split("_"))):
                mz = float(t_2[col_name_mz].iloc[i].split("_")[j])
                rt = float(t_2[col_name_rt].iloc[i].split("_")[j])
                nearest_id = mz_list.index(nearest_value(mz_list, mz))
                if MS1_match(mz, rt, mz_list, rt_list, col_name_mz, col_name_rt, ppm_threshold, rt_tolerance):
                    second_flag = 1
                    continue
            if second_flag == 0:
                flag = 1
        if flag == 1:
            flagged_rows.append(i)
    return t_2.drop(labels=flagged_rows, axis=0)

if __name__ == '__main__':

    logger.info('merging peak tables')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input_1', help="define the location of the first input peak table;", \
        default="table_1.csv", dest="input_1", required=False)
    parser.add_argument(
        '-j', '--input_2', help="define the location of the second input peak table;", \
        default="table_2.csv", dest="input_2", required=False)
    parser.add_argument(
        '-m', '--col_name_mz', help="define the column name of mz values;", \
        default="mz", dest="col_name_mz", required=False)
    parser.add_argument(
        '-r', '--col_name_rt', help="define the column name of retention time values;", \
        default="rt", dest="col_name_rt", required=False)
    parser.add_argument(
        '-p', '--ppm_threshold', help="define ppm threshold for merging;", \
        default="10", dest="ppm_threshold", required=False)
    parser.add_argument(
        '-t', '--rt_tolerance', help="define the retention time tolerance for merging;", \
        default="0.3", dest="rt_tolerance", required=False)
    parser.add_argument(
        '-o', '--output', help="define the location of output file;", \
        default="peak_merged.csv", required=False)

    args = parser.parse_args()

    table_1 = pd.read_csv(args.input_1, delimiter="\t")
    table_2 = pd.read_csv(args.input_2, delimiter="\t")
    peak_merged = overlapping_peaks(table_1, table_2, args.col_name_mz, args.col_name_rt, \
                                    args.ppm_threshold, args.rt_tolerance)
    peak_merged.to_csv(args.output, index=False)
