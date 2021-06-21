# THIS SOURCE CODE IS SUPPLIED “AS IS” WITHOUT WARRANTY OF ANY KIND, AND ITS
# AUTHOR AND THE JOURNAL OF MACHINE LEARNING RESEARCH (JMLR) AND JMLR’S
# PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL WARRANTIES, INCLUDING BUT
# NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE, AND ANY WARRANTIES OR NON INFRINGEMENT. THE USER ASSUMES
# ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND NEITHER THE
# AUTHOR NOR JMLR, NOR JMLR’S PUBLISHERS AND DISTRIBUTORS, WILL BE LIABLE FOR
# DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting the generality
# of the foregoing, neither the author, nor JMLR, nor JMLR’s publishers and
# distributors, warrant that the Source Code will be error-free, will operate
# without interruption, or will meet the needs of the user.

import numpy as np
import pandas as pd


COLUMNS = np.array(['TIMESTAMP_START', 'MONTH', 'SW_IN_F', 'TA_F', 'NEE_VUT_USTAR50', 'VPD_F', 'H_F_MDS', 'LE_F_MDS'])
VAR_NAMES = np.array(['DATETIME', 'MONTH', 'Rg', 'T', 'NEE', 'VPD', 'H', 'LE'])


def parse_date(x):
    date_str = x['TIMESTAMP_START'].astype(str)
    x['YEAR'] = int(date_str[:4])
    x['MONTH'] = int(date_str[4:6])
    x['DATE'] = int(date_str[:8])
    return x


def process_file(df, save_file, start=201404010000, end=201503312330):
    # Only keep observations from given time period
    df = df[(df.TIMESTAMP_START >= start) & (df.TIMESTAMP_START <= end)]
    df = df.apply(lambda x: parse_date(x), axis=1)

    columns = COLUMNS
    var_names = VAR_NAMES

    # Compute the daily maximum potential downward shortwave radiation
    df = df.merge(df[['DATE', 'SW_IN_POT']].groupby('DATE').max(), left_on='DATE', right_on='DATE', suffixes=('', '_MAX'))

    # Add noise to radiation when the potential is 0
    df.loc[df.SW_IN_POT == 0, 'SW_IN_F'] = np.random.normal(0, 0.000001, sum(df.SW_IN_POT == 0))

    # Filter out observations where potential radiation less than four fifths of the daily maximum
    df = df[df.SW_IN_POT >= (4 / 5) * df.SW_IN_POT_MAX]

    # Filter out gapfilled observations
    df = df[(df.NEE_VUT_USTAR50_QC == 0) & ((df.H_F_MDS_QC == 0) | (df.H_F_MDS_QC == -9999)) & ((df.LE_F_MDS_QC == 0) | (df.LE_F_MDS_QC == -9999))]

    # Filter columns
    df = df[columns]

    # Rename columns
    df = df.rename(columns=dict(zip(columns, var_names)))

    # Save to file without row names and without date
    df.drop(columns=['DATETIME']).to_csv(save_file, index=False)

    return df


df = pd.read_csv('../../FLX_FI-Hyy_FLUXNET2015_FULLSET_HH_1996-2018_beta-3.csv') # Path to file with full FLUXNET database.

process_file(df, './FI_HYY_Preprocessed_2014.csv')
process_file(df, './FI_HYY_Preprocessed_2013-2015.csv', start=201301010000, end=201512312330)
