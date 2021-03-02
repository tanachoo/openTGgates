# Author: yoshi
# Date: 10/08/2020
# Updated: 03/02/2021
# Project: openTGgates
# Script: curate MAS5-processed csv file to generate dataset for BN input
# Run: python make_dataset.py --species <rat/human> --exptype <vivo/vitro> --reptype <single/repeat> --output <../data/TESTdataset_acetaminophen.txt>

import pandas as pd
import numpy as np
import os
import glob
import argparse

def main():

    ## load attribute file
    master_file = '../data/Open-tggates_AllAttribute.tsv'
    print(f'[LOAD] attritute file: {master_file}')
    att = pd.read_table(master_file, sep='\t', header=0)
    att_ = att.loc[:, ['BARCODE',
                       'INDIVIDUAL_ID',
                       'ORGAN',
                       'COMPOUND_NAME',
                       'COMPOUND Abbr.',
                       'SPECIES',
                       'TEST_TYPE',
                       'SIN_REP_TYPE',
                       'SACRI_PERIOD',
                       'DOSE_LEVEL']]


    ## Select Rat/Human and vivo/vitro
    if args.species == 'rat':
        att_ = att_[att_['SPECIES'] == 'Rat'] # Select 'Rat' experiment data

        if args.exptype == 'vivo':
            att_ = att_[att_['TEST_TYPE'] == 'in vivo'] # Select 'in vivo' design

            if args.reptype == 'single':
                att_ = att_[att_['SIN_REP_TYPE'] == 'Single'] # Select 'Single' administration design
                att_ = att_[att_['ORGAN'] == 'Liver'] # Select 'Liver' design
                att_ = att_[att_['BARCODE'] != 'No ChipData'] # Drop unavailable barcode data
                print(f'[EXP DESIGN]: Rat / in vivo / liver / single')

            elif args.reptype == 'repeat':
                att_ = att_[att_['SIN_REP_TYPE'] == 'Repeat'] # Select 'Repeat' administration design
                att_ = att_[att_['ORGAN'] == 'Liver'] # Select 'Liver' design
                att_ = att_[att_['BARCODE'] != 'No ChipData'] # Drop unavailable barcode data
                print(f'[EXP DESIGN]: Rat / in vivo / liver / repeat')

        elif args.exptype == 'vitro':
            att_ = att_[att_['TEST_TYPE'] == 'in vitro'] # Select 'in vitro' design
            att_ = att_[att_['BARCODE'] != 'No ChipData'] # Drop unavailable barcode data
            print(f'[EXP DESIGN]: Rat / vitro')

    elif args.species == 'human':
        att_ = att_[att_['SPECIES'] == 'Human'] # Select Human experiment data
        att_ = att_[att_['BARCODE'] != 'No ChipData'] # Drop unavailable barcode data
        print(f'[EXP DESIGN]: Human / vitro')


    ## generate a label column, such as 'CBP_2hr_Control_1'('Drug_TreatmentTime_Dose_replicates')
    labellist=[]
    for k, l, m, n in zip(att_['COMPOUND Abbr.'], att_['SACRI_PERIOD'], att_['DOSE_LEVEL'], att_['INDIVIDUAL_ID']):
        compound = k
        time = l.replace(' ', '')
        dose = m
        replicateID = str(n)
        label = compound + '_' + time + '_' + dose + '_' + replicateID
        labellist.append(label)
    att_['LABEL'] = labellist # Add the label column 


    ### Here, select arbitrary conditions (Dose/Time/Drug) ###
    ### If you do NOT select any criteria, you'll get all the number of samples ###
    '''
    ## If you use drug txt dataset...Run code below
    drugfile = '../data/Shared126DrugList_RatVivoVitro_HumanVitro.txt'
    print(f'[LOAD]: {drugfile}')
    druglist = []
    with open(drugfile, 'r') as f:
        for line in f:
            drug = line.strip()
            druglist.append(drug)
    print(f'Drug: {len(druglist)}')

    ## ex.
    #att_ = att_[att_['COMPOUND_NAME'].isin(druglist)]
    att_ = att_[att_['SACRI_PERIOD'] == '24 hr']
    att_ = att_[(att_['DOSE_LEVEL'] == 'Control') | (att_['DOSE_LEVEL'] == 'High')]
    print(f'TimePoint: 24hr')
    print(f'Dose: Control, High')
    '''
    # att_ = att_[att_['COMPOUND_NAME'] == 'acetaminophen']


    # Prep required barcodes list, such as barcodelist=['3016032003', '3016032004',..]
    barcodelist = list(att_['BARCODE'])


    ## Prep barcode-label dict {'3016032003': 'GF_2hr_Low_1', '3016032018': 'GF_24hr_Low_1',...}
    barcode2label_mapping = {}
    for k, l in zip(att_['BARCODE'], att_['LABEL']):
        barcode2label_mapping[k] = l


    ## Prepare for loading microarray data
    if args.species == 'rat':
        if args.exptype == 'vivo':
            if args.reptype == 'single':
                dirpath = '../data/Rat_vivo_liver_single_csv/*/celfiles/*csv'
            elif args.reptype == 'repeat':
                dirpath = '../data/Rat_vivo_liver_repeat_csv/*/celfiles/*csv'
        elif args.exptype == 'vitro':
            dirpath = '../data/Rat_vitro_liver_csv/*/celfiles/*csv'
    elif args.species == 'human':
        dirpath = '../data/Human_vitro_liver_csv/*/celfiles/*csv'

    print(f'[LOAD] microarray data directory: {dirpath}')
    dirlist = glob.glob(dirpath)
    datalist=[]
    for k in barcodelist:
        barcode_ = '00' + k + '.csv'
        for l in dirlist:
            filename = os.path.basename(l)
            if barcode_ in filename:
                # print(f'file path: {l}')
                data = pd.read_table(l, sep=',', header=0)
                data_ = data.loc[:, ['probe_set_id', 'median_normalized_signal_intensity', 'mas5calls']]
                data_ = data_[data_['mas5calls'] == 'P'] # Select credible signal value
                data_ = data_.loc[:, ['probe_set_id', 'median_normalized_signal_intensity']]
                data_= data_.set_index('probe_set_id')
                data_s = data_['median_normalized_signal_intensity'] 
                data_s = data_s.rename(barcode2label_mapping[k]) # rename to identical label
                # print(f'matrix: {data_s.shape}')
                datalist.append(data_s)

    # Generate a gene-by-sample matrix
    #matrix = pd.concat(datalist, join='inner', axis=1) # remove probes with at least one 'NA'
    matrix = pd.concat(datalist, join='outer', axis=1) # keep all probes with at least one 'NA'
    print(f'original matrix post concat: {matrix.shape}')

    ## load probe annotation file to convert probeID to GeneSymbol
    if args.species == 'rat':
        anot_file = '../data/GPL1355-10794.txt' # Rat annotation file
        print(f'[LOAD] Rat annotation file: {anot_file}')
    elif args.species == 'human':
        anot_file = '../data/GPL570-55999.txt' # Human annotation file
        print(f'[LOAD] Human annotation file: {anot_file}')

    anot = pd.read_table(anot_file, sep='\t', skiprows=16)
    anot_ = anot.loc[:, ['ID', 'SPOT_ID', 'Gene Symbol']]
    #print(f'{anot_.isnull().sum()}')
    anot_ = anot_.fillna({'Gene Symbol': 'NA'}) # Convert blank 'Gene Symbol' to NA
    anot_naSymbol = anot_[anot_['Gene Symbol'] == 'NA'] # extract GeneSymbol with NA probe
    noGeneSymbol_list = list(anot_naSymbol['ID'])
    anot_control = anot_[anot_['SPOT_ID'] == '--Control'] # extract control probe
    controlProbe_list = list(anot_control['ID'])

    ## prep probe-gene dict
    probe2gene_mapping = {} #dict = {probe: Gene symbol, ...}
    for k, l in zip(anot_['ID'], anot_['Gene Symbol']):
        probe2gene_mapping[k] = l

    ## Remove unrequired probe from matrix
    for i in noGeneSymbol_list:
        if i in list(matrix.index):
            matrix = matrix.drop(i, axis=0)
    print(f'matrix post removal of noGeneSymbol probe: {matrix.shape}')

    for i in controlProbe_list:
        if i in list(matrix.index):
            matrix = matrix.drop(i, axis=0)
    print(f'matrix post removal of controlProbe probe: {matrix.shape}')

    ## obtain final probe list
    probelist = list(matrix.index)

    ## Convert probeID --> GeneSymbol
    probe2gene=[]
    for i in probelist:
        gene = probe2gene_mapping[i]
        probe2gene.append(gene)

    matrix.index = probe2gene # Relabel index name to Gene
    matrix_ = matrix.reset_index()
    matrix_ = matrix_.rename(columns={'index': 'Gene'})
    matrix_ = matrix_.groupby('Gene').mean() # take average for duplicated gene signal
    print(f'matrix shape post conversion of duplicated same gene: {matrix_.shape}')
    matrix_log = np.log2(matrix_ + 1) # take log2
    print(f'final matrix: {matrix_log.shape}')

    ## Export dataset as txt file
    print(f'[SAVE]: {args.output}')
    with open(args.output, 'w') as f:
        matrix_log.to_csv(f, sep='\t', header=True, index=True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--species', type=str, help="rat / human")
    parser.add_argument('--exptype', type=str, help="vivo / vitro")
    parser.add_argument('--reptype', type=str, help="single / repeat")
    parser.add_argument('--output', type=str, help="Set output filename")
    args = parser.parse_args()

    main()

