# Author: yoshi
# Date: 9/24/2020
# Updated: 10/08/2020
# Project: open TG gates
# Script: curate MAS5-processed csv file to generate matrix

import pandas as pd
import numpy as np
import os
import glob

def main():

    ## load attribute file
    file1 = '../data/Open-tggates_AllAttribute.tsv'
    print(f'[LOAD]: {file1}')
    att = pd.read_table(file1, sep='\t', header=0)
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

    att_ = att_[att_['SPECIES'] == 'Human'] # Select Human experiment data
    att_ = att_[att_['BARCODE'] != 'No ChipData'] # Drop unavailable barcode data

    ## generate a label column, such as 'CBP_2hr_Control'(Drug_TreatmentTime_Dose)
    labellist=[]
    for k, l, m, n in zip(att_['COMPOUND Abbr.'], att_['SACRI_PERIOD'], att_['DOSE_LEVEL'], att_['INDIVIDUAL_ID']):
        compound = k
        time = l.replace(' ', '')
        dose = m
        replicateID = str(n)
        label = compound + '_' + time + '_' + dose + '_' + replicateID
        labellist.append(label)
    #print(f'labellist: {len(labellist)}')
    #print(labellist)
    att_['LABEL'] = labellist # Add the label column 


    ### Select arbitrary dataset HERE!!! ###
    drugfile = '../data/Shared126DrugList_RatVivoVitro_HumanVitro.txt'
    print(f'[LOAD]: {drugfile}')
    druglist = []
    with open(drugfile, 'r') as f:
        for line in f:
            drug = line.strip()
            druglist.append(drug)
    print(f'Drug: {len(druglist)}')

    ## For example, use 1 code below to pick 8hr and 24hr ##
    #att_ = att_[(att_['SACRI_PERIOD'] == '8 hr') | (att_['SACRI_PERIOD'] == '24 hr')]
    att_ = att_[att_['COMPOUND_NAME'].isin(druglist)]
    att_ = att_[att_['SACRI_PERIOD'] == '24 hr']
    att_ = att_[(att_['DOSE_LEVEL'] == 'Control') | (att_['DOSE_LEVEL'] == 'High')]
    print(f'TimePoint: 24hr')
    print(f'Dose: Control, High')

    barcodelist = list(att_['BARCODE']) # Prep required barcodes list, such as barcodelist=['3016032003', '3016032004',..]

    ## prep barcode-label dict {'3016032003': 'GF_2hr_Low', '3016032018': 'GF_24hr_Low',...}
    barcode2label_mapping = {}
    for k, l in zip(att_['BARCODE'], att_['LABEL']):
        barcode2label_mapping[k] = l


    ## prepare for loading microarray data
    file2 = '../data/Human_vitro_liver_csv/*/celfiles/*csv'
    print(f'[LOAD] microarray data directory: {file2}')
    #dirlist = glob.glob('../data/Human_in_vitro/*/celfiles/*csv')
    dirlist = glob.glob(file2)
    datalist=[]
    for k in barcodelist:
        barcode_ = '00' + k + '.csv'
        for l in dirlist:
            filename = os.path.basename(l)
            if barcode_ in filename:
                print(f'file path: {l}')
                data = pd.read_table(l, sep=',', header=0)
                data_ = data.loc[:, ['probe_set_id', 'median_normalized_signal_intensity', 'mas5calls']]
                data_ = data_[data_['mas5calls'] == 'P'] # Select credible signal value
                data_ = data_.loc[:, ['probe_set_id', 'median_normalized_signal_intensity']]
                data_= data_.set_index('probe_set_id')
                data_s = data_['median_normalized_signal_intensity'] 
                data_s = data_s.rename(barcode2label_mapping[k]) # rename to identical label
                print(f'matrix: {data_s.shape}')
                datalist.append(data_s)
            #else:
            #    print(f'No matched file')

    # Generate a gene-by-sample matrix
    matrix = pd.concat(datalist, join='inner', axis=1) # remove probes with at least one 'NA'
    print(f'original matrix post concat: {matrix.shape}')

    ## load probe annotation file to convert probeID to GeneSymbol
    file3 = '../data/GPL570-55999.txt'
    print(f'[LOAD]: {file3}')
    anot = pd.read_table(file3, sep='\t', skiprows=16)
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

    outputfile = '../data/HumanVitroLiver_matrixTEST.txt'
    print(f'[SAVE]: {outputfile}')
    with open(outputfile, 'w') as f:
        matrix_log.to_csv(f, sep='\t', header=True, index=True)


if __name__ == '__main__':
    main()
