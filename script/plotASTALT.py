# Author: yoshi
# Date: 10/5/2020
# Updated: 10/6/2020
# Project: openTGgates
# Script: plot AST/ALT relation

import pandas as pd
import matplotlib.pyplot as plt

def main():

    filename1 = '../data/Open-tggates_AllAttribute.tsv'
    print(f'[LOAD]: {filename1}')
    data = pd.read_table(filename1, sep='\t', header=0)
    data_ = data[data['SPECIES'] == 'Rat'] # pick 'Rat'
    data_ = data_[data_['TEST_TYPE'] == 'in vivo'] # pick 'in vitro'
    data_ = data_[data_['SIN_REP_TYPE'] == 'Single'] # pick 'single' administration
    data_ = data_[data_['ORGAN'] == 'Liver'] # pick 'liver'
    #data_ = data_[data_['COMPOUND_NAME'] == 'colchicine']

    ## dose course
    data_ctl = data_[data_['DOSE_LEVEL'] == 'Control']
    data_low = data_[data_['DOSE_LEVEL'] == 'Low']
    data_mid = data_[data_['DOSE_LEVEL'] == 'Middle']
    data_high = data_[data_['DOSE_LEVEL'] == 'High']

    
    ## Fix timepoint (inc/ all timepoints) ##
    ## AST-ALT plot with dose dependent

    # high
    for i, (name, group) in enumerate(data_high.groupby('SACRI_PERIOD')):
        plt.scatter(group['AST(IU/L)'], group['ALT(IU/L)'], alpha=0.4, s=10, c=[plt.get_cmap('tab10').colors[i]], marker='h', label=name)
    # middle
    for i, (name, group) in enumerate(data_mid.groupby('SACRI_PERIOD')):
        plt.scatter(group['AST(IU/L)'], group['ALT(IU/L)'], alpha=0.4, s=10, c=[plt.get_cmap('tab10').colors[i]], marker='^', label=name)
    # low
    for i, (name, group) in enumerate(data_low.groupby('SACRI_PERIOD')):
        plt.scatter(group['AST(IU/L)'], group['ALT(IU/L)'], alpha=0.4, s=10, c=[plt.get_cmap('tab10').colors[i]], marker='D', label=name)
    # control
    for i, (name, group) in enumerate(data_ctl.groupby('SACRI_PERIOD')):
        plt.scatter(group['AST(IU/L)'], group['ALT(IU/L)'], alpha=0.4, s=10, c=[plt.get_cmap('tab10').colors[i]], marker='o', label=name)


    # high
    plt.scatter(data_high['AST(IU/L)'], data_high['ALT(IU/L)'], alpha=0.4, s=10, marker='h', c='magenta', label='High')
    for x, y, name in zip(data_high['AST(IU/L)'], data_high['ALT(IU/L)'], data_high['COMPOUND_NAME']):
        plt.text(x, y, name, alpha=0.8, color='magenta', fontsize=3)
    # mid
    plt.scatter(data_mid['AST(IU/L)'], data_mid['ALT(IU/L)'], alpha=0.4, s=10, marker='^', c='orange', label='Middle')
    for x, y, name in zip(data_mid['AST(IU/L)'], data_mid['ALT(IU/L)'], data_mid['COMPOUND_NAME']):
        plt.text(x, y, name, alpha=0.8, color='orange', fontsize=3)
    # low
    plt.scatter(data_low['AST(IU/L)'], data_low['ALT(IU/L)'], alpha=0.4, s=10, marker='D', c='blue', label='Low')
    for x, y, name in zip(data_low['AST(IU/L)'], data_low['ALT(IU/L)'], data_low['COMPOUND_NAME']):
        plt.text(x, y, name, alpha=0.8, color='blue', fontsize=3)
    # control
    plt.scatter(data_ctl['AST(IU/L)'], data_ctl['ALT(IU/L)'], alpha=0.4, s=10, marker='o', c='green', label='Control')
    for x, y, name in zip(data_ctl['AST(IU/L)'], data_ctl['ALT(IU/L)'], data_ctl['COMPOUND_NAME']):
        plt.text(x, y, name, alpha=0.8, color='green', fontsize=3)


    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('AST(IU/L)')
    plt.ylabel('ALT(IU/L)')
    #plt.title('Rat Vivo, Dose Dependend Test')
    plt.legend(loc='lower right')
    filename2 = '../data/scatterplot_ASTALT_RatVivo_DoseDependent_TEST.png'
    plt.savefig(filename2, dpi=600, format='png')
    print(f'[SAVE]: {filename2}')
    plt.clf()



    ## timepoint course
    data_3hr = data_[data_['SACRI_PERIOD'] == '3 hr']  
    data_6hr = data_[data_['SACRI_PERIOD'] == '6 hr']  
    data_9hr = data_[data_['SACRI_PERIOD'] == '9 hr']  
    data_24hr = data_[data_['SACRI_PERIOD'] == '24 hr']  

    ## Fix dose level ##
    ## AST-ALT plot with time dependent

    # 24hr
    plt.scatter(data_24hr['AST(IU/L)'], data_24hr['ALT(IU/L)'], alpha=0.4, s=10, marker='h', c='magenta', label='24hr')
    for x, y, name in zip(data_24hr['AST(IU/L)'], data_24hr['ALT(IU/L)'], data_24hr['COMPOUND_NAME']):
        plt.text(x, y, name, alpha=0.8, color='magenta', fontsize=3)
    # 9hr
    plt.scatter(data_9hr['AST(IU/L)'], data_9hr['ALT(IU/L)'], alpha=0.4, s=10, marker='^', c='orange', label='9hr')
    for x, y, name in zip(data_9hr['AST(IU/L)'], data_9hr['ALT(IU/L)'], data_9hr['COMPOUND_NAME']):
        plt.text(x, y, name, alpha=0.8, color='orange', fontsize=3)
    # 6hr
    plt.scatter(data_6hr['AST(IU/L)'], data_6hr['ALT(IU/L)'], alpha=0.4, s=10, marker='D', c='blue', label='6hr')
    for x, y, name in zip(data_6hr['AST(IU/L)'], data_6hr['ALT(IU/L)'], data_6hr['COMPOUND_NAME']):
        plt.text(x, y, name, alpha=0.8, color='blue', fontsize=3)
    # 3hr
    plt.scatter(data_3hr['AST(IU/L)'], data_3hr['ALT(IU/L)'], alpha=0.4, s=10, marker='o', c='green', label='3hr')
    for x, y, name in zip(data_3hr['AST(IU/L)'], data_3hr['ALT(IU/L)'], data_3hr['COMPOUND_NAME']):
        plt.text(x, y, name, alpha=0.8, color='green', fontsize=3)

    plt.xscale('log')
    plt.yscale('log')
    plt.title('Rat Vivo, Time Dependent Test')
    plt.xlabel('AST(IU/L)')
    plt.ylabel('ALT(IU/L)')
    plt.legend(loc='lower right')
    filename3 = '../data/scatterplot_ASTALT_RatVivo_TimeDependent_TEST1.png'
    plt.savefig(filename3, dpi=600, format='png')
    print(f'[SAVE]: {filename3}')
    plt.clf()




    data_24hr_ = data_[data_['SACRI_PERIOD'] == '3 hr']
    data_24hr_ctl = data_24hr_[data_24hr_['DOSE_LEVEL'] == 'Control']
    data_24hr_low = data_24hr_[data_24hr_['DOSE_LEVEL'] == 'Low']
    data_24hr_middle = data_24hr_[data_24hr_['DOSE_LEVEL'] == 'Middle']
    data_24hr_high = data_24hr_[data_24hr_['DOSE_LEVEL'] == 'High']

    plt.scatter(data_24hr_high['AST(IU/L)'], data_24hr_high['ALT(IU/L)'], alpha=0.4, s=20, marker='o', c='magenta', label='High')
    plt.scatter(data_24hr_middle['AST(IU/L)'], data_24hr_middle['ALT(IU/L)'], alpha=0.4, s=20, marker='o', c='orange', label='Middle')
    plt.scatter(data_24hr_low['AST(IU/L)'], data_24hr_low['ALT(IU/L)'], alpha=0.4, s=20, marker='o', c='blue', label='Low')
    plt.scatter(data_24hr_ctl['AST(IU/L)'], data_24hr_ctl['ALT(IU/L)'], alpha=0.4, s=20, marker='o', c='green', label='Control')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Rat in vivo, TimePoint: 3hr')
    plt.xlabel('AST(IU/L)')
    plt.ylabel('ALT(IU/L)')
    plt.legend(loc='lower right')
    filename4 = '../data/scatterplot_AST_ALT_ratvivo_3hr_dose.png'
    plt.savefig(filename4, dpi=300, format='png')
    print(f'[SAVE]: {filename4}')
    plt.clf()


if __name__ == '__main__':
    main()
