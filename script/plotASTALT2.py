# Author: yoshi
# Date: 10/21/2020
# Update:
# Project: oepnTGgates
# Script: to plot AST/ALT with dose depenency

import pandas as pd
import matplotlib.pyplot as plt

def main():
    filename1 = '../data/Open-tggates_AllAttribute.tsv'
    print(f'[LOAD]: {filename1}')
    data = pd.read_table(filename1, sep='\t', header=0)

    data_ = data[data['SPECIES'] == 'Rat']
    data_ = data_[data_['TEST_TYPE'] == 'in vivo']
    data_ = data_[data_['SIN_REP_TYPE'] == 'Single']
    data_ = data_[data_['ORGAN'] == 'Liver']

    data_ = data_.loc[:,['COMPOUND_NAME',
                         'SACRI_PERIOD',
                         'DOSE',
                         'DOSE_UNIT',
                         'DOSE_LEVEL',
                         'AST(IU/L)',
                         'ALT(IU/L)']]

    dosetable = data_.loc[:,['COMPOUND_NAME','DOSE_UNIT']]
    dosetable_ = dosetable[~dosetable.duplicated()] 

    data_3hr = data_[data_['SACRI_PERIOD'] == '3 hr']
    data_6hr = data_[data_['SACRI_PERIOD'] == '6 hr']
    data_9hr = data_[data_['SACRI_PERIOD'] == '9 hr']
    data_24hr = data_[data_['SACRI_PERIOD'] == '24 hr']

    #druglist=['colchicine'] # tentative code
    for k, l  in zip(dosetable_['COMPOUND_NAME'], dosetable_['DOSE_UNIT']):
    #for k in druglist: # tentative code
        drugname = k
        doseunit = l
        #doseunit = 'mg/kg' # tentative code

        data_3hr_ = data_3hr[data_3hr['COMPOUND_NAME'] == k]
        data_6hr_ = data_6hr[data_6hr['COMPOUND_NAME'] == k]
        data_9hr_ = data_9hr[data_9hr['COMPOUND_NAME'] == k]
        data_24hr_ = data_24hr[data_24hr['COMPOUND_NAME'] == k]

        ## plot AST/ALT
        fig = plt.figure(tight_layout=True)
        fig.suptitle(drugname, fontsize=10)

        ## 3hr (upper left) ##
        ax1_1 = fig.add_subplot(2,2,1)
        for i, (name, group) in enumerate(data_3hr_.groupby('DOSE_LEVEL')):
            #group.plot.scatter(x='DOSE', y='AST(IU/L)', alpha=1, s=7, c=[plt.get_cmap('Set1').colors[i]], marker='o', label=name, ax=ax1_1)
            group.plot.scatter(x='DOSE', y='AST(IU/L)', alpha=0.7, s=10, c='blue', marker='o', ax=ax1_1)
        xlabel1_1 = 'Dose (' + doseunit + ')'
        title1_1 = '3 hr'
        ylabel1_1 = 'AST(IU/L)'
        ax1_1.set_xlabel(xlabel1_1, fontsize=7)
        ax1_1.set_ylabel(ylabel1_1, fontsize=7)
        ax1_1.set_title(title1_1, fontsize=8)
        ax1_1.tick_params(axis='x', labelsize=5)
        ax1_1.tick_params(axis='y', labelsize=5)
        #ax1_1.legend(loc='upper left', fontsize=4, title='AST').get_title().set_fontsize(4)
        ax1_1.legend(loc='upper left', fontsize=5, labels=['AST'])

        ax1_2 = ax1_1.twinx()
        for i, (name, group) in enumerate(data_3hr_.groupby('DOSE_LEVEL')):
            #group.plot.scatter(x='DOSE', y='ALT(IU/L)', alpha=0.5, s=7, c=[plt.get_cmap('Set2').colors[i]], marker='D', label=name, ax=ax1_2)
            group.plot.scatter(x='DOSE', y='ALT(IU/L)', alpha=0.7, s=10, c='deeppink', marker='D', ax=ax1_2)
        ylabel1_2 = 'ALT(IU/L)'
        ax1_2.set_ylabel(ylabel1_2, fontsize=7)
        ax1_2.tick_params(axis='y', labelsize=5)
        #ax1_2.legend(loc='lower right', fontsize=4, title='ALT').get_title().set_fontsize(4)
        ax1_2.legend(loc='lower right', fontsize=5, labels=['ALT'])

        ## 6hr (upper right) ##
        ax2_1 = fig.add_subplot(2,2,2)
        for i, (name, group) in enumerate(data_6hr_.groupby('DOSE_LEVEL')):
            #group.plot.scatter(x='DOSE', y='AST(IU/L)', alpha=1, s=7, c=[plt.get_cmap('Set1').colors[i]], label=name, ax=ax2_1)
            group.plot.scatter(x='DOSE', y='AST(IU/L)', alpha=0.7, s=10, c='blue', marker='o', ax=ax2_1)
        xlabel2_1 = 'Dose (' + doseunit + ')'
        title2_1 = '6 hr'
        ylabel2_1 = 'AST(IU/L)'
        ax2_1.set_xlabel(xlabel2_1, fontsize=7)
        ax2_1.set_ylabel(ylabel2_1, fontsize=7)
        ax2_1.set_title(title2_1, fontsize=8)
        ax2_1.tick_params(axis='x', labelsize=5)
        ax2_1.tick_params(axis='y', labelsize=5)
        #ax2_1.legend(loc='upper left', fontsize=4, title='AST').get_title().set_fontsize(4)
        ax2_1.legend(loc='upper left', fontsize=5, labels=['AST'])

        ax2_2 = ax2_1.twinx()
        for i, (name, group) in enumerate(data_6hr_.groupby('DOSE_LEVEL')):
            #group.plot.scatter(x='DOSE', y='ALT(IU/L)', alpha=0.5, s=7, c=[plt.get_cmap('Set2').colors[i]], marker='D', label=name, ax=ax2_2)
            group.plot.scatter(x='DOSE', y='ALT(IU/L)', alpha=0.7, s=10, c='deeppink', marker='D', ax=ax2_2)
        ylabel2_2 = 'ALT(IU/L)'
        ax2_2.set_ylabel(ylabel2_2, fontsize=7)
        ax2_2.tick_params(axis='y', labelsize=5)
        #ax2_2.legend(loc='lower right', fontsize=4, title='ALT').get_title().set_fontsize(4)
        ax2_2.legend(loc='lower right', fontsize=5, labels=['ALT'])

        ## 9hr (lower left) ##
        ax3_1 = fig.add_subplot(2,2,3)
        for i, (name, group) in enumerate(data_9hr_.groupby('DOSE_LEVEL')):
            #group.plot.scatter(x='DOSE', y='AST(IU/L)', alpha=1, s=7, c=[plt.get_cmap('Set1').colors[i]], label=name, ax=ax3_1)
            group.plot.scatter(x='DOSE', y='AST(IU/L)', alpha=0.7, s=10, c='blue', marker='o', ax=ax3_1)
        xlabel3_1 = 'Dose (' + doseunit + ')'
        title3_1 = '9 hr'
        ylabel3_1 = 'AST(IU/L)'
        ax3_1.set_xlabel(xlabel3_1,fontsize=7)
        ax3_1.set_ylabel(ylabel3_1, fontsize=7)
        ax3_1.set_title(title3_1, fontsize=8)
        ax3_1.tick_params(axis='x', labelsize=5)
        ax3_1.tick_params(axis='y', labelsize=5)
        #ax3_1.legend(loc='upper left', fontsize=4, title='AST').get_title().set_fontsize(4)
        ax3_1.legend(loc='upper left', fontsize=5, labels=['AST'])

        ax3_2 = ax3_1.twinx()
        for i, (name, group) in enumerate(data_9hr_.groupby('DOSE_LEVEL')):
            #group.plot.scatter(x='DOSE', y='ALT(IU/L)', alpha=0.5, s=7, c=[plt.get_cmap('Set2').colors[i]], marker='D', label=name, ax=ax3_2)
            group.plot.scatter(x='DOSE', y='ALT(IU/L)', alpha=0.7, s=10, c='deeppink', marker='D', ax=ax3_2)
        ylabel3_2 = 'ALT(IU/L)'
        ax3_2.set_ylabel(ylabel3_2, fontsize=7)
        ax3_2.tick_params(axis='y', labelsize=5)
        #ax3_2.legend(loc='lower right', fontsize=4, title='ALT').get_title().set_fontsize(4)
        ax3_2.legend(loc='lower right', fontsize=4, labels=['ALT'])

        ## 24hr (lower right) ##
        ax4_1 = fig.add_subplot(2,2,4)
        for i, (name, group) in enumerate(data_24hr_.groupby('DOSE_LEVEL')):
            #group.plot.scatter(x='DOSE', y='AST(IU/L)', alpha=1, s=7, c=[plt.get_cmap('Set1').colors[i]], label=name, ax=ax4_1)
            group.plot.scatter(x='DOSE', y='AST(IU/L)', alpha=0.7, s=10, c='blue', marker='o', ax=ax4_1)
        xlabel4_1 = 'Dose (' + doseunit + ')'
        title4_1 = '24 hr'
        ylabel4_1 = 'AST(IU/L)'
        ax4_1.set_xlabel(xlabel4_1,fontsize=7)
        ax4_1.set_ylabel(ylabel4_1,fontsize=7)
        ax4_1.set_title(title4_1, fontsize=8)
        ax4_1.tick_params(axis='x', labelsize=5)
        ax4_1.tick_params(axis='y', labelsize=5)
        #ax4_1.legend(loc='upper left', fontsize=4, title='AST').get_title().set_fontsize(4)
        ax4_1.legend(loc='upper left', fontsize=5, labels=['AST'])

        ax4_2 = ax4_1.twinx()
        for i, (name, group) in enumerate(data_24hr_.groupby('DOSE_LEVEL')):
            #group.plot.scatter(x='DOSE', y='ALT(IU/L)', alpha=0.5, s=7, c=[plt.get_cmap('Set2').colors[i]], marker='D', label=name, ax=ax4_2)
            group.plot.scatter(x='DOSE', y='ALT(IU/L)', alpha=0.7, s=10, c='deeppink', marker='D', ax=ax4_2)
        ylabel4_2 = 'ALT(IU/L)'
        ax4_2.set_ylabel(ylabel4_2, fontsize=7)
        ax4_2.tick_params(axis='y', labelsize=5)
        #ax4_2.legend(loc='lower right', fontsize=4, title='ALT').get_title().set_fontsize(4)
        ax4_2.legend(loc='lower right', fontsize=5, labels=['ALT'])

        # Save
        filename2 = '../data/ASTALTplot/scatterplot_ASTALT_RatVivo_DoseDependent_'+ drugname + '.png'
        plt.savefig(filename2, dpi=300, format='png')
        print(f'[SAVE]: {filename2}')
        plt.clf()


if __name__ == '__main__':
    main()
