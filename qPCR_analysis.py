#!/usr/bin/env python

# sample command - >python Akutagawa_Jon_BME263_Graduate_Assignment_Week5.py Splice_Locations.bed Mus_musculus.GRCm38.dna.primary_assembly.fa


import xlrd
import sys
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches

plt.style.use('BME163.mplstyle')

qPCR_excel = pd.read_excel(sys.argv[1], sheetname=2, skiprows=35, header=0, parse_cols="B,D,E,I", na_values='Undetermined')

print (qPCR_excel.iloc[0,1])
print (qPCR_excel.columns.values)
#for index, row in qPCR_excel.iterrows():
    #print (row['Sample Name'],row['CT'])

def generateMeanList ():
    mean_dict = dict()
    old_sample_name = qPCR_excel.iloc[0,1]
    current_value_list = list()
    for index, row in qPCR_excel.iterrows():
        current_value = row['CT']
        if not pd.isnull(current_value):
            current_sample_name = row['Sample Name']
            if old_sample_name == current_sample_name:
                current_value_list.append(float(current_value))
                #print (current_value)
                old_sample_name = current_sample_name
            else:
                mean_dict[old_sample_name]=(sum(current_value_list)/len(current_value_list))
                current_value_list = [current_value]
                old_sample_name = current_sample_name
        else:
            print ('NAN HOMIE')

    mean_dict[old_sample_name]=(sum(current_value_list)/len(current_value_list))
    return (mean_dict)

def generateStandards (mean_dict):
    standards_dict = dict()
    mean_dict_keys = list(mean_dict.keys())
    gene_standards = list()
    housekeeping_standards = list()
    #print (mean_dict_keys)
    for i in range(0,len(mean_dict_keys)):
        #print (mean_dict_keys[i])
        if 'NCCIT' in mean_dict_keys[i]:
            concentration = np.log10(float(mean_dict_keys[i].split()[1]))
            primer = mean_dict_keys[i].split()[2]
            if primer == 'HLNC1':
                gene_standards.append((concentration, mean_dict[mean_dict_keys[i]]))
            else:
                housekeeping_standards.append((concentration, mean_dict[mean_dict_keys[i]]))
            #standards_dict[mean_dict_keys[i]] =
    #for i in range(0,len(standards_dict))
    #print (gene_standards, housekeeping_standards)

    gene_standards_y,gene_standards_x = zip(*gene_standards)
    housekeeping_standards_y,housekeeping_standards_x = zip(*housekeeping_standards)
    #housekeeping_standards_x = [x[0] for x in housekeeping_standards]

    gene_standards_coefficients = np.polyfit(gene_standards_x,gene_standards_y,1)
    housekeeping_standards_coefficients = np.polyfit(housekeeping_standards_x,housekeeping_standards_y,1)

    print ('yo')
    print (housekeeping_standards_x)
    print (housekeeping_standards_y)
    x_pos = np.arange(20,35,0.25)
    g_y_pos = gene_standards_coefficients[0] * x_pos + gene_standards_coefficients[1]
    hk_y_pos = housekeeping_standards_coefficients[0] * x_pos + housekeeping_standards_coefficients[1]
    plt.figure(1,figsize=(2,2))

    panel1=plt.axes([0.2,0.2,1.5/2,1.5/2])
    #panel2=plt.axes([2.1,0.15,1.5/4,1.5/2])

    panel1.scatter(gene_standards_x,gene_standards_y, \
                   #s=np.array(gene_length_list)/1000,\
                   s=8,\
                   marker='.',\
                   alpha=1,\
                   facecolor=(0,0,0),\
                   edgecolor='none',\
                   linewidth=0.1)


    panel1.scatter(housekeeping_standards_x,housekeeping_standards_y, \
                   #s=np.array(gene_length_list)/1000,\
                   s=8,\
                   marker='.',\
                   alpha=1,\
                   facecolor=(0,0,0),\
                   edgecolor='none',\
                   linewidth=0.1)

    panel1.plot(x_pos,g_y_pos,alpha = 0.5,label='HLNC1',color='blue')
    panel1.plot(x_pos,hk_y_pos,alpha = 0.5, label='ACTB',color='orange')
    plt.legend(fontsize = 4)
    panel1.set_xlabel('CT')
    panel1.set_ylabel('log10 (dilution factor)')
    gene = 'y = '+str(round(gene_standards_coefficients[0],4))+ ' * x + ' + str(round(gene_standards_coefficients[1],4))
    hk = 'y = '+str(round(housekeeping_standards_coefficients[0],4))+ ' * x + ' + str(round(housekeeping_standards_coefficients[1],4))
    panel1.text(22,-2,gene,color = 'blue',fontsize=4)
    panel1.text(26,0,hk,color = 'orange',fontsize=4)
    panel1.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on',\
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    """panel2.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on',\
                       right='off', labelright='off',\
                       top='off', labeltop='off')"""

    plt.savefig('standards_graphs.png')

    return (gene_standards_coefficients, housekeeping_standards_coefficients)

def calculateRVQ (qPCR_excel, gene_standards_coefficients,housekeeping_standards_coefficients):
    RVQ_list = list()
    for index, row in qPCR_excel.iterrows():
        if row['Target Name'] == 'HLNC1':
            score = gene_standards_coefficients[0] * row['CT'] + gene_standards_coefficients[1]
            RVQ_list.append(score)
        else:
            score = housekeeping_standards_coefficients[0] * row['CT'] + housekeeping_standards_coefficients[1]
            RVQ_list.append(score)
    #print (gene_standards_coefficients)
    #print (RVQ_list)
    qPCR_excel['RVQ'] = RVQ_list
    #qPCR_excel['10^RVQ'] = qPCR_excel['RVQ'].apply(10^)
    qPCR_excel['10^RVQ'] = np.power(10,qPCR_excel['RVQ'])
    #print (qPCR_excel.loc[:,['Sample Name','CT']])
    return (qPCR_excel)

#def f (x):
#    return

def normalize (qPCR_excel):
    print(qPCR_excel)
    #mean_qPCR_excel = qPCR_excel.groupby(['Sample Name','Target Name']).apply(np.mean)
    #mean_qPCR_excel = row['']
    #mean_qPCR_excel.columns.droplevel(0)
    sorted_excel = qPCR_excel.sort_values(['Target Name','Sample Name'])
    print (sorted_excel)
    #print (qPCR_excel.sort(['Target Name','Sample Name']))
    #print (len(qPCR_excel))


    old_name = sorted_excel.iloc[[0]]['Sample Name'].item()
    #print (old_name)
    #if old_name == 'Cy3 2.5 ACTB':
    #    print ('dope')

    count = 0
    replicate = 0
    sum_score = 0
    normalized_score_list = list()
    sample_dict = dict()
    norm_score_array = np.ones(96)

    for index, row in sorted_excel.iterrows():
        if count >= (len(qPCR_excel)/2):
            normalized_score_list.append((old_name,sample_dict))
            #print ((old_name,sample_dict))
            break
        #print (old_name)
        sample_name = row['Sample Name']
        #print (sample_name)
        housekeeping_score = sorted_excel.iloc[[count]]['10^RVQ'].values
        gene_score = sorted_excel.iloc[[count+48]]['10^RVQ'].values
        normalized_score = gene_score / housekeeping_score


        print ('sup')
        print (gene_score)
        print (housekeeping_score)

        print (normalized_score)
        norm_score_array[count+48]=(float(normalized_score))
        count += 1
        if old_name == sample_name:
            #print (row['Sample Name'])
            #print ('cool')
            #print (row['10^RVQ'])
            sample_dict[replicate] = normalized_score
            #score_array = normalized_score
            #np.r_(sample_np_array, normalized_score)
            #print (sample_np_array)
            sum_score += normalized_score
            replicate += 1
            old_name = sample_name
        else:
            #mean_score = np.mean(sample_dict.values())
            #sample_dict[replicate+1]=mean_score
            #print (sum_score)

            normalized_score_list.append((old_name,sample_dict))

            #normalized_score_list.append((old_name,sample_np_array))
            #print ('whatup')
            #print ((old_name,sample_dict))
            #print ((old_name,sample_np_array))
            #sample_np_array = np.array([])
            sample_dict=dict()
            replicate = 0
            sample_dict[replicate] = normalized_score
            replicate += 1
            old_name = sample_name

        #normalized_score_list.append((sample_name,sample_dict))
        #print (row)
    print (norm_score_array)
    print (len(norm_score_array))

    sorted_excel['Normalized 10^RVQ'] = norm_score_array
    print (sorted_excel)
    sorted_excel.to_csv('qPCR_analysis.csv', sep='\t')
        #print (count)
    #print (normalized_score_list)
    return(normalized_score_list)
    #mean_qPCR_excel.columns = ['Sample Name','Target Name','CT','RVQ','10^RVQ']
    #mean_qPCR_excel['Sample Name'] = mean_qPCR_excel['Sample Name'].replace(regex=True,inplace=True,to_replace=r'\D',value=r'')
    #print (mean_qPCR_excel.sort['Target','Sample Name'])

def printGraphs (normalized_scores):
    plt.figure(2,figsize=(5,5))

    panel1=plt.axes([0.15,0.2,3.75/5,3.75/5])
    #panel2=plt.axes([0.1,2,3.75/15,3.75/5])
    #panel3=plt.axes([4.15,0.2,3.75/15,3.75/5])
    #panel4=plt.axes([6.15,0.2,3.75/15,3.75/5])
    name_list = list()
    mean_list = list()
    for tup in normalized_scores:
        #print ('yo')
        name = tup[0].split()[0] + ' ' + tup[0].split()[1]
        name_list.append(name)
        sample_array = np.array(list(tup[1].values()))
        #print (sample_array)
        mean_score = np.nanmean(sample_array)
        if mean_score:
            mean_list.append(mean_score)
    print (mean_list)
    #print (name_list)

    locations = np.arange(len(mean_list))

    count = 0
    for i in range(0,len(mean_list)-4):
        left = locations[i]-.35
        width = 0.7
        height = mean_list[i]
        """
        if count < 2:
            rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                        linewidth=0.1,\
                                        facecolor=(0.5,0.5,0.5),\
                                        edgecolor=(0,0,0))
            count += 1
        else:
            rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                        linewidth=0.1,\
                                        facecolor=(0.5,0.5,0.5),\
                                        edgecolor=(0,0,0))
            count = 0
        """
        rectangle1=mplpatches.Rectangle((left,0),width,height,\
                                    linewidth=0.1,\
                                    facecolor=(0.5,0.5,0.5),\
                                    edgecolor=(0,0,0))
        panel1.add_patch(rectangle1)

    panel1.set_xlim([5,9])
    panel1.set_ylim([0,1])
    panel1.xaxis.set_ticks(np.arange(6,len(mean_list)-7,1))
    panel1.set_xticklabels(name_list[6:10],rotation=90)
    #panel1.set_xticks([name_list])
    #panel1.margins(0.25)

    panel1.set_xlabel('Sample')
    panel1.set_ylabel('Normalized 10^RVQ (Sample/Housekeeping)')
    panel1.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on',\
                       right='off', labelright='off',\
                       top='off', labeltop='off')

    plt.savefig('HLNC1-2_normalized_graphs.png')


def main():

    mean_dict = generateMeanList()
    gene_standards_coefficients, housekeeping_standards_coefficients = generateStandards(mean_dict)
    qPCR_excel_analysis = calculateRVQ(qPCR_excel,gene_standards_coefficients,housekeeping_standards_coefficients)
    normalized_scores = normalize (qPCR_excel_analysis)
    printGraphs(normalized_scores)

if __name__ == "__main__":
    main()

"""
# using xlrd
qPCR_excel = xlrd.open_workbook(sys.argv[1])
print (qPCR_excel.sheet_names())
results_sheet = qPCR_excel.sheet_by_index(2)

well_column = results_sheet.col_slice(colx=1,
                                start_rowx=36,
                                end_rowx=132)


sample_name_column = results_sheet.col_slice(colx=3,
                                start_rowx=36,
                                end_rowx=132)

CT_column = results_sheet.col_slice(colx=8,
                                start_rowx=36,
                                end_rowx=132)
"""

#for cell in well_column:
#    print (cell.value)

#for cell in CT_column:
    #print (cell.value)
"""
def generateMeanList ():
    mean_list = list()
    old_sample_name = sample_name_column[0].value
    current_value_list = list()
    for i in range(1,len(CT_column)):
        if row['CT'] == 'Undetermined':
            current_value = 'CT_column[i].value'
        else:
            current_value = CT_column[i].value
        current_sample_name = sample_name_column[i].value
        if old_sample_name == current_sample_name:
            #current_value_list.append(float(current_value))
            print (current_value)
            old_sample_name = current_sample_name
        else:
            #mean_list.append(sum(current_value_list)/len(current_value_list))
            current_value_list = [current_value]
            old_sample_name = current_sample_name

generateMeanList()
"""
