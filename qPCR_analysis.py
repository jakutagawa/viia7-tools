#!/usr/bin/env python

# sample command - >python qPCR_analysis.py 2017-03-02_jon_edited.xls

#import xlrd
import sys
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.cm as cm

# stylesheet for graphs
plt.style.use('bw_graph.mplstyle')

# read in xls/xlsx file
def readExcelFile (filename):
    #print (filename)
    qPCR_excel = pd.read_excel(filename, sheetname=2, skiprows=35, header=0, parse_cols="B,D,E,I", na_values='Undetermined')
    #return 0
    return qPCR_excel

def generateMeanList (myExcelFile):
    mean_dict = dict()
    old_sample_name = myExcelFile.iloc[0,1]
    current_value_list = list()
    for index, row in myExcelFile.iterrows():
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

def generateStandards (mean_dict, standard_name, gene_name, housekeeping_name):
    standards_dict = dict()
    mean_dict_keys = list(mean_dict.keys())
    gene_standards = list()
    #gene2_standards = list()
    housekeeping_standards = list()
    #print (mean_dict)
    #print (gene_name)
    for i in range(0,len(mean_dict_keys)):
        #print (mean_dict_keys[i])
        if standard_name in mean_dict_keys[i]:
            concentration = np.log10(float(mean_dict_keys[i].split()[1]))
            primer = mean_dict_keys[i].split()[2]
            #print (primer)
            if primer == gene_name:
                gene_standards.append((concentration, mean_dict[mean_dict_keys[i]]))
            #elif primer == 'OCT4':
                #gene2_standards.append((concentration, mean_dict[mean_dict_keys[i]]))
            elif primer == housekeeping_name:
                housekeeping_standards.append((concentration, mean_dict[mean_dict_keys[i]]))
            #standards_dict[mean_dict_keys[i]] =
    #for i in range(0,len(standards_dict))
    #print (gene_standards, housekeeping_standards)


    HLNC1_standards = [(-0.69897000433601875, 25.29866663614909), (-2.6989700043360187, 33.68900108337402), (0.0, 22.18500010172526), (-1.6989700043360187, 30.017666498819988)]
    ACTB_standards = [(-0.69897000433601875, 26.45466677347819), (-1.6989700043360187, 30.544666926066082), (-2.6989700043360187, 34.641000747680664), (0.0, 23.806000391642254)]

    if not gene_standards:
        gene_standards = HLNC1_standards
    if not housekeeping_standards:
        housekeeping_standards = ACTB_standards

    #print ('here are the standards')
    #print (gene_standards)
    #print (housekeeping_standards)

    gene_standards_y,gene_standards_x = zip(*gene_standards)
    housekeeping_standards_y,housekeeping_standards_x = zip(*housekeeping_standards)

    #housekeeping_standards_x = [x[0] for x in housekeeping_standards]

    # use numpy.polyfit to fit a linear curve to standards
    gene_standards_coefficients = np.polyfit(gene_standards_x,gene_standards_y,1)
    housekeeping_standards_coefficients = np.polyfit(housekeeping_standards_x,housekeeping_standards_y,1)

    # set
    x_pos = np.arange(20,40,0.25)
    g_y_pos = gene_standards_coefficients[0] * x_pos + gene_standards_coefficients[1]
    hk_y_pos = housekeeping_standards_coefficients[0] * x_pos + housekeeping_standards_coefficients[1]

    fig_1 = plt.figure(figsize=(4,2))
    panel1=plt.axes([0.1,0.2,1.5/4,1.5/2])
    panel2=plt.axes([0.6,0.2,1.5/4,1.5/2])

    panel1.scatter(gene_standards_x,gene_standards_y, \
                   #s=np.array(gene_length_list)/1000,\
                   s=20,\
                   marker='.',\
                   alpha=1,\
                   facecolor=(0,0,0),\
                   edgecolor='none',\
                   linewidth=0.1)

    panel2.scatter(housekeeping_standards_x,housekeeping_standards_y, \
                   #s=np.array(gene_length_list)/1000,\
                   s=20,\
                   marker='.',\
                   alpha=1,\
                   facecolor=(0,0,0),\
                   edgecolor='none',\
                   linewidth=0.1)

    panel1.plot(x_pos,g_y_pos,alpha = 0.5,label=gene_name,color='blue')
    panel2.plot(x_pos,hk_y_pos,alpha = 0.5, label=housekeeping_name,color='orange')

    panel1.legend(fontsize = 4)
    panel2.legend(fontsize = 4)
    panel1.set_xlabel('CT')
    panel1.set_ylabel('log$\mathregular{_{10}}$ (dilution factor)')
    gene = 'y = '+str(round(gene_standards_coefficients[0],4))+ ' * x + ' + str(round(gene_standards_coefficients[1],4))
    hk = 'y = '+str(round(housekeeping_standards_coefficients[0],4))+ ' * x + ' + str(round(housekeeping_standards_coefficients[1],4))
    panel1.text(22,-2,gene,color = 'blue',fontsize=4)
    panel2.text(26,0,hk,color = 'orange',fontsize=4)
    panel1.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on',\
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panel2.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on',\
                       right='off', labelright='off',\
                       top='off', labeltop='off')

    plt.savefig('standards_graphs.png')

    return (gene_standards_coefficients, housekeeping_standards_coefficients)
    #return 0,0

def calculateRVQ (qPCR_excel, gene_standards_coefficients,housekeeping_standards_coefficients, gene_name, housekeeping_name):
    RVQ_list = list()
    #sorted_excel = qPCR_excel[(qPCR_excel['Target Name']==gene_name) or (pdata['Target Name']==housekeeping_name)]
    gene_name_rows = qPCR_excel[qPCR_excel['Target Name']==gene_name]
    housekeeping_name_rows = qPCR_excel[qPCR_excel['Target Name']==housekeeping_name]
    filtered_excel = gene_name_rows.append(housekeeping_name_rows)
    #print (gene_name)
    print (filtered_excel)
    for index, row in filtered_excel.iterrows():
        if row['Target Name'] == gene_name:
            score = gene_standards_coefficients[0] * row['CT'] + gene_standards_coefficients[1]
            RVQ_list.append(score)
        elif row['Target Name'] == housekeeping_name:
            score = housekeeping_standards_coefficients[0] * row['CT'] + housekeeping_standards_coefficients[1]
            RVQ_list.append(score)
    #print (gene_standards_coefficients)
    print (RVQ_list)
    print (len(RVQ_list))
    print (len(filtered_excel))
    #print (len(RVQ_list))
    filtered_excel['RVQ'] = RVQ_list

    filtered_excel['10^RVQ'] = np.power(10,filtered_excel['RVQ'])
    #print (filtered_excel)
    #print (qPCR_excel.loc[:,['Sample Name','CT']])

    return (filtered_excel)


#def f (x):
#    return

def normalize (qPCR_excel, gene_name, housekeeping_name, standard_name):
    #print ('were here')
    #print (qPCR_excel)
    #mean_qPCR_excel = qPCR_excel.groupby(['Sample Name','Target Name']).apply(np.mean)
    #mean_qPCR_excel = row['']
    #mean_qPCR_excel.columns.droplevel(0)
    sorted_excel = qPCR_excel.sort_values(['Target Name','Sample Name'])
    #print (sorted_excel)
    #sorted_excel = sorted_excel['']

    category_list = list()
    valid_entries = 0
    for index, row in sorted_excel.iterrows():
        #if row['Target Name'] == gene_name:
        if standard_name in row['Sample Name']:
            category_list.append('other')
        elif row['Target Name'] == housekeeping_name:
            category_list.append('housekeeping')
            valid_entries += 1
        elif row['Target Name'] == gene_name:
            category_list.append('gene')
            valid_entries += 1

    #print (qPCR_excel.sort(['Target Name','Sample Name']))
    #print (len(qPCR_excel))
    sorted_excel['category'] = category_list
    #print (len(sorted_excel))
    print (sorted_excel)
    filtered_excel = sorted_excel[sorted_excel['category'] != 'other']
    print (filtered_excel)


    old_name = filtered_excel.iloc[[0]]['Sample Name'].item()
    #print (old_name)
    #if old_name == 'Cy3 2.5 ACTB':
    #    print ('dope')

    count = 0
    replicate = 0
    sum_score = 0
    normalized_score_list = list()
    sample_dict = dict()
    norm_score_array = np.ones(valid_entries)

    for index, row in filtered_excel.iterrows():
        if count >= (valid_entries/2):
            normalized_score_list.append((old_name,sample_dict))
            #print ((old_name,sample_dict))
            break
        #print (old_name)
        sample_name = row['Sample Name']
        #print (sample_name)
        housekeeping_score = filtered_excel.iloc[[count]]['10^RVQ'].values
        gene_score = filtered_excel.iloc[[count+valid_entries/2]]['10^RVQ'].values
        normalized_score = gene_score / housekeeping_score

        norm_score_array[count+valid_entries/2]=(float(normalized_score))
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

    filtered_excel['Normalized 10^RVQ'] = norm_score_array
    print (sorted_excel)
    sorted_excel.to_csv('qPCR_analysis.csv', sep='\t')
        #print (count)
    #print (normalized_score_list)
    return(normalized_score_list)
    #mean_qPCR_excel.columns = ['Sample Name','Target Name','CT','RVQ','10^RVQ']
    #mean_qPCR_excel['Sample Name'] = mean_qPCR_excel['Sample Name'].replace(regex=True,inplace=True,to_replace=r'\D',value=r'')
    #print (mean_qPCR_excel.sort['Target','Sample Name'])

def printGraphs (normalized_scores, gene_name):
    plt.figure(2,figsize=(5,5))

    panel1=plt.axes([0.15,0.2,3.75/5,3.75/5])
    #panel2=plt.axes([0.1,2,3.75/15,3.75/5])
    #panel3=plt.axes([4.15,0.2,3.75/15,3.75/5])
    #panel4=plt.axes([6.15,0.2,3.75/15,3.75/5])
    name_list = list()
    group_list = list()
    mean_list = list()

    for tup in normalized_scores:
        #print ('yo')
        name = tup[0].split()[0] + ' ' + tup[0].split()[1]
        print (name)

        sample_array = np.array(list(tup[1].values()))
        # find the mean and disregard missing points
        mean_score = np.nanmean(sample_array)
        if not ('NCCIT' in name):
            if mean_score:
                name_list.append(name)
                group_list.append(name.split()[0])
                mean_list.append(mean_score)
    print ('check this out')
    print (mean_list)
    print (name_list)

    locations = np.arange(len(mean_list))

    unique_groups = set(group_list)
    unique_groups_list = sorted(list(unique_groups) )            # == set(['a', 'b', 'c'])
    unique_groups_count = len(unique_groups)
    print (unique_groups_list)
    print (name_list)
    print (unique_groups_count)
    #plt.jet()
    #plt.set_cmap('jet')
    colors_subsection = np.linspace(0, 1, len(unique_groups))
    colors = [cm.jet(x) for x in colors_subsection]
    #print (colors)

    count = 0
    group_count = 0
    temp_group = unique_groups_list[0]

    for i in range(0,len(mean_list)):
        left = locations[i]-.35
        width = 0.7
        height = mean_list[i]

        # iterate through group colors
        if not (unique_groups_list[group_count] in name_list[i]):
            group_count += 1

        rectangle1=mplpatches.Rectangle((left,0),width,height,\
                                    linewidth=0.1,\
                                    facecolor=colors[group_count],\
                                    edgecolor=(0,0,0))
        panel1.add_patch(rectangle1)

    #panel1.set_xlim([-1,9])
    #panel1.set_ylim([0,1])
    #panel1.xaxis.set_ticks(np.arange(6,len(mean_list)-7,1))
    #panel1.set_xticklabels(name_list[6:10],rotation=90)

    y_max = max(mean_list) * 1.25
    panel1.set_xlim([-1,len(mean_list)])
    panel1.set_ylim([0,y_max])
    panel1.xaxis.set_ticks(np.arange(0,len(mean_list),1))
    panel1.set_xticklabels(name_list,rotation=90)
    #panel1.set_xticks([name_list])
    #panel1.margins(0.25)
    #print (gene_name)
    #raph_title = str(gene_name) + ' Expression'
    panel1.set_title(str(gene_name) + ' Expression')
    panel1.set_xlabel('Sample')
    panel1.set_ylabel('Normalized 10^RVQ (Sample/Housekeeping)')
    panel1.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on',\
                       right='off', labelright='off',\
                       top='off', labeltop='off')

    plt.savefig('normalized_graphs.png')


class CommandLine() :
    '''
    modified from David Bernick

    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']
    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import argparse

        self.parser = argparse.ArgumentParser(description = 'Program prolog - Specify the input excel file, gene names, and standards.',
                                              epilog = 'Program epilog - parameters of infile and outfile must be given.',
                                              add_help = True, #default is True
                                              prefix_chars = '-',
                                              usage = '%(prog)s -f | -g | -hk | -s'
                                              )
        self.parser.add_argument('-f', '--excelfile', dest='location',
                                 default='2017-03-12_jon.xls', action='store', type=str,
                                 help='specify the filename')
        self.parser.add_argument('-g', '--geneOfInterest', dest='gene_name',
                                 default='HLNC1', action='store', type=str,
                                 help='specify the gene of interest')
        self.parser.add_argument('-hk', '--housekeepingGene', default='ACTB',
                                 dest='housekeeping_name', action='store', type=str,
                                 help='specify the housekeeping gene')
        self.parser.add_argument('-s', '--standard', default='standard',
                                 dest='standard_name', action='store', type=str,
                                 help='specify tag to identify standards')


        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

        #print (self.args)

def main(myCommandLine=None):

    myCommandLine = CommandLine(myCommandLine)
    #qPCR_excel = pd.read_excel(sys.argv[1], sheetname=2, skiprows=35, header=0, parse_cols="B,D,E,I", na_values='Undetermined')
    #print (myCommandline.args)
    #print (myCommandLine.args.location)
    myExcelFile = readExcelFile(myCommandLine.args.location)
    mean_dict = generateMeanList(myExcelFile)
    #print (mean_dict)

    gene_standards_coefficients, housekeeping_standards_coefficients = generateStandards(mean_dict, myCommandLine.args.standard_name, myCommandLine.args.gene_name, myCommandLine.args.housekeeping_name)
    #print (gene_standards_coefficients)
    #print (housekeeping_standards_coefficients)
    qPCR_excel_analysis = calculateRVQ(myExcelFile, gene_standards_coefficients,housekeeping_standards_coefficients, myCommandLine.args.gene_name, myCommandLine.args.housekeeping_name)
    normalized_scores = normalize (qPCR_excel_analysis, myCommandLine.args.gene_name, myCommandLine.args.housekeeping_name, myCommandLine.args.standard_name)
    #print (len(normalized_scores))
    printGraphs(normalized_scores, myCommandLine.args.gene_name)

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
