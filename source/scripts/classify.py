import cPickle
import argparse
import numpy as np
import pandas as pd
from sklearn import ensemble

parser = argparse.ArgumentParser(description='ML module for Pegasus')
parser.add_argument('-i','--infile', help='path to the transcript annotation report', required=True)
parser.add_argument('-m','--modelfile', help='path to theserialized gradient boosting classifier object', required=True)
parser.add_argument('-o','--outfile', help='path to the desired final Pegasus report', required=True)
parser.add_argument('-l','--logdir', help='directory for writing log file', required=True)
args = vars(parser.parse_args())

report_in = args['infile']
classifier = args['modelfile']
report_out = args['outfile']
log_dir = args['logdir']

f_log = open(log_dir + '/log_classify', 'w')
f_log.write('START\n')
f_log.write('classifier located at {0}\n'.format(classifier))
f_log.write('reading input file from {0}\n'.format(report_in))
x_input = pd.read_csv(report_in, header=0, sep='\t')
f_log.write('successfully loaded input file with shape {0}\n'.format(x_input.shape))
clf = cPickle.load(open(classifier))
f_log.write('successfully loaded classifier')

features = ['5p-kinase', '3p-kinase', 'in-frame', 'premature stop codon', 
            '5p-intron', '3p-intron', '5p-BPR intergenic', '5p-BPR CDS', '3p-BPR intergenic', '3p-BPR CDS']

mini = pd.DataFrame(np.zeros((len(x_input),len(features))), columns=features)

for row in x_input.index:
    cols = list(x_input.ix[row,:])
    chr1 = str(cols[5]).strip()
    chr2 = str(cols[6]).strip()
    kinase = str(cols[22]).strip()
    frame = str(cols[25]).strip()
    aa_seq = str(cols[30]).strip()
    transcript_5p = str(cols[31]).strip()
    transcript_3p = str(cols[32]).strip()
    BPR_5p = str(cols[33]).strip()
    BPR_3p = str(cols[34]).strip()
    
    if kinase in ['5p_KINASE', 'BOTH_KINASE']:
        mini.ix[row,'5p-kinase'] = 1
    if kinase in ['3p_KINASE', 'BOTH_KINASE']:
        mini.ix[row,'3p-kinase'] = 1
    if frame == 'InFrame':
        mini.ix[row,'in-frame'] = 1
    if aa_seq.find('*') > -1:
        mini.ix[row,'premature stop codon'] = 1
    if transcript_5p.find('Intron') > -1:
        mini.ix[row,'5p-intron'] = 1
    if transcript_3p.find('Intron') > -1:
        mini.ix[row,'3p-intron'] = 1
    if BPR_5p.find('Intergenic') > -1:
        mini.ix[row,'5p-BPR intergenic'] = 1
    if BPR_5p.find('CDS') > -1:
        mini.ix[row,'5p-BPR CDS'] = 1
    if BPR_3p.find('Intergenic') > -1:
        mini.ix[row,'3p-BPR intergenic'] = 1
    if BPR_3p.find('CDS') > -1:
        mini.ix[row,'3p-BPR CDS'] = 1

domains = x_input[x_input.columns[39:]]
full = pd.merge(mini, domains, left_index=True, right_index=True)
X = full.values
f_log.write('successfully built feature space matrix with shape {0}\n'.format(X.shape))
y = clf.predict_proba(X)[:,1]
f_log.write('successfully applied classifier to feature space')
df1 = pd.DataFrame(y, index=x_input.index, columns=['DriverScore'])
df2 = x_input[x_input.columns[:39]]
df_final = pd.merge(df1, df2, left_index=True, right_index=True)
df_final.to_csv(report_out, sep='\t', index=False, float_format='%1.4f')
f_log.write('final report written to {0}\n'.format(report_out))
f_log.write('END')
f_log.close()
