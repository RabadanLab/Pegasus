from datetime import datetime
import cPickle
import numpy as np
import pandas as pd
from sklearn import ensemble

print '------------------------------'
print '[ loading training data ] {}'.format(datetime.utcnow())

x_train_pos = pd.read_csv('../resources/data_training/X_chimerDB.txt', header=0, sep='\t')
x_train_neg_chimera = pd.read_csv('../resources/data_training/X_chimerDB_frameshifted.txt', header=0, sep='\t')
neg_first = 0
neg_last = 1500
x_train_neg_RLN = pd.read_csv('../resources/data_training/X_reactiveLN.txt', header=0, sep='\t')
x_train_neg = pd.concat([x_train_neg_chimera, x_train_neg_RLN.ix[neg_first:neg_last,:]], axis=0, ignore_index=True)

print '[ {} positive training examples, {} negative training examples ] {}'.format(x_train_pos.shape[0], x_train_neg.shape[0], datetime.utcnow())

features = ['5p-kinase', '3p-kinase', 'in-frame', 'premature stop codon', 
            '5p-intron', '3p-intron', '5p-BPR intergenic', '5p-BPR CDS', '3p-BPR intergenic', '3p-BPR CDS']
mini_train_pos = pd.DataFrame(np.zeros((len(x_train_pos),len(features))), columns=features)
mini_train_neg = pd.DataFrame(np.zeros((len(x_train_neg),len(features))), columns=features)

for row in x_train_pos.index:
    cols = list(x_train_pos.ix[row,:])
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
        mini_train_pos.ix[row,'5p-kinase'] = 1
    if kinase in ['3p_KINASE', 'BOTH_KINASE']:
        mini_train_pos.ix[row,'3p-kinase'] = 1
    if frame == 'InFrame':
        mini_train_pos.ix[row,'in-frame'] = 1
    if aa_seq.find('*') > -1:
        mini_train_pos.ix[row,'premature stop codon'] = 1
    if transcript_5p.find('Intron') > -1:
        mini_train_pos.ix[row,'5p-intron'] = 1
    if transcript_3p.find('Intron') > -1:
        mini_train_pos.ix[row,'3p-intron'] = 1
    if BPR_5p.find('Intergenic') > -1:
        mini_train_pos.ix[row,'5p-BPR intergenic'] = 1
    if BPR_5p.find('CDS') > -1:
        mini_train_pos.ix[row,'5p-BPR CDS'] = 1
    if BPR_3p.find('Intergenic') > -1:
        mini_train_pos.ix[row,'3p-BPR intergenic'] = 1
    if BPR_3p.find('CDS') > -1:
        mini_train_pos.ix[row,'3p-BPR CDS'] = 1

for row in x_train_neg.index:
    cols = list(x_train_neg.ix[row,:])
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
        mini_train_neg.ix[row,'5p-kinase'] = 1
    if kinase in ['3p_KINASE', 'BOTH_KINASE']:
        mini_train_neg.ix[row,'3p-kinase'] = 1
    if frame == 'InFrame':
        mini_train_neg.ix[row,'in-frame'] = 1
    if aa_seq.find('*') > -1:
        mini_train_neg.ix[row,'premature stop codon'] = 1
    if transcript_5p.find('Intron') > -1:
        mini_train_neg.ix[row,'5p-intron'] = 1
    if transcript_3p.find('Intron') > -1:
        mini_train_neg.ix[row,'3p-intron'] = 1
    if BPR_5p.find('Intergenic') > -1:
        mini_train_neg.ix[row,'5p-BPR intergenic'] = 1
    if BPR_5p.find('CDS') > -1:
        mini_train_neg.ix[row,'5p-BPR CDS'] = 1
    if BPR_3p.find('Intergenic') > -1:
        mini_train_neg.ix[row,'3p-BPR intergenic'] = 1
    if BPR_3p.find('CDS') > -1:
        mini_train_neg.ix[row,'3p-BPR CDS'] = 1

domains = x_train_pos.columns[39:]
domain_train_pos = x_train_pos[domains]
domain_train_neg = x_train_neg[domains]
full_train_pos = pd.merge(mini_train_pos, domain_train_pos, left_index=True, right_index=True)
full_train_neg = pd.merge(mini_train_neg, domain_train_neg, left_index=True, right_index=True)
X_train = np.concatenate((full_train_pos.values, full_train_neg.values), axis=0)
y_train = np.squeeze(np.concatenate((np.ones((len(full_train_pos),1)), np.zeros((len(full_train_neg),1)))))

print '[ fitting classifiers ] {}'.format(datetime.utcnow())
#svc = svm.SVC(C=1.0, kernel='rbf', gamma='auto', probability=True).fit(X_train, y_train)
gbc = ensemble.GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=5).fit(X_train, y_train)
rfc = ensemble.RandomForestClassifier(n_estimators=100).fit(X_train, y_train)

print '[ serializing fitted models ] {}'.format(datetime.utcnow())
#with open('models/trained_model_svc.pk','wb') as fh_svc:
#    cPickle.dump(svc, fh_svc)
with open('models/trained_model_gbc.pk','wb') as fh_gbc:
    cPickle.dump(gbc, fh_gbc)
with open('models/trained_model_rfc.pk','wb') as fh_rfc:
    cPickle.dump(rfc, fh_rfc)

print '------------------------------'
