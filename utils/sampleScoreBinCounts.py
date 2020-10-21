import pandas as pd
import sys,os,io,re
import numpy as np
import itertools

path=sys.argv[1]
outName=sys.argv[2]
thresh=int(sys.argv[3])

bin_names=['low_out','not_out','high_out','no_score']
event_types=['alt_3prime','alt_5prime','exon_skip','intron_retention','mutex_exons']

cats=list(itertools.product(bin_names,event_types))

filelist=os.listdir(path)
out_files=[]
for file in filelist:
    if file.startswith(outName) and file.endswith('.bisbeeOutlier.csv'):
        out_files.append(file)

columns=pd.read_csv(path + "/" + out_files[0],nrows=1).columns
data_col=columns[np.where(columns=='event_jid')[0][0]+1:]
out_counts=pd.DataFrame(index=pd.MultiIndex.from_tuples(cats,names=['bins','event_type']),columns=data_col,data=0)

for file in out_files:
    curr_data=pd.read_csv(path + "/" + file)
    curr_data['event_type']=curr_data.event_id.apply(lambda x: x.split('_')[0]+'_'+x.split('_')[1])
    for col in data_col:
        bin_data=pd.cut(curr_data[col],[-np.inf,-thresh,thresh,np.inf],labels=bin_names[0:-1])
        bin_data=bin_data.cat.add_categories(bin_names[-1]).fillna(bin_names[-1])
        curr_counts=pd.crosstab(bin_data,curr_data.event_type).stack()
        out_counts.loc[curr_counts.index,col]+=curr_counts

#out_counts['total']=total_counts
out_counts.transpose().to_csv(outName + '.bisbeeOutlier.sampleScoreBinCounts.csv')
