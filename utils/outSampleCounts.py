import pandas as pd
import sys,os,io,re
import numpy as np

path=sys.argv[1]
outName=sys.argv[2]
thresh=int(sys.argv[3])
anno_file=sys.argv[4]

anno_table=pd.read_csv(anno_file)
anno_col=["event_cat","group_increased_alt","aa_change_type","effect_cat"]
anno_col=list(np.intersect1d(anno_col,anno_table.columns))

filelist=os.listdir(path)
out_files=[]
for file in filelist:
    if file.startswith(outName) and file.endswith('.bisbeeOutlier.csv'):
        out_files.append(file)

columns=pd.read_csv(path + "/" + out_files[0],nrows=1).columns
data_col=columns[np.where(columns=='event_jid')[0][0]+1:]
total_counts=anno_table.pivot_table(index=anno_col,values='event_jid',aggfunc=len)
total_counts.index = ["_".join(v) for v in total_counts.index.values]
out_counts=pd.DataFrame(index=total_counts.index,columns=data_col,data=0)
score_counts=pd.DataFrame(index=total_counts.index,columns=data_col,data=0)

for file in out_files:
    curr_data=pd.read_csv(path + "/" + file)
    curr_data=pd.merge(anno_table,curr_data,on='event_jid',how='inner')
    curr_data[anno_col].fillna('NA')
    curr_counts=pd.DataFrame(curr_data.pivot_table(index=anno_col,values=data_col,aggfunc=lambda x: sum(abs(x)>thresh)))
    curr_counts.index = ["_".join(v) for v in curr_counts.index.values]
    out_counts=out_counts.add(curr_counts,fill_value=0)
    curr_counts=pd.DataFrame(curr_data.pivot_table(index=anno_col,values=data_col,aggfunc=lambda x: sum(np.isnan(x)==False)))
    curr_counts.index = ["_".join(v) for v in curr_counts.index.values]
    score_counts=score_counts.add(curr_counts,fill_value=0)


#out_counts['total']=total_counts
out_counts.transpose().to_csv(anno_file.replace('anno.csv','samplePassCounts.csv'))
score_counts.transpose().to_csv(anno_file.replace('anno.csv','sampleScoreCounts.csv'))
