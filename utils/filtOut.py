import pandas as pd
import sys,os,io,re
import numpy as np
import time

path=sys.argv[1]
outName=sys.argv[2]
thresh=float(sys.argv[3])
sample_count=int(sys.argv[4])
if len(sys.argv)>5:
    sample_file=sys.argv[5]
    select_group=sys.argv[6]
if len(sys.argv)>7:
    exclude_group=sys.argv[7]
    exclude_count=int(sys.argv[8])

filelist=os.listdir(path)
out_files=[]
for file in filelist:
 if file.startswith(outName) and file.endswith('.bisbeeOutlier.csv'):
  out_files.append(file)

colnames=pd.read_csv(path + "/" + out_files[0],nrows=1).columns
info_col=["contig","gene","strand","event_id","confirmed","event_jid"]
info_col=list(np.intersect1d(info_col,colnames))
#info_col=colnames.values[0:4]
#info_col=np.append(info_col,'event_jid')
try:
    sample_table=pd.read_table(sample_file,names=['samples','group'])
    sample_table.samples=sample_table.samples.apply(lambda x: x.replace('-','.'))
except:
    sample_names=colnames[range(np.where(colnames=='event_jid')[0][0]+1,len(colnames))]
    sample_table=pd.DataFrame({'samples':sample_names})
    sample_table['group']='test'
    select_group='test'

group_list=sample_table.group.unique()
print(group_list)
select_col=info_col
for group in group_list:
    select_col=np.append(select_col,[group + '_high_counts',group + '_low_counts',group + '_median_score',group + '_mean_score',group + '_abs_max_score'])
filt_events=pd.DataFrame(columns=select_col)
select_samples=list(sample_table.samples[sample_table.group==select_group])
filt_data=pd.DataFrame(columns=info_col + select_samples)
count_col=select_col[list(map(lambda x: x.endswith('counts'),select_col))]
for file in out_files:
    t0= time.time()
    print(file)
    curr_events=pd.read_csv(path + "/" + file,usecols=np.append(info_col,sample_table.samples.unique()))
    for group in group_list:
        sample_names=sample_table.samples[sample_table.group==group]
        curr_events[group + '_high_counts']=(curr_events[sample_names]>thresh).sum(1)
        curr_events[group + '_low_counts']=(curr_events[sample_names]<(-thresh)).sum(1)
        curr_events[group + '_median_score']=curr_events[sample_names].median(axis=1)
        curr_events[group + '_mean_score']=curr_events[sample_names].mean(axis=1)
        curr_events[group + '_abs_max_score']=curr_events[sample_names].abs().max(axis=1)
    sIdx=curr_events.loc[:,[select_group + '_high_counts',select_group + '_low_counts']].max(axis=1)>=sample_count
    try:
        idx=sIdx & (curr_events.loc[:,[exclude_group + '_high_counts',exclude_group + '_low_counts']].max(axis=1)<=exclude_count)
    except:
        idx=sIdx
    filt_events=filt_events.append(curr_events.loc[idx,select_col])
    filt_data=filt_data.append(curr_events.loc[idx,info_col + select_samples])
    t1 = time.time() - t0
    print("Time elapsed: ", t1 - t0,flush=True)

filt_events.insert(column=select_group + '_total_out',loc=5,
                   value=filt_events.loc[:,[select_group + '_high_counts',select_group + '_low_counts']].sum(axis=1))
filt_events=filt_events.sort_values(by=[select_group + "_total_out",select_group + "_abs_max_score"],ascending=False)
filt_events['event_type']=filt_events.event_id.apply(lambda x: x.split('_')[0]+'_'+x.split('_')[1])

filt_events.insert(column='PSI_higher',loc=6,value='Ref')
filt_events.insert(column='PSI_lower',loc=7,value='Ref')
filt_events.loc[filt_events[select_group + '_high_counts']>=sample_count,'PSI_higher']=select_group
filt_events.loc[filt_events[select_group + '_low_counts']>=sample_count,'PSI_lower']=select_group

name=outName + '.bisbeeOutlier.thresh' + str(thresh) + '.min' + str(sample_count) + select_group
try:
    name=name + '.max' + str(exclude_count) + exclude_group + '.csv'
except:
    name=name + '.csv'

filt_events.to_csv(name)
filt_data.to_csv(name.replace('.csv','.scores.csv'))
summary=filt_events.pivot_table(index=['event_type'],columns=['PSI_higher'],values=['event_jid'],aggfunc=len,fill_value=0)
summary.to_csv(name.replace('.csv','.summary.csv'))
