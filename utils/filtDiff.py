import pandas as pd
import sys,os,io,re
import numpy as np

path=sys.argv[1]
diffName=sys.argv[2]
thresh=float(sys.argv[3])

filelist=os.listdir(path)
diff_files=[]
for file in filelist:
 if file.startswith(diffName) and file.endswith('.bisbeeDiff.csv'):
  diff_files.append(file)

colnames=pd.read_csv(path + "/" + diff_files[0],nrows=1).columns
select_col=colnames[0:4]
select_col=select_col.append(colnames[range(np.where(colnames=='event_jid')[0][0],len(colnames))])
print(colnames)
print(select_col)

filt_events=pd.DataFrame(columns=select_col)
for file in diff_files:
    print(file)
    curr_events=pd.read_csv(path + "/" + file,usecols=select_col)
    if np.isnan(thresh):
        filt_events=filt_events.append(curr_events)
    else:
        filt_events=filt_events.append(curr_events[curr_events.ll_ratio>thresh])

filt_events=filt_events.sort_values(by="ll_ratio",ascending=False)
filt_events['event_type']=filt_events.event_id.apply(lambda x: x.split('_')[0]+'_'+x.split('_')[1])

group_names=filt_events.columns[list(map(lambda x: x.startswith('fitPSI'),filt_events.columns))][0:2]
filt_events.insert(column='PSI_higher',loc=5,value=None)
filt_events.insert(column='PSI_lower',loc=6,value=None)
filt_events.loc[filt_events[group_names[0]]>filt_events[group_names[1]],'PSI_higher']=group_names[0].replace('fitPSI_','')
filt_events.loc[filt_events[group_names[0]]<filt_events[group_names[1]],'PSI_higher']=group_names[1].replace('fitPSI_','')
filt_events.loc[filt_events[group_names[0]]>filt_events[group_names[1]],'PSI_lower']=group_names[1].replace('fitPSI_','')
filt_events.loc[filt_events[group_names[0]]<filt_events[group_names[1]],'PSI_lower']=group_names[0].replace('fitPSI_','')

filt_events.to_csv(diffName + '.bisbeeDiff.thresh' + str(thresh) + '.csv')
summary=filt_events.pivot_table(index=['event_type'],columns=['PSI_higher'],values=['event_jid'],aggfunc=len,fill_value=0)
summary.to_csv(diffName + '.bisbeeDiff.thresh' + str(thresh) + '.summary.csv')
