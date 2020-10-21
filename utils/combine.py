import pandas as pd
import sys,os,io,re
import numpy as np

path=sys.argv[1]
prefix=sys.argv[2]
suffix=sys.argv[3]
data_col=sys.argv[4].split(',')

anno_col=["event_cat","group_increased_alt","ref_seq_header","alt_seq_header",
          "aa_change_type","effect_cat"]

filelist=os.listdir(path)
data_files={}
for file in filelist:
 if file.startswith(prefix) and file.endswith(suffix):
  data_files[file.replace(prefix,"").replace(suffix,"")]=file

name,file=data_files.popitem()
data_col_names=list(map(lambda x: x+"_"+name,data_col))
data=pd.read_csv(file,index_col="event_jid")
data=data.rename(columns=dict(zip(data_col,data_col_names)))

for name,file in data_files.items():
    new_data=pd.read_csv(file,index_col="event_jid")
    new_data_col_names=list(map(lambda x: x+"_"+name,data_col))
    data_col_names=data_col_names + new_data_col_names
    new_data=new_data.rename(columns=dict(zip(data_col,new_data_col_names)))
    data=data.combine_first(new_data)

anno_col_names=np.intersect1d(anno_col,data.columns)
other_col=np.setdiff1d(data.columns,np.append(anno_col_names,data_col_names))

col_ord=list(np.append(np.append(anno_col_names,data_col_names),other_col))
print(col_ord)
data.loc[:,col_ord].to_csv(prefix + 'comb' + suffix)
