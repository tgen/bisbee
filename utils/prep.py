import pandas as pd
import sys, os, re
import h5py
import numpy as np

events_file=sys.argv[1]
event_type=sys.argv[2]
outname=sys.argv[3]

event_dict={'A3':'alt_3prime','A5':'alt_5prime','ES':'exon_skip','IR':'intron_retention','MUT':'mutex_exons'}

f=h5py.File(events_file,'r')
gene_idx=f['gene_idx']
conf=f['confirmed']
gene_names=np.array(list(map(lambda x: x.decode('UTF-8'),f['gene_names'])))
gene_chr=np.array(list(map(lambda x: x.decode('UTF-8'),f['gene_chr'])))
gene_strand=np.array(list(map(lambda x: x.decode('UTF-8'),f['gene_strand'])))
sample_names=list(map(lambda x: x.decode('UTF-8'),f['strains']))

data=pd.DataFrame({'gene':gene_names[gene_idx],'strand':gene_strand[gene_idx],'contig':gene_chr[gene_idx]})
data['event_id']=list(map(lambda x: event_dict[event_type] + '_' + str(x),range(1,len(gene_idx)+1)))
data['confirmed']=conf

event_pos=f['event_pos']
if event_type=='A3':
    data['exon_const_start']=event_pos[:,0]+1
    data['exon_const_end']=event_pos[:,1]
    data['exon_alt1_start']=event_pos[:,4]+1
    data['exon_alt1_end']=event_pos[:,5]
    data['exon_alt2_start']=event_pos[:,6]+1
    data['exon_alt2_end']=event_pos[:,7]
    data.loc[data['strand']=='-','exon_const_start']=event_pos[data['strand']=='-',4]+1
    data.loc[data['strand']=='-','exon_const_end']=event_pos[data['strand']=='-',5]
    data.loc[data['strand']=='-','exon_alt1_start']=event_pos[data['strand']=='-',0]+1
    data.loc[data['strand']=='-','exon_alt1_end']=event_pos[data['strand']=='-',1]
    data.loc[data['strand']=='-','exon_alt2_start']=event_pos[data['strand']=='-',2]+1
    data.loc[data['strand']=='-','exon_alt2_end']=event_pos[data['strand']=='-',3]
    iso1_str=data.apply(lambda x: str(x['exon_const_end'])+'j'+str(x['exon_alt2_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_const_end'])+'j'+str(x['exon_alt1_start']),axis=1)
    iso1_str[data['strand']=='-']=data[data['strand']=='-'].apply(lambda x: str(x['exon_alt2_end'])+'j'+str(x['exon_const_start']),axis=1)
    iso2_str[data['strand']=='-']=data[data['strand']=='-'].apply(lambda x: str(x['exon_alt1_end'])+'j'+str(x['exon_const_start']),axis=1)

if event_type=='A5':
    data['exon_const_start']=event_pos[:,0]+1
    data['exon_const_end']=event_pos[:,1]
    data['exon_alt1_start']=event_pos[:,4]+1
    data['exon_alt1_end']=event_pos[:,5]
    data['exon_alt2_start']=event_pos[:,6]+1
    data['exon_alt2_end']=event_pos[:,7]
    data.loc[data['strand']=='+','exon_const_start']=event_pos[data['strand']=='+',4]+1
    data.loc[data['strand']=='+','exon_const_end']=event_pos[data['strand']=='+',5]
    data.loc[data['strand']=='+','exon_alt1_start']=event_pos[data['strand']=='+',0]+1
    data.loc[data['strand']=='+','exon_alt1_end']=event_pos[data['strand']=='+',1]
    data.loc[data['strand']=='+','exon_alt2_start']=event_pos[data['strand']=='+',2]+1
    data.loc[data['strand']=='+','exon_alt2_end']=event_pos[data['strand']=='+',3]
    iso1_str=data.apply(lambda x: str(x['exon_const_end'])+'j'+str(x['exon_alt2_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_const_end'])+'j'+str(x['exon_alt1_start']),axis=1)
    iso1_str[data['strand']=='+']=data[data['strand']=='+'].apply(lambda x: str(x['exon_alt2_end'])+'j'+str(x['exon_const_start']),axis=1)
    iso2_str[data['strand']=='+']=data[data['strand']=='+'].apply(lambda x: str(x['exon_alt1_end'])+'j'+str(x['exon_const_start']),axis=1)

if event_type=='ES':
    data['exon_pre_end']=event_pos[:,1]
    data['exon_start']=event_pos[:,2]+1
    data['exon_end']=event_pos[:,3]
    data['exon_aft_start']=event_pos[:,4]+1
    iso1_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon_start'])+'_'+str(x['exon_end'])+'j'+str(x['exon_aft_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon_aft_start']),axis=1)

if event_type=='IR':
    data['exon1_end']=event_pos[:,1]
    data['exon2_start']=event_pos[:,2]+1
    iso1_str='NONE'
    iso2_str=data.apply(lambda x: str(x['exon1_end'])+'j'+str(x['exon2_start']),axis=1)

if event_type=='MUT':
    event_coord=np.reshape(event_pos,(len(gene_idx),8),(0,2,1))
    data['exon_pre_end']=event_coord[:,1]
    data['exon1_start']=event_coord[:,2]+1
    data['exon1_end']=event_coord[:,3]
    data['exon2_start']=event_coord[:,4]+1
    data['exon2_end']=event_coord[:,5]
    data['exon_aft_start']=event_coord[:,6]+1
    iso1_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon1_start'])+'_'+str(x['exon1_end'])+'j'+str(x['exon_aft_start']),axis=1)
    iso2_str=data.apply(lambda x: str(x['exon_pre_end'])+'j'+str(x['exon2_start'])+'_'+str(x['exon2_end'])+'j'+str(x['exon_aft_start']),axis=1)

data['event_jid']=data['contig'] + ':g.' + iso1_str + '>' + iso2_str + "[spl" + event_type + "]"
iso1=pd.DataFrame(np.transpose(f['iso1']),columns=list(map(lambda x: x + '_iso1',sample_names)))
iso2=pd.DataFrame(np.transpose(f['iso2']),columns=list(map(lambda x: x + '_iso2',sample_names)))
data=pd.concat([data,iso1,iso2],axis=1)
data.to_csv(outname,index=False)
