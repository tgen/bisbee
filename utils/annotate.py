import pandas as pd
import sys,os,io,re
import numpy as np

event_file=sys.argv[1]
path=sys.argv[2]
protName=sys.argv[3]

filt_events=pd.read_csv(event_file)
filelist=os.listdir(path)
effect_files=[]
for file in filelist:
 if file.startswith(protName + '.') and file.endswith('.effects.csv'):
  effect_files.append(file)

colnames=pd.read_csv(path + "/" + effect_files[0],nrows=1).columns
select_col=colnames[range(np.where(colnames=='event_jid')[0][0],len(colnames))]

filt_effects=pd.DataFrame(columns=select_col)
for file in effect_files:
    curr_effects=pd.read_csv(path + "/" + file,usecols=select_col)
    filt_effects=filt_effects.append(curr_effects[curr_effects.event_jid.isin(filt_events.event_jid)])

filt_effects=filt_effects.drop_duplicates()
filt_events=pd.merge(filt_events,filt_effects,on="event_jid",how="left")

top_file=path + "/" + protName + ".top.peptides.csv"
top_pept=pd.read_csv(top_file)
#top_pept['event_jid']=top_pept.effectId.apply(lambda x: x.partition('_')[2])

filt_events=pd.merge(filt_events,top_pept,on="event_jid",how="left")
if not ('event_type' in filt_events.columns):
    filt_events['event_type']=filt_events.event_id_x.apply(lambda x: x.split('_')[0]+'_'+x.split('_')[1])

filt_events.insert(column='event_cat',loc=5,value=None)
filt_events.loc[filt_events.event_type=='alt_3prime','event_cat']='Alt3'
filt_events.loc[filt_events.event_type=='alt_5prime','event_cat']='Alt5'
filt_events.loc[filt_events.event_type=='mutex_exons','event_cat']='MutEx'
filt_events.loc[filt_events.event_type=='exon_skip','event_cat']=filt_events.loc[filt_events.event_type=='exon_skip','PSI_higher'].apply(lambda x: str(x) + 'ExonInc')
filt_events.loc[filt_events.event_type=='intron_retention','event_cat']=filt_events.loc[filt_events.event_type=='intron_retention','PSI_higher'].apply(lambda x: str(x) + 'IntronRet')

if 'PSI_higher' in filt_events.columns:
    filt_events.insert(column='group_increased_alt',loc=5,value=filt_events.PSI_higher.copy())
    filt_events.loc[(filt_events.refIsoform=='iso1'),'group_increased_alt']=filt_events.loc[(filt_events.refIsoform=='iso1'),'PSI_lower']

filt_events.insert(column='ref_seq_header',loc=7,value=None)
filt_events.loc[(filt_events.refIsoform=='iso1'),'ref_seq_header']=filt_events.loc[(filt_events.refIsoform=='iso1'),'iso1_header']
filt_events.loc[(filt_events.refIsoform=='iso2'),'ref_seq_header']=filt_events.loc[(filt_events.refIsoform=='iso2'),'iso2_header']
filt_events.insert(column='alt_seq_header',loc=8,value=None)
filt_events.loc[(filt_events.refIsoform=='iso1'),'alt_seq_header']=filt_events.loc[(filt_events.refIsoform=='iso1'),'iso2_header']
filt_events.loc[(filt_events.refIsoform=='iso2'),'alt_seq_header']=filt_events.loc[(filt_events.refIsoform=='iso2'),'iso1_header']

print(filt_events.head())
filt_events.insert(column='aa_change_type',loc=9,value='Other')
filt_events.loc[filt_events.coding_transcript_effect.apply(lambda x: x.startswith('Novel'))  & filt_events.novelSeq & filt_events.novelPept,'aa_change_type']='Novel'
filt_events.loc[(filt_events.coding_transcript_effect=='IsoformSwitch')  & (filt_events.novelSeq==False) & (filt_events.novelPept==False),'aa_change_type']='Canonical'
filt_events.loc[filt_events.coding_transcript_effect=='Unknown','aa_change_type']='ND'
filt_events.loc[filt_events.coding_transcript_effect=='NonCoding','aa_change_type']='None'
filt_events.loc[filt_events.aa_effect=='Silent','aa_change_type']='None'
filt_events.loc[filt_events.aa_effect=='ProteinLoss','aa_change_type']='None'

filt_events.insert(column='effect_cat',loc=10,value="NA")
filt_events.loc[filt_events.orf_effect=='InFrame','effect_cat']=filt_events.loc[filt_events.orf_effect=='InFrame','aa_effect']
filt_events.loc[filt_events.orf_effect=='PrematureStop','effect_cat']='FrameDisruption'
filt_events.loc[filt_events.orf_effect=='PostStop','effect_cat']='FrameDisruption'
filt_events.loc[filt_events.orf_effect=='StopLost','effect_cat']='ProteinLoss'
filt_events.loc[filt_events.orf_effect=='StartLost','effect_cat']='ProteinLoss'
filt_events.loc[filt_events.coding_transcript_effect=='NonCoding','effect_cat']='NonCoding'

filt_events.to_csv(event_file.replace('.csv','.anno.csv'))
summary=filt_events.fillna('NA').pivot_table(index=['event_cat','aa_change_type','effect_cat'],columns=['group_increased_alt'],values=['event_jid'],aggfunc=len,fill_value=0)
summary.to_csv(event_file.replace('.csv','.anno.summary.csv'))
