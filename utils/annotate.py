import pandas as pd
import sys,os,io,re
import numpy as np

path=sys.argv[1]
event_file=sys.argv[2]
protName=sys.argv[3]

filt_events=pd.read_csv(event_file)
filelist=os.listdir(path)
effect_files=[]
for file in filelist:
 if file.startswith(protName) and file.endswith('.effects.csv'):
  effect_files.append(file)

colnames=pd.read_csv(effect_files[0],nrows=1).columns
select_col=colnames[range(np.where(colnames=='event_jid')[0][0],len(colnames))]

filt_effects=pd.DataFrame(columns=select_col)
for file in effect_files:
    curr_effects=pd.read_csv(file,usecols=select_col)
    filt_effects=filt_effects.append(curr_effects[curr_effects.event_jid.isin(filt_events.event_jid)])

filt_events=pd.merge(filt_events,filt_effects,on="event_jid",how="left")

top_file=protName + ".top.peptides.csv"
top_pept=pd.read_csv(top_file)
top_pept['event_jid']=top_pept.effectId.apply(lambda x: x.partition('_')[2])

filt_events=pd.merge(filt_events,top_pept,on="event_jid",how="left")
if not ('event_type' in filt_events.columns):
    filt_events['event_type']=filt_events.event_id_x.apply(lambda x: x.split('_')[0]+'_'+x.split('_')[1])

filt_events.insert(column='event_cat',loc=5,value=None)
filt_events.loc[filt_events.event_type=='alt_3prime','event_cat']='Alt'
filt_events.loc[filt_events.event_type=='alt_5prime','event_cat']='Alt'
filt_events.loc[filt_events.event_type=='mutex_exons','event_cat']='MutEx'
filt_events.loc[(filt_events.event_type=='exon_skip') & (filt_events.wtIsoform=='iso1'),'event_cat']='ExonInc'
filt_events.loc[(filt_events.event_type=='exon_skip') & (filt_events.wtIsoform=='iso2'),'event_cat']='ExonSkip'
filt_events.loc[(filt_events.event_type=='intron_retention') & (filt_events.wtIsoform=='iso2'),'event_cat']='IntronRet'
filt_events.loc[(filt_events.event_type=='intron_retention') & (filt_events.wtIsoform=='iso1'),'event_cat']='IntronExc'

if 'PSI_higher' in filt_events.columns:
    filt_events.insert(column='group_higher',loc=5,value=filt_events.PSI_higher.copy())
    filt_events.loc[(filt_events.wtIsoform=='iso1'),'group_higher']=filt_events.loc[(filt_events.wtIsoform=='iso1'),'PSI_lower']

filt_events.insert(column='aa_change_type',loc=7,value='Other')
filt_events.loc[filt_events.coding_transcript_effect.apply(lambda x: x.startswith('Novel'))  & filt_events.novelSeq & filt_events.novelPept,'aa_change_type']='Novel'
filt_events.loc[(filt_events.coding_transcript_effect=='IsoformSwitch')  & (filt_events.novelSeq==False) & (filt_events.novelPept==False),'aa_change_type']='Canonical'
filt_events.loc[filt_events.coding_transcript_effect=='Unknown','aa_change_type']='ND'
filt_events.loc[filt_events.coding_transcript_effect=='NonCoding','aa_change_type']='None'
filt_events.loc[filt_events.aa_effect=='Silent','aa_change_type']='None'
filt_events.loc[filt_events.aa_effect=='ProteinLoss','aa_change_type']='None'

filt_events.insert(column='effect_cat',loc=8,value="NA")
filt_events.loc[filt_events.orf_effect=='InFrame','effect_cat']=filt_events.loc[filt_events.orf_effect=='InFrame','aa_effect']
filt_events.loc[filt_events.orf_effect=='PrematureStop','effect_cat']='FrameDisruption'
filt_events.loc[filt_events.orf_effect=='PostStop','effect_cat']='FrameDisruption'
filt_events.loc[filt_events.orf_effect=='StopLost','effect_cat']='ProteinLoss'
filt_events.loc[filt_events.orf_effect=='StartLost','effect_cat']='ProteinLoss'
filt_events.loc[filt_events.coding_transcript_effect=='NonCoding','effect_cat']='NonCoding'

filt_events.to_csv(event_file.replace('.csv','.anno.csv'))
summary=filt_events.pivot_table(index=['event_cat','aa_change_type','effect_cat'],columns=['group_higher'],values=['effect_jid'],aggfunc=len,fill_value=0)
summary.to_csv(event_file.replace('.csv','.anno.summary.csv'))
