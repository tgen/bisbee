import Bio.SeqIO
import pandas as pd
import sys,os,io,re
import csv

path=sys.argv[1]
name=sys.argv[2]

pept_file=path + "/" +name + ".all.peptides.csv"
if os.path.exists(pept_file):
 os.remove(pept_file)
top_file=path + "/" + name + ".top.peptides.csv"
if os.path.exists(top_file):
 os.remove(top_file)

filelist=os.listdir(path)
peptfiles=[]
for file in filelist:
 if file.startswith(name) and file.endswith('.peptides.csv'):
  peptfiles.append(file)

seq_table=pd.read_csv(path + "/" + name + ".seqTable.csv")

peptTable=pd.concat([pd.read_csv(path + "/" + f) for f in peptfiles],ignore_index=True)
peptTable["novelPept"]=peptTable.altPept.isin(peptTable.refPept)==False
enst=peptTable.effectId.apply(lambda x: x.split('_')[0])
junc=peptTable.effectId.apply(lambda x: x.split('g.')[1].split('[')[0])
iso1_str=junc.apply(lambda x: x.split('>')[0])
iso2_str=junc.apply(lambda x: x.split('>')[1])
iso1_id="en|" + enst + ":g." + junc
iso2_id="en|" + enst + ":g." + iso2_str + ">" + iso1_str
peptTable.loc[peptTable.refIsoform=="iso1","iso1_id"]="en|" + enst[peptTable.refIsoform=="iso1"]
peptTable.loc[peptTable.refIsoform=="iso2","iso2_id"]="en|" + enst[peptTable.refIsoform=="iso2"]
peptTable.loc[peptTable.refIsoform=="iso1","iso2_id"]=iso1_id[peptTable.refIsoform=="iso1"]
peptTable.loc[peptTable.refIsoform=="iso2","iso1_id"]=iso2_id[peptTable.refIsoform=="iso2"]

peptTable=peptTable.merge(seq_table,left_on="iso1_id",right_on="seq_id",how="left")
peptTable.rename(columns={"seq_header":"iso1_header"},inplace=True)
peptTable.drop(["seq_id"],inplace=True,axis=1)
peptTable=peptTable.merge(seq_table,left_on="iso2_id",right_on="seq_id",how="left")
peptTable.rename(columns={"seq_header":"iso2_header"},inplace=True)
peptTable["novelSeq"]=peptTable.novelSeq_x | peptTable.novelSeq_y
peptTable.drop(["seq_id","novelSeq_x","novelSeq_y"],inplace=True,axis=1)

peptTable["event_jid"]=peptTable.effectId.apply(lambda x: x.partition('_')[2])
peptTable.sort_values(by=["event_jid","novelPept","novelSeq","insSeqLen","refSeqLen"],ascending=[True,False,True,False,False],inplace=True)
topPeptTable=peptTable.drop_duplicates(subset="event_jid")
peptTable.topEffect=peptTable.effectId.isin(topPeptTable.effectId)
topPeptTable.topEffect=True

peptTable.to_csv(pept_file)
topPeptTable.to_csv(top_file)

fasta=path + "/" + name + ".unique.fasta"
records=Bio.SeqIO.to_dict(Bio.SeqIO.parse(fasta,format="fasta"))
keys=set(records.keys())
iso1=keys.intersection(topPeptTable.iso1_header)
iso2=keys.intersection(topPeptTable.iso2_header)
top=iso1.union(iso2)

top_dict=dict()
for seq_id in top:
 top_dict[seq_id]=records[seq_id]

Bio.SeqIO.write(top_dict.values(),name + ".top.fasta","fasta")
