import Bio.SeqIO
import pandas as pd
import sys,os,io,re

path=sys.argv[1]
name=sys.argv[2]

out_file=path + "/" + '.'.join([name,"unique.fasta"])
if os.path.exists(out_file):
 os.remove(out_file)

out_table=path + "/" + '.'.join([name,"seqTable.csv"])
if os.path.exists(out_table):
 os.remove(out_table)

filelist=os.listdir(path)
fastalist=[]
for file in filelist:
 if file.startswith(name) and file.endswith('.fasta'):
  fastalist.append(file)

unique_seq=dict()

for fasta in fastalist:
 for record in Bio.SeqIO.parse(path + "/" + fasta,"fasta"):
  header=record.description
  seq=str(record.seq)
  fields=header.split(' ')
  seq_id=fields[0]
  gene=re.sub(string=fields[1],pattern='GN=',repl='')
  event_id=re.sub(string=fields[2],pattern='ID=',repl='')
  jid=re.sub(string=fields[3],pattern='JID=',repl='')
  pc=re.sub(string=fields[4],pattern='PC=',repl='')
  #event_info=event_id + "_" + pc
  if seq in unique_seq.keys():
   info=unique_seq[seq]
   info["ids"].add(seq_id)
   info["genes"].add(gene)
   info["events"][event_id]=pc
   info["jids"].add(jid)
   unique_seq[seq]=info
  else:
   info={"ids":{seq_id},"genes":{gene},"events":{event_id: pc},"jids":{jid}}
   unique_seq[seq]=info

out_fasta=open(out_file,"a+")
out_csv=open(out_table,"a+")
seq_table=pd.DataFrame(columns=["seq_header","seq_id","novelSeq"])
seq_table.to_csv(out_csv,header=True)
#i=0
for seq,info in unique_seq.items():
 if len(seq)==0:
  continue
 wt_ids=[]
 for seq_id in info["ids"]:
  if seq_id.find('g')==-1:
   wt_ids.append(seq_id)
 if len(wt_ids)>0:
  novel=False
  id_str='en|' + re.sub(string='|'.join(wt_ids),pattern='en\|',repl='')
 else:
  novel=True
  id_str=id_str='en|' + re.sub(string='|'.join(info["ids"]),pattern='en\|',repl='')
 seq_table=pd.DataFrame({"seq_header": id_str, "seq_id": list(info["ids"]),"novelSeq": novel})
 seq_table=seq_table[["seq_header","seq_id","novelSeq"]]
 genes='|'.join(info["genes"])
 events=''
 for event_id,pc in info["events"].items():
  events+=event_id + "_" + pc +'|'
 jids='|'.join(info["jids"])
 desc='GN='+ genes + " ID=" + events + " JID=" + jids + " NS=" + str(novel)
 out_handle=io.StringIO()
 record=Bio.SeqRecord.SeqRecord(id=id_str,description=desc,seq=Bio.Seq.Seq(seq))
 Bio.SeqIO.write(record,out_handle,"fasta")
 out_fasta.write(out_handle.getvalue())
 seq_table.to_csv(out_csv,header=False)

out_fasta.close()
out_csv.close()
