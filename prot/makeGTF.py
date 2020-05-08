import pyensembl
import sys,os
import pandas as pd
import utils as bb

event_file=sys.argv[1]
ensemble_release=int(sys.argv[2])

events_table=pd.read_csv(event_file)
ensembl=pyensembl.EnsemblRelease(ensemble_release)
gtf=pd.DataFrame(columns=["seqname","feature","start","end","strand","attribute"])

for index,row in events_table.iterrows():
    event_coords=bb.jid_to_coords(events_table.loc[index,"event_jid"])
    try:
        tid=row["effectId"].split('_')[0]
    except:
        continue
    transcript=ensembl.transcript_by_id(tid)
    isoform=bb.get_matching_isoform(transcript,event_coords)
    attribute='gene_id ' + transcript.gene_id + '; transcript_id ' + row["effectId"]
    exon_coord=pd.DataFrame(transcript.exon_intervals,columns=["start","end"])
    new_coord=bb.get_new_coord(isoform,event_coords,exon_coord)
    gtf=gtf.append({"seqname":transcript.contig,"feature":"transcript","start":new_coord.start.min(),"end":new_coord.end.max(),"strand":transcript.strand,"attribute":attribute},ignore_index=True)
    new_coord["seqname"]=transcript.contig
    new_coord["feature"]="exon"
    new_coord["strand"]=transcript.strand
    new_coord["attribute"]=list(map(lambda x:'gene_id ' + transcript.gene_id + '; transcript_id ' + row["effectId"] + "; exon_number " + str(x),range(1,new_coord.shape[0]+1)))
    gtf=gtf.append(new_coord)

gtf["source"]=event_file
gtf["score"]=0
gtf["frame"]="."

gtf=gtf[["seqname","source","feature","start","end","score","strand","frame","attribute"]]
gtf.to_csv(event_file.replace("csv","gtf"),sep="\t",header=False,index=False)
