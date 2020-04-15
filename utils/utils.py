
import pyensembl
import Bio.SeqIO
import Bio.Seq
import pandas as pd
import sys
import re
from Bio import pairwise2


def get_transcript_adj_exons(ensembl,gene_id,exon_coord):
 try:
  transcript_ids=ensembl.transcript_ids_of_gene_id(gene_id)
 except:
  print('Warning: ' + gene_id + ' not found')
  transcript_ids=[]
 transcript_list=[]
 for tid in transcript_ids:
   transcript=ensembl.transcript_by_id(tid)
   transcript.exon_intervals.sort()
   if exon_coord[0] in transcript.exon_intervals:
    idx=transcript.exon_intervals.index(exon_coord[0])
    if exon_coord == transcript.exon_intervals[idx:idx+len(exon_coord)]:
      transcript_list.append(transcript)
 return transcript_list

def has_coding_transcript(transcript_list):
 has_coding=False
 for transcript in transcript_list:
  if transcript.biotype=='protein_coding':
   has_coding=True
 return has_coding

def get_transcript_contain_exons(ensembl,gene_id,exon_coord):
 try:
  transcript_ids=ensembl.transcript_ids_of_gene_id(gene_id)
 except:
  print('Warning: ' + gene_id + ' not found')
  transcript_ids=[]
 transcript_list=[]
 for tid in transcript_ids:
   transcript=ensembl.transcript_by_id(tid)
   if set(exon_coord).issubset(set(transcript.exon_intervals)):
     transcript_list.append(transcript)
 return transcript_list

def find_overlap(rangeDF,coord):
 return (coord[0]<=rangeDF.loc[:,'end']) & (coord[1]>=rangeDF.loc[:,'start'])

def make_seq_from_coord(ref,contig,coordDF,strand):
  seq=''
  if not contig in ref:
   contig="chr"+str(contig)
   if not contig in ref:
    print(contig + "not found in ref")
    return seq
  for index,row in coordDF.iterrows():
   if strand=='+':
    seq=seq+str(ref[contig].seq[int(row.start)-1:int(row.end)])
   else:
    seq=str(ref[contig].seq[int(row.start)-1:int(row.end)].reverse_complement())+seq
  return seq

def find_seq_diff(seq1,seq2):
 align_list=pairwise2.align.globalms(seq1,seq2,2,-1,-10,0)
 if len(align_list)==0:
  seq1_diff_pos=(0,len(seq1)-1)
  seq2_diff_pos=(0,len(seq2)-1)
 elif seq1==seq2:
  seq1_diff_pos=(float('nan'),float('nan'))
  seq2_diff_pos=(float('nan'),float('nan'))
 else:
  align_data=pairwise2.format_alignment(*align_list[0]).split('\n')
  #print(align_data)
  seq_comp=pd.DataFrame({"seq1": list(align_data[0]), "seq2": list(align_data[2])})
  seq_comp=seq_comp.assign(seq1_pos= seq_comp.seq1.isin(list(Bio.Seq.Alphabet.IUPAC.IUPACProtein.letters)).cumsum())
  seq_comp=seq_comp.assign(seq2_pos= seq_comp.seq2.isin(list(Bio.Seq.Alphabet.IUPAC.IUPACProtein.letters)).cumsum())
  seq_comp=seq_comp.assign(match=seq_comp.seq1==seq_comp.seq2)
 #print(seq_comp.head())
 #print(seq_comp[seq_comp.match==False])
  first_mismatch=seq_comp.loc[seq_comp.match==False].index[0]
  if first_mismatch==0:
    seq1_diff_pos=(-1,max(seq_comp.seq1_pos[seq_comp.match==False]))
    seq2_diff_pos=(-1,max(seq_comp.seq2_pos[seq_comp.match==False]))
  else:
   seq1_diff_pos=(seq_comp.seq1_pos[first_mismatch-1],max(seq_comp.seq1_pos[seq_comp.match==False]))
   seq2_diff_pos=(seq_comp.seq2_pos[first_mismatch-1],max(seq_comp.seq2_pos[seq_comp.match==False]))

 return seq1_diff_pos, seq2_diff_pos

def is_coding_effect(transcript,effect_coord):
 coding_effect=dict()
 if transcript.biotype!='protein_coding' or not transcript.contains_stop_codon or not transcript.contains_start_codon:
  coding_effect['type']='NoncodingOrIncompleteTranscript'
  coding_effect['is_coding']=False
 else:
  coding_coord=pd.DataFrame(transcript.coding_sequence_position_ranges,columns=['start','end'])
  coding_effect['coord']=coding_coord.append(pd.Series({'start':min(transcript.stop_codon_positions),'end':max(transcript.stop_codon_positions)}),ignore_index=True).sort_values(by="start")
  coding_effect['overlap']=find_overlap(coding_effect['coord'],effect_coord)
  if effect_coord[1]<coding_coord.start.min() or effect_coord[0]>coding_coord.end.max():
   coding_effect['type']='UTR'
   coding_effect['is_coding']=False
  else:
   coding_effect['is_coding']=True
 return coding_effect

def make_prot_from_coord(transcript,coord,ref):
  trans=Bio.Seq.translate(make_seq_from_coord(ref,transcript.contig,coord,transcript.strand))
  stop_count=trans.count('*')
  if trans.startswith('M'):
   start_lost=False
  else:
   start_lost=True
  if stop_count==0:
   stop_lost=True
  else:
   stop_lost=False

  if stop_lost and not start_lost:
   if transcript.strand=='+':
    coord.iloc[-1,1]=transcript.exon_intervals[-1][1]
   else:
    coord.iloc[0,0]=transcript.exon_intervals[0][0]
    trans=Bio.Seq.translate(make_seq_from_coord(ref,transcript.contig,coord,transcript.strand))
    stop_count=trans.count('*')

  if start_lost or stop_count==0:
   prot_seq=''
  else:
   prot_seq=trans.split('*')[0]+'*'

  if start_lost:
   effect='StartLost'
  elif stop_count==0:
   effect='StopLost'
  elif stop_lost:
   effect='PostStop'
  elif stop_count==1 and trans.endswith('*'):
   effect='InFrame'
  else:
   effect='PrematureStop'
  return prot_seq, effect

def get_event_coords(event_info,event_type,spladder_version):
 event_coords=pd.DataFrame(columns=("isoform","start","end"))
 if event_type=="IR":
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon1_end"],"end":event_info["exon2_start"]},ignore_index=True)
 elif event_type=="ES":
  event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_pre_end"],"end":event_info["exon_aft_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_pre_end"],"end":event_info["exon_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_end"],"end":event_info["exon_aft_start"]},ignore_index=True)
 elif event_type=="MUT":
  event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_pre_end"],"end":event_info["exon1_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon1_end"],"end":event_info["exon_aft_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_pre_end"],"end":event_info["exon2_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon2_end"],"end":event_info["exon_aft_start"]},ignore_index=True)
 elif (event_type=="A3" and event_info["strand"]=="+") or (event_type=="A5" and event_info["strand"]=="-"):
  if spladder_version==1:
   event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_const_end"],"end":event_info["exon_alt1_start"]},ignore_index=True)
   event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_const_end"],"end":event_info["exon_alt2_start"]},ignore_index=True)
  else:
   event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_const_end"],"end":event_info["exon_alt1_start"]},ignore_index=True)
   event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_const_end"],"end":event_info["exon_alt2_start"]},ignore_index=True)
 elif (event_type=="A3" and event_info["strand"]=="-") or (event_type=="A5" and event_info["strand"]=="+"):
  if spladder_version==1:
   event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_alt1_end"],"end":event_info["exon_const_start"]},ignore_index=True)
   event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_alt2_end"],"end":event_info["exon_const_start"]},ignore_index=True)
  else:
   event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_alt1_end"],"end":event_info["exon_const_start"]},ignore_index=True)
   event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_alt2_end"],"end":event_info["exon_const_start"]},ignore_index=True)
 return event_coords

def find_matching_transcripts(ensembl,gene_id,event_coords):
 try:
  transcript_ids=ensembl.transcript_ids_of_gene_id(gene_id)
 except:
  print('Warning: ' + gene_id + ' not found')
  transcript_ids=[]
 transcript_table=pd.DataFrame(columns=["coding","matching_isoform"],index=transcript_ids)
 for tid in transcript_ids:
  transcript=ensembl.transcript_by_id(tid)
  exons=pd.DataFrame(transcript.exon_intervals,columns=["start","end"]).sort_values(by="start")
  event_region=[event_coords.start.min(),event_coords.end.max()]
  exons=exons[find_overlap(exons,event_region)]
  junctions=pd.DataFrame(columns=["start","end"])
  junctions["end"]=exons.start[1:].values
  junctions["start"]=exons.end[0:-1].values
  if transcript.biotype!='protein_coding' or not transcript.contains_stop_codon or not transcript.contains_start_codon:
   transcript_table.loc[tid,"coding"]=False
  else:
   transcript_table.loc[tid,"coding"]=True
  #print(junctions)
  if junctions.equals(event_coords.loc[event_coords.isoform=="iso1",["start","end"]].reset_index(drop=True).astype(int)):
   transcript_table.loc[tid,"matching_isoform"]="iso1"
  elif junctions.equals(event_coords.loc[event_coords.isoform=="iso2",["start","end"]].reset_index(drop=True).astype(int)):
   transcript_table.loc[tid,"matching_isoform"]="iso2"
  else:
   transcript_table.loc[tid,"matching_isoform"]="none"
 return transcript_table
