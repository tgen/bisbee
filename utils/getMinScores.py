import pandas as pd
import numpy as np
import sys, os


scoresFile1=sys.argv[1]
scoresFile2=sys.argv[2]
outname=sys.argv[3]
stats=pd.DataFrame()

info_col=["contig","gene","strand","event_id","confirmed","event_jid"]
colnames=pd.read_csv(scoresFile1,nrows=1).columns
sample_names=colnames[range(np.where(colnames=='event_jid')[0][0]+1,len(colnames))]
info=pd.read_csv(scoresFile1,index_col="event_jid",usecols=info_col)
scores1=pd.read_csv(scoresFile1,index_col="event_jid",usecols=np.append("event_jid",sample_names))
scores2=pd.read_csv(scoresFile2,index_col="event_jid",usecols=np.append("event_jid",sample_names))
minScores=pd.DataFrame(index=scores1.index,columns=scores1.columns)
use1=scores1.abs()<=scores2.abs()
use2=scores1.abs()>scores2.abs()
minScores[use1]=scores1[use1]
minScores[use2]=scores2[use2]
minScores=info.join(minScores)
minScores.insert(column="event_jid",loc=5,value=minScores.index)
minScores.to_csv(outname + ".minScores.bisbeeOutlier.csv",index=False)
