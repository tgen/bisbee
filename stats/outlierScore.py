import sys, pandas
from scipy.stats import betabinom
import numpy as np
from itertools import compress
from operator import add

fit_file = sys.argv[1]
counts_file =sys.argv[2]
outname = sys.argv[3]+".bisbeeOutlier.csv"

psi_stats = pandas.read_csv(fit_file)

data = pandas.read_csv(counts_file)
iso1_idx = [col for col in data.columns if 'iso1' in col]
iso2_idx = [col for col in data.columns if 'iso2' in col]
sample_names= [x.replace("_iso1","") for x in iso1_idx]

sample_count=len(sample_names)
event_count = data.shape[0]

print("sampleCount: "+str(sample_count))
print("eventCount: "+str(event_count))

iso1=data[iso1_idx]
iso2=data[iso2_idx]
info_idx = [col for col in data.columns if 'iso1' not in col and 'iso2' not in col]

idx= psi_stats.reset_index().set_index('event_jid').loc[data["event_jid"], 'index'].values
alpha = psi_stats.alpha[idx] #
beta = psi_stats.beta[idx]
scores = pandas.DataFrame(index=range(0,len(alpha)), columns = sample_names)

idx1= (alpha>=beta)
idx1= list(compress(range(len(idx1)),idx1))
idx2= (alpha<beta)
idx2= list(compress(range(len(idx2)),idx2))

if(len(idx1) >0):
	scores.iloc[idx1]=list(map(lambda x: np.log(np.maximum(betabinom.cdf(iso1.iloc[x],list(map(add, iso1.iloc[x],iso2.iloc[x])),alpha[x],beta[x]),sys.float_info.min)),idx1))
if(len(idx2)>0):
	scores.iloc[idx2]=list(map(lambda x: -np.log(np.maximum(betabinom.cdf(iso2.iloc[x],list(map(add, iso1.iloc[x],iso2.iloc[x])),beta[x],alpha[x]),sys.float_info.min)),idx2))
iso1.columns =iso2.columns
scores[iso1+iso2 ==0]=float("nan")

dfout=pandas.concat([data[info_idx],scores], axis=1)
dfout.to_csv(outname)
