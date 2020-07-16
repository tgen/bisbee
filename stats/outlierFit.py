import sys
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
from scipy.optimize import minimize
from scipy.stats import betabinom

sample_file_exists = False;

counts_file = sys.argv[1]
max_beta = int(sys.argv[2])
outname = str(sys.argv[3])+'.bisbeeFit.csv'

if (len(sys.argv)>4):
	sampleFile = sys.argv[4]
	sample_file_exists = True;
    
data = pd.read_csv(counts_file)
iso1_idx = [col for col in data.columns if 'iso1' in col]
iso2_idx = [col for col in data.columns if 'iso2' in col]
sample_names = pd.DataFrame([x.replace("_iso1","") for x in iso1_idx])

if (sample_file_exists):
    inc_samples = pd.read_csv(sampleFile, header=None)
    sample_idx = sample_names.reset_index().set_index(0).loc[inc_samples[0], 'index'].values
else:
	sample_idx = range(0,len(sample_names))

sample_count = len(sample_idx)
event_count = data.shape[0]
print('sampleCount: '+str(sample_count))
print('eventCount: '+str(event_count))


iso1=data[iso1_idx]
iso2=data[iso2_idx]
info_idx=range(0,data.columns.get_loc("event_jid")+1)
list_addition=['alpha','beta','depth_min','depth_q05','depth_q25','depth_median','depth_q75','depth_q95','depth_max','psi_min','psi_q05','psi_q25','psi_median','psi_q75','psi_q95','psi_max','log-likelihood']
cnames = data.loc[:,:"event_jid"].columns.tolist()
cnames.extend(list_addition)
iso2_adjusted=iso2
iso2_adjusted.columns=iso1.columns
isosum = iso1+iso2_adjusted
depthQ=pd.DataFrame(list(map(lambda x: np.percentile(isosum.iloc[x], [0, 5, 25, 50, 75, 95, 100]), range(0,isosum.shape[0]))))
psi=iso1/isosum
meanPSI=psi.mean(axis=1,skipna=True)
sumDepth = isosum.sum(axis=1, skipna=True)
hIdx = meanPSI>0.5
hIdx=hIdx[hIdx].index.tolist()
mCounts=iso1
mCounts.iloc[hIdx] = iso2.iloc[hIdx]


def nLL(y, x, mCounts, isosum, max_beta): 
	return -sum(betabinom.logpmf(mCounts.iloc[x], isosum.iloc[x], (1-1/max_beta)/(1+np.exp(y[0]))+1/max_beta,(max_beta-1)/(1+np.exp(y[1]))+1))
def mle_bb_mode0(x, mCounts, isosum, max_beta):    
	x0=np.array([0,0]) 
	try:
		mleRes = minimize(nLL, x0, args = (x, mCounts, isosum, max_beta), method = "nelder-mead")
	except:
		try:
			mleRes = minimize(nLL, x0, method="bfgs")
		except:
			mleRes = None
	return mleRes

a= [1]*len(meanPSI)
a = pd.DataFrame(a)
a = a.astype(float)
b = list(map(lambda x: np.minimum(x, max_beta),sumDepth+1))
b = pd.DataFrame(b)
b = b.astype(float)

ll = pd.DataFrame([None]*len(meanPSI))
fIdx = mCounts.sum(axis=1)>0
fIdx = pd.DataFrame(fIdx[fIdx].index.tolist())


if(len(fIdx)>0):
	param = [None] * len(fIdx)
	for i in range(len(fIdx)):
		param[i] = mle_bb_mode0(fIdx.iloc[i][0], mCounts, isosum, max_beta)

coefs = pd.DataFrame(None, index = range(len(fIdx)), columns = ['v1','v2'])
for i in range(len(fIdx)):
    if param[i] != None:
        coefs['v1'][i] = param[i]['x'][0]
        coefs['v2'][i] = param[i]['x'][1]
    else:
        coefs['v1'][i] = None
        coefs['v2'][i] = None
    
idx = coefs['v1'].map(lambda x: x is not None)
idx=idx[idx].index.tolist()

if (len(idx)>0):
    for i in idx:
        a[0][fIdx[0][i]] = (1-1/max_beta)/(1+np.exp(coefs["v1"][i]))+1/max_beta
        b[0][fIdx[0][i]] = (max_beta-1)/(1+np.exp(coefs["v2"][i]))+1

alpha = a.copy()
alpha.columns = ['alpha']
alpha['alpha'][hIdx] = b[0][hIdx]
beta = b.copy()
beta.columns = ['beta']
beta['beta'][hIdx] = a[0][hIdx]

if(len(fIdx)>0):
    nIdx = coefs['v1'].map(lambda x: x is None)
    nIdx = nIdx[nIdx].index.tolist()
    if len(nIdx)>0:
        for i in range(len(nIdx)):
            alpha['alpha'][fIdx[0][nIdx[i]]] = None
            beta['beta'][fIdx[0][nIdx[i]]] = None
        
psiQ=pd.DataFrame(list(map(lambda x: np.nanpercentile(iso1.iloc[x]/isosum.iloc[x], [0, 5, 25, 50, 75, 95, 100]), range(0,isosum.shape[0]))))
dfout=pd.concat([data.iloc[:,info_idx],alpha,beta,depthQ,psiQ,ll], axis=1)
dfout.to_csv(outname)

