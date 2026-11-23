"""Summarize the three PC specifications on matched method-specific test scopes."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from prepare import PCS

root,out=map(Path,sys.argv[1:])
rda=[]
lfmm=[]
overlap=[]
sets={}
for pc in PCS:
    l=root/f'lfmm_pc{pc}'
    r=root/f'rda_pc{pc}'
    d=pd.read_csv(l/'diagnostics.tsv',sep='\t')
    pairs=pd.read_csv(l/'discoveries.tsv',sep='\t')
    lfmm.append(d)
    global_=pd.read_csv(r/'global.tsv',sep='\t')
    axes=pd.read_csv(r/'axes.tsv',sep='\t')
    tails=pd.read_csv(r/'exploratory_candidates.tsv',sep='\t')
    rda.append(global_.assign(eligible_axes=axes.eligible.sum(),exploratory_loci=len(tails)))
    sets[('LFMM2',pc)]=set(zip(pairs.locus,pairs.environment))
    sets[('RDA',pc)]=set(tails.locus)
    inter=set(pairs.locus)&set(tails.locus)
    overlap.append({'method':'LFMM2 / RDA','pc_a':pc,'pc_b':pc,'unit':'locus',
                    'intersection':len(inter),'union':len(set(pairs.locus)|set(tails.locus)),
                    'comparison':'descriptive cross-method support'})
for method in ['LFMM2','RDA']:
    for a,b in [(3,5),(3,10),(5,10)]:
        x,y=sets[method,a],sets[method,b]
        overlap.append({'method':method,'pc_a':a,'pc_b':b,'unit':'locus-environment' if method=='LFMM2' else 'locus',
                        'intersection':len(x&y),'union':len(x|y),'comparison':'PC sensitivity'})
tables=pd.DataFrame(overlap)
tables['jaccard']=tables.intersection/tables.union.replace(0,np.nan)
pd.concat(rda).to_csv(out/'rda_comparison.tsv',sep='\t',index=False)
pd.concat(lfmm).to_csv(out/'lfmm_comparison.tsv',sep='\t',index=False)
tables.to_csv(out/'overlap.tsv',sep='\t',index=False,na_rep='NA')
