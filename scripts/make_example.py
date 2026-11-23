"""Generate a deterministic, fictional cohort with group-level environment."""
from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument('directory',type=Path)
args=ap.parse_args()
out=args.directory
out.mkdir(exist_ok=True,parents=True)
rng=np.random.default_rng(987654)
n_groups=28
per_group=4
n=n_groups*per_group
labels=[f'DEMO_GROUP_{i+1:02d}' for i in range(n_groups)]
sample_ids=[f'DEMO_SAMPLE_{i+1:03d}' for i in range(n)]
metadata=pd.DataFrame(rng.normal(size=(n_groups,5)),columns=['bio3','bio15','bio18','wind','srad'])
metadata.insert(0,'group',labels)
metadata['lon']=rng.uniform(60,120,n_groups)
metadata['lat']=rng.uniform(8,48,n_groups)
metadata['block']=[f'B{i//7+1}' for i in range(n_groups)]
samples=pd.DataFrame({'sample':sample_ids,'group':np.repeat(labels,per_group)})
probability=rng.uniform(.12,.55,(n,160))
probability[:,:8]=np.clip(probability[:,:8]+.17*np.repeat(metadata.bio3.to_numpy(),per_group)[:,None],.02,.98)
genotype=rng.binomial(2,probability).astype(float)
genotype[:3,140]=np.nan
genotype[:1,150:]=np.nan
pd.DataFrame(genotype,index=sample_ids,columns=[f'DEMO_SV_{i+1:04d}' for i in range(160)]).to_csv(
    out/'dosage.tsv',sep='\t',index_label='sample',na_rep='NA')
samples.to_csv(out/'samples.tsv',sep='\t',index=False)
metadata.to_csv(out/'groups.tsv',sep='\t',index=False)
config=dict(kind='synthetic',genotype_format='dosage_tsv',genotype='dosage.tsv',samples='samples.tsv',
            groups='groups.tsv',environments=['bio3','bio15','bio18','wind','srad'],lfmm_k=5,
            permutations=39,seed=20261123,min_called_per_group=3,adjust_geography=True,
            rda_blocked=True,alpha=.05,loading_z=3)
(out/'run.json').write_text(json.dumps(config,indent=2)+'\n')
