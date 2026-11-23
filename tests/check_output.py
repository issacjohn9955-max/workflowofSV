"""Check the synthetic cohort, method scopes and restricted permutations."""
from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('run',type=Path)
root=parser.parse_args().run
state=json.loads((root/'progress.json').read_text())
assert state['status']=='complete' and len(state['completed'])==8
prepared=root/'prepared'
groups=pd.read_csv(prepared/'group_design.tsv',sep='\t')
cfg=json.loads((prepared/'settings.json').read_text())
assert (cfg['individuals'],cfg['coordinate_groups'],cfg['lfmm_markers'],cfg['rda_markers'])==(112,28,149,159)
last_shuffle=None
for pc in (3,5,10):
    lf=root/f'lfmm_pc{pc}'
    rd=root/f'rda_pc{pc}'
    table=pd.read_csv(lf/'associations.tsv.gz',sep='\t')
    assert len(table)==149*5 and table.pc_count.eq(pc).all() and table.latent_k.eq(5).all()
    assert table.status.eq('tested').all() and table.q.between(0,1).all()
    fit=pd.read_csv(rd/'global.tsv',sep='\t').iloc[0]
    assert fit.groups==28 and fit.loci==159 and fit.residual_df==28-1-pc-5-2
    shuffles=pd.read_csv(rd/'permutations.tsv.gz',sep='\t').to_numpy()-1
    assert np.array_equal(groups.block.to_numpy()[shuffles],np.broadcast_to(groups.block.to_numpy(),shuffles.shape))
    if last_shuffle is not None:
        assert np.array_equal(shuffles,last_shuffle)
    last_shuffle=shuffles
print('PASS: 8 stages, frozen scopes, PC degrees of freedom and identical blocked group permutations.')
