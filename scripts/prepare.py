"""Align a frozen cohort, compute one SV PCA and prepare method-specific matrices."""
import argparse
import gzip
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd

PCS = [3, 5, 10]


def load_config(path):
    config = json.loads(path.read_text())
    for key in ('genotype', 'samples', 'groups'):
        config[key] = path.parent / config[key]
    if config['kind'] not in ('observed', 'synthetic'):
        raise ValueError('kind must be observed or synthetic.')
    if config['genotype_format'] not in ('dosage_tsv', 'vcf_bcf'):
        raise ValueError('genotype_format must be dosage_tsv or vcf_bcf.')
    variables = config['environments']
    if not variables or len(set(variables)) != len(variables) or any(
        re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', v) is None or
        v.lower() in {'group','sample','lon','lat','block','y','true','false','na','inf','nan',
                      'if','else','repeat','while','function','for','in','next','break','null'} or
        re.fullmatch(r'PC[0-9]+',v,re.I) for v in variables
    ):
        raise ValueError('Use distinct simple non-reserved environment column names.')
    for name, minimum in [('lfmm_k',1),('permutations',19),('min_called_per_group',1),('seed',1)]:
        if type(config[name]) is not int or config[name]<minimum:
            raise ValueError(f'{name} must be an integer >= {minimum}.')
    for key in ('adjust_geography','rda_blocked'):
        if not isinstance(config[key],bool):
            raise ValueError(f'{key} must be a JSON boolean.')
    for key in ('alpha','loading_z'):
        if not np.isfinite(config[key]) or config[key]<=0:
            raise ValueError(f'{key} must be positive and finite.')
    if config['alpha']>=1:
        raise ValueError('alpha must be less than one.')
    return config


def read_ids(path, columns):
    tab = pd.read_csv(path, sep='\t', dtype={'sample':str,'group':str,'block':str})
    if not set(columns)<=set(tab.columns) or tab.empty:
        raise ValueError(f'Missing required columns or rows in {path.name}.')
    tab = tab[columns]
    if tab.isna().any().any():
        raise ValueError(f'Missing required values in {path.name}.')
    for name in set(columns)&{'sample','group','block'}:
        if not tab[name].map(lambda x: bool(x.strip()) and x==x.strip() and '\n' not in x).all():
            raise ValueError(f'Invalid ID in {path.name}: {name}.')
    return tab


def dosage_table(path, ids):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as stream:
        header = stream.readline().rstrip('\r\n').split('\t')
    if header[0]!='sample' or len(set(header))!=len(header) or any(not x.strip() for x in header):
        raise ValueError('Dosage TSV needs a sample column and unique, nonempty marker IDs.')
    frame = pd.read_csv(path,sep='\t',dtype={'sample':str})
    if frame['sample'].duplicated().any() or not set(ids)<=set(frame['sample']):
        raise ValueError('Dosage table has duplicated IDs or omits a requested sample.')
    dose = frame.set_index('sample').loc[ids].to_numpy(dtype=np.float32)
    observed = dose[~np.isnan(dose)]
    if not np.isin(observed,[0,1,2]).all():
        raise ValueError('Dosage TSV accepts only observed diploid integers 0/1/2 or NA.')
    return dose, None, pd.DataFrame({'locus':header[1:]})


def query_bcf(path, ids, executable):
    samples = subprocess.check_output([executable,'query','-l',str(path)],text=True).splitlines()
    if len(set(samples))!=len(samples) or not set(ids)<=set(samples):
        raise ValueError('BCF sample IDs are duplicated or a requested sample is absent.')
    idx = [samples.index(s) for s in ids]
    dosage, ploidy, catalog = [], [], []
    query = [executable,'query','-f',r'%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT[\t%GT]\n',str(path)]
    with subprocess.Popen(query,stdout=subprocess.PIPE,text=True) as proc:
        try:
            for number,line in enumerate(proc.stdout,1):
                chrom,pos,end,record,ref,alt,*gt = line.rstrip('\n').split('\t')
                if len(gt)!=len(samples):
                    raise ValueError(f'GT column count differs at record {number}.')
                if alt == '.':
                    continue
                alts = alt.split(',')
                calls = np.full((len(alts),len(idx)), np.nan, dtype=np.float32)
                copy_number = np.zeros(len(idx), dtype=np.uint8)
                for row,i in enumerate(idx):
                    parts = re.split(r'[/|]',gt[i])
                    if '.' in parts:
                        continue
                    alleles = list(map(int,parts))
                    if not 1<=len(alleles)<=2 or not all(0<=a<=len(alts) for a in alleles):
                        raise ValueError(f'Unsupported called GT at record {number}.')
                    copy_number[row]=len(alleles)
                    calls[:,row]=np.bincount(alleles,minlength=len(alts)+1)[1:]*2/len(alleles)
                for j,alt_seq in enumerate(alts,1):
                    catalog.append((f'SV{number:09d}_A{j}',chrom,pos,end,record,ref,alt_seq,j))
                    dosage.append(calls[j-1])
                    ploidy.append(copy_number)
        except BaseException:
            proc.kill()
            raise
        if proc.wait()!=0:
            raise RuntimeError('bcftools query failed.')
    if not dosage:
        raise ValueError('No ALT allele columns found.')
    markers=pd.DataFrame(catalog,columns=['locus','chrom','position','end','record_id','ref','alt','alt_index'])
    return np.stack(dosage,axis=1),np.stack(ploidy,axis=1),markers


def pca(dose):
    n,m=dose.shape
    gram=np.zeros((n,n))
    valid_marker=np.zeros(m,dtype=bool)
    for start in range(0,m,2048):
        block=dose[:,start:start+2048].astype(float)
        counts=np.isfinite(block).sum(axis=0)
        means=np.divide(np.nansum(block,axis=0),counts,out=np.zeros(block.shape[1]),where=counts>0)
        block-=means
        block=np.nan_to_num(block,nan=0.0)
        variance=np.divide((block*block).sum(axis=0),counts,out=np.zeros(len(counts)),where=counts>0)
        keep=(counts>=2)&(variance>1e-12)
        valid_marker[start:start+len(keep)]=keep
        z=block[:,keep]/np.sqrt(variance[keep])
        gram+=z@z.T
    if valid_marker.sum()<10:
        raise ValueError('Fewer than 10 variable markers are available for PCA.')
    values,vectors=np.linalg.eigh(gram/valid_marker.sum())
    values=np.maximum(values[::-1],0)
    if len(values)<11 or values[9]<=values[0]*1e-10:
        raise ValueError('The genotype matrix supplies fewer than 10 estimable PCs.')
    pcs=vectors[:,::-1][:,:10]*np.sqrt(values[:10])
    for j in range(10):
        if pcs[np.abs(pcs[:,j]).argmax(),j]<0:
            pcs[:,j]*=-1
    return pcs,values,valid_marker


def write_matrix(values, rows, columns, path):
    pd.DataFrame(values,index=rows,columns=columns).to_csv(path,sep='\t',index=True,index_label='id')


def make_group_frequencies(dose, ploidy, sample_groups, groups, min_called):
    freqs=np.full((len(groups),dose.shape[1]),np.nan)
    for row,label in enumerate(groups):
        mask=sample_groups==label
        d=dose[mask]
        valid=np.isfinite(d)
        called=valid.sum(axis=0)
        n=np.where(valid,2 if ploidy is None else ploidy[mask],0)
        denominator=n.sum(axis=0)
        numerator=np.nansum(d*n/2,axis=0)
        freqs[row]=np.divide(numerator,denominator,out=np.full(len(called),np.nan),
                             where=(called>=min_called)&(denominator>0))
    return freqs


def prepare(config, out):
    variables=config['environments']
    geo=['lon','lat'] if config['adjust_geography'] else []
    blocks=['block'] if config['rda_blocked'] else []
    samples=read_ids(config['samples'],['sample','group'])
    groups=read_ids(config['groups'],['group',*variables,*geo,*blocks]).set_index('group')
    if samples['sample'].duplicated().any() or not groups.index.is_unique:
        raise ValueError('Frozen samples and group metadata require unique keys.')
    if set(samples['group'])!=set(groups.index):
        raise ValueError('Sample groups and group table must contain exactly the same group IDs.')
    groups=groups.sort_index()
    numeric=groups[variables+geo].to_numpy(float)
    if not np.isfinite(numeric).all():
        raise ValueError('Environment and geographic columns must be finite numeric values.')
    if geo and (not groups.lon.between(-180,180).all() or not groups.lat.between(-90,90).all()):
        raise ValueError('Longitude or latitude is outside its valid range.')
    counts=samples.groupby('group').size().reindex(groups.index)
    if (counts<config['min_called_per_group']).any():
        raise ValueError('A group has fewer individuals than min_called_per_group.')
    if len(groups)-1-10-len(variables)-len(geo)<1:
        raise ValueError('PC10 partial RDA requires positive residual degrees of freedom.')
    if len(samples)-1-10-len(variables)-len(geo)-config['lfmm_k']<1:
        raise ValueError('PC10 LFMM2 requires positive residual degrees of freedom.')
    ids=samples['sample'].tolist()
    if config['genotype_format']=='vcf_bcf':
        dose,ploidy,catalog=query_bcf(config['genotype'],ids,config.get('bcftools','bcftools'))
    else:
        dose,ploidy,catalog=dosage_table(config['genotype'],ids)
    pcs,variance,variable=pca(dose)
    pc_names=[f'PC{i}' for i in range(1,11)]
    pc_frame=pd.DataFrame(pcs,index=ids,columns=pc_names)
    pcs_group=pc_frame.groupby(samples.set_index('sample')['group']).mean().reindex(groups.index)
    f=make_group_frequencies(dose,ploidy,samples['group'].to_numpy(),groups.index,config['min_called_per_group'])
    keep_rda=np.isfinite(f).all(axis=0)&(np.var(np.nan_to_num(f),axis=0)>1e-12)
    keep_lfmm=np.isfinite(dose).all(axis=0)&variable
    if keep_lfmm.sum()<max(3,config['lfmm_k']) or keep_rda.sum()<2:
        raise ValueError('Too few eligible loci for LFMM2 or RDA.')
    out.mkdir(parents=True,exist_ok=True)
    catalog.assign(pca_variable=variable,lfmm_complete_variable=keep_lfmm,rda_complete_group_variable=keep_rda
                   ).to_csv(out/'marker_catalog.tsv.gz',sep='\t',index=False)
    write_matrix(dose[:,keep_lfmm],ids,catalog.locus[keep_lfmm],out/'lfmm_dosage.tsv.gz')
    write_matrix(f[:,keep_rda],groups.index,catalog.locus[keep_rda],out/'rda_frequency.tsv.gz')
    group_meta=groups.join(pcs_group).assign(n_individuals=counts)
    group_meta.to_csv(out/'group_design.tsv',sep='\t',index_label='id')
    subject_meta=samples.set_index('sample').join(groups,on='group').join(pc_frame)
    subject_meta.to_csv(out/'individual_design.tsv',sep='\t',index_label='id')
    pd.DataFrame({'component':np.arange(1,len(variance)+1),'eigenvalue':variance,
                  'fraction_total':variance/variance.sum()}).to_csv(out/'pca_eigenvalues.tsv',sep='\t',index=False)
    settings={k:v for k,v in config.items() if k not in ['genotype','samples','groups']}
    settings.update({'individuals':len(ids),'coordinate_groups':len(groups),'allele_columns':len(catalog),
                     'pca_markers':int(variable.sum()),'lfmm_markers':int(keep_lfmm.sum()),
                     'rda_markers':int(keep_rda.sum()),'pc_counts':PCS})
    (out/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    print(f'{len(ids)} individuals, {len(groups)} groups, LFMM2 {keep_lfmm.sum()} loci, RDA {keep_rda.sum()} loci')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('config',type=Path)
    ap.add_argument('output',type=Path)
    args=ap.parse_args()
    prepare(load_config(args.config.resolve()),args.output)
