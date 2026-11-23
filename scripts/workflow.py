"""Run the SV sensitivity workflow with checked, reusable stage checkpoints."""
from datetime import datetime
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys

from prepare import PCS, load_config

HERE=Path(__file__).resolve().parent


def digest(path):
    with path.open('rb') as stream:
        return sha256_file(stream)


def sha256_file(stream):
    h=sha256()
    for block in iter(lambda:stream.read(2**20),b''):
        h.update(block)
    return h.hexdigest()


def json_write(path,data):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data,indent=2)+'\n')
    temporary.replace(path)


def products(directory):
    return {p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(directory.iterdir())
            if p.is_file() and p.name!='checkpoint.json'}


def log(message):
    print(f'[{datetime.now().astimezone().isoformat(timespec="seconds")}] {message}',flush=True)


def software(config):
    rscript=os.getenv('SV_RSCRIPT','Rscript')
    packages=subprocess.check_output([rscript,'-e',
        'p<-Sys.getenv("SV_R_LIBRARY"); if(nzchar(p)) .libPaths(c(p,.libPaths())); '
        'cat(R.version.string, "\\n"); for(x in c("LEA","vegan","permute","jsonlite")) '
        'cat(x,as.character(packageVersion(x)),"\\n")'],text=True)
    fingerprint={'R':packages,'python':sys.version,'rscript':shutil.which(rscript)}
    if config['genotype_format']=='vcf_bcf':
        fingerprint['bcftools']=subprocess.check_output([config.get('bcftools','bcftools'),'--version'],text=True)
    import numpy
    import pandas
    fingerprint.update(numpy=numpy.__version__,pandas=pandas.__version__)
    return rscript,fingerprint


def run(path,out):
    settings=load_config(path)
    out.mkdir(parents=True,exist_ok=True)
    with (out/'workflow.lock').open('a') as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'Another workflow uses {out}.')
        rscript,versions=software(settings)
        fingerprints={key:{'path':str(settings[key].resolve()),'sha256':digest(settings[key])}
                      for key in ['genotype','samples','groups']}
        provenance={'config_sha256':digest(path),'files':fingerprints,'software':versions,
                    'scripts':{p.name:digest(p) for p in sorted([*HERE.glob('*.py'),*HERE.glob('*.R')])}}
        receipt=out/'run_fingerprint.json'
        if receipt.exists() and json.loads(receipt.read_text())!=provenance:
            raise ValueError('Inputs, code, configuration or software changed. Choose a new output directory.')
        json_write(receipt,provenance)
        python=sys.executable
        stages=[('prepared',[python,str(HERE/'prepare.py'),str(path)])]
        for pc in PCS:
            for method in ['lfmm','rda']:
                stages.append((f'{method}_pc{pc}',[rscript,str(HERE/'models.R'),method,str(pc),str(out/'prepared')]))
        stages.append(('summary',[python,str(HERE/'summarize.py'),str(out)]))
        total=len(stages)
        finished=[]
        receipts={}
        for name,command in stages:
            target=out/name
            status={'status':'running','stage':name,'completed':finished,'total':total}
            json_write(out/'progress.json',status)
            marker=target/'checkpoint.json'
            try:
                existing=json.loads(marker.read_text()) if marker.is_file() else {}
            except (ValueError,OSError):
                existing={}
            parents=[] if name=='prepared' else (list(receipts)[1:] if name=='summary' else ['prepared'])
            inputs={parent:receipts[parent] for parent in parents}
            if existing and existing=={'files':products(target),'inputs':inputs}:
                log(f'REUSE {name}')
            else:
                log(f'START {name}')
                try:
                    with TemporaryDirectory(prefix=f'.{name}.',dir=out) as tmp:
                        build=Path(tmp)/'products'
                        build.mkdir()
                        subprocess.run([*command,str(build)],check=True)
                        manifest=products(build)
                        if not manifest or any(entry['bytes']==0 for entry in manifest.values()):
                            raise RuntimeError(f'Empty output detected in {name}.')
                        json_write(build/'checkpoint.json',{'files':manifest,'inputs':inputs})
                        if target.exists():
                            shutil.rmtree(target)
                        build.rename(target)
                except Exception:
                    json_write(out/'progress.json',{**status,'status':'failed'})
                    raise
                log(f'DONE {name}')
            receipts[name]=digest(marker)
            finished.append(name)
        json_write(out/'progress.json',{'status':'complete','stage':None,'completed':finished,'total':total})
        log(f'COMPLETE {total}/{total}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    run(args.config.resolve(),args.output.resolve())
