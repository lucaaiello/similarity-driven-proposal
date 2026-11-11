"""Prepare source/input files for a fresh cluster workspace; never run R or submit jobs."""
from pathlib import Path
import argparse,hashlib,tarfile
root=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser()
parser.add_argument('--output',default='dist/real_data_analysis_bundle.tgz')
parser.add_argument('--manifest-only',action='store_true')
args=parser.parse_args()
manifest=root/'Real Data Analysis/Cluster/source_files.sha256'
files=[]
for p in (root/'Real Data Analysis').rglob('*'):
    if not p.is_file() or p==manifest:continue
    rel=p.relative_to(root)
    if any(x in {'Results','Logs','__pycache__'} for x in rel.parts):continue
    if p.name=='submitted_jobs.txt' or p.suffix in {'.pyc','.log','.dll','.so','.o'}:continue
    files.append(p)
files.append(root/'Simulations/Figures_and_Tables/plot_theme.R')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest.write_text(''.join(f'{sha(p)}  {p.relative_to(root).as_posix()}\n' for p in sorted(files)),encoding='utf-8')
if args.manifest_only:
    print('Live cluster source/input manifest prepared; no R, compilation, MCMC or submission.')
else:
    target=(root/args.output).resolve()
    assert target.is_relative_to(root/'dist') and not target.exists(),'Choose a new filename inside dist/'
    target.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(target,'w:gz',format=tarfile.GNU_FORMAT) as tar:
        for p in sorted(files+[manifest]):
            info=tar.gettarinfo(str(p),arcname=p.relative_to(root).as_posix())
            info.mode=0o755 if p.suffix in {'.sh','.sbatch'} else 0o644
            info.uid=info.gid=0;info.uname=info.gname=''
            with p.open('rb') as f:tar.addfile(info,f)
    with tarfile.open(target) as tar:
        for p in files+[manifest]:
            with tar.extractfile(p.relative_to(root).as_posix()) as f:
                assert hashlib.sha256(f.read()).hexdigest()==sha(p)
    target.with_suffix(target.suffix+'.sha256').write_text(f'{sha(target)}  {target.name}\n',encoding='ascii')
    print('Cluster source/input bundle verified:',target)
