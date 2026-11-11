"""Build a self-contained, verified journal archive, without Git metadata/LFS pointers.

Python standard library only. Does not stage, commit, upload or run R.
"""
from pathlib import Path
import argparse,csv,hashlib,json,zipfile
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--output',default='dist/reproduction_package.zip')
args=parser.parse_args()
target=(ROOT/args.output).resolve()
if not target.is_relative_to(ROOT/'dist'):raise SystemExit('Output must be within this repository/dist.')
if target.exists():raise SystemExit('Choose a new archive filename; existing archive preserved.')
def selected():
    files=[]
    for name in ['New_simulations','Simulations','Real Data Analysis','JCGS_submission','reproduce','docs']:
        for p in (ROOT/name).rglob('*'):
            if not p.is_file():continue
            r=p.relative_to(ROOT)
            if any(x in {'.rcpp_cache','__pycache__'} or x.startswith('sourceCpp-') for x in r.parts):continue
            if p.suffix.lower() in {'.dll','.so','.o','.obj','.pyc','.log','.fls','.blg','.fdb_latexmk','.out'}:continue
            if name=='JCGS_submission' and p.stem in {'response_letter','cover_letter'}:continue
            files.append(p)
    files += [ROOT/n for n in ['README.md','LICENSING.md','CITATION.cff','.gitignore','.gitattributes','asa_jcgs.bst']]
    return sorted(set(files))
files=selected()
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
manifest={p.relative_to(ROOT).as_posix():{'bytes':p.stat().st_size,'sha256':digest(p)} for p in files}
for p in files:
    if p.suffix=='.rds':
        with p.open('rb') as f:prefix=f.read(100)
        if prefix.startswith(b'version https://git-lfs.github.com/spec/'):raise SystemExit(f'LFS pointer: run git lfs pull: {p}')
target.parent.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    z.writestr('ARCHIVE_MANIFEST_SHA256.json',json.dumps(manifest,indent=2))
with zipfile.ZipFile(target) as z:
    for name,entry in manifest.items():
        h=hashlib.sha256()
        with z.open(name) as f:
            for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
        if h.hexdigest()!=entry['sha256']:raise RuntimeError('Archive verification failed: '+name)
archive_hash=digest(target)
target.with_suffix(target.suffix+'.sha256').write_text(f'{archive_hash}  {target.name}\n',encoding='ascii')
print(f'VERIFIED JOURNAL ARCHIVE: {len(files)} files; {target.stat().st_size/1024**2:.1f} MiB\n{target}')
