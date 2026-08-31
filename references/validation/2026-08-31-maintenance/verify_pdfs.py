from pathlib import Path
from PIL import Image, ImageStat
import subprocess,sys
run=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
files=sorted(run.rglob('*.pdf'))
assert files
for i,p in enumerate(files):
 assert p.stat().st_size>1000,p
 subprocess.run(['pdfinfo',str(p)],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,check=True)
 target=out/f'{i:02d}'
 subprocess.run(['pdftoppm','-f','1','-singlefile','-r','54','-png',str(p),str(target)],stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,check=True)
 im=Image.open(target.with_suffix('.png')).convert('RGB')
 assert max(ImageStat.Stat(im).stddev)>1,p
 print(f'{i:02d}: {p.relative_to(run)} | nonblank first page')
print(f'PASS: {len(files)} PDFs have nonempty, renderable, nonblank first pages.')
