from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
import json,re,hashlib
B=Path('/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/review/resolve_comments_20260919');N={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math','v':'urn:schemas-microsoft-com:vml'};q=lambda t:'{'+N['w']+'}'+t
oldpath=B/'Repli-ATAC-seq-20260705_before_resolution.docx';newpath=B/'Repli-ATAC-seq-20260705_resolved.docx'
def read(p):
 with ZipFile(p) as z:
  assert z.testzip() is None
  return {n:z.read(n) for n in z.namelist()}
a=read(oldpath);b=read(newpath);ar=E.fromstring(a['word/document.xml']);br=E.fromstring(b['word/document.xml']);aps=ar.findall('.//w:body/w:p',N);bps=br.findall('.//w:body/w:p',N)
assert len(aps)==len(bps)==343
text=lambda p:''.join(p.xpath('.//w:t/text()',namespaces=N))
audit=json.loads((B/'resolution_audit.json').read_text());expected={v['paragraph']:v['after'] for v in audit['changes']}
for i,(ap,bp) in enumerate(zip(aps,bps),1):
 assert E.tostring(ap.find('w:pPr',N))==E.tostring(bp.find('w:pPr',N)),i
 if i not in [149,289]:assert text(bp)==expected.get(i,text(ap)),i
for i in [132,134,166]:assert text(aps[i-1])==text(bps[i-1])
assert ar.xpath('.//m:t/text()',namespaces=N)==br.xpath('.//m:t/text()',namespaces=N)
assert ar.xpath('.//w:instrText/text()',namespaces=N)==br.xpath('.//w:instrText/text()',namespaces=N)
for x,y in zip(ar.findall('.//w:tbl',N),br.findall('.//w:tbl',N)):assert E.tostring(x)==E.tostring(y)
media=[n for n in a if n.startswith('word/media/')];assert all(a[n]==b[n] for n in media)
assert [n for n in a if a[n]!=b[n]]==['word/document.xml','word/comments.xml']
cr=E.fromstring(b['word/comments.xml']);ids={int(c.get(q('id'))) for c in cr}
assert ids==set(audit['remaining_comment_ids']) and len(ids)==15
for tag in ['commentRangeStart','commentRangeEnd','commentReference']:
 got=[int(n.get(q('id'))) for n in br.findall('.//w:'+tag,N)];assert len(got)==len(ids) and set(got)==ids
for i,p in enumerate(bps,1):
 if p.find('w:commentRangeStart',N) is not None:
  assert all(r.find('w:rPr/w:color',N).get(q('val'))=='FF0000' for r in p.findall('.//w:r',N) if r.find('w:t',N) is not None),i
 else:
  if i in [235,238,249,258,272,293,156,157,330,318,303,314]:
   assert not p.xpath('.//w:color[@w:val="FF0000"]',namespaces=N),i
full='\n'.join(text(p) for p in bps[:301]);abbrevs=['RT','TFs','OCRs','WT','WRT','GO','TEs','TPM','TSS','PCA','NIP','ZH11','EdU','AF488','DEGs','DARs','FDR','MAF'];counts={s:full.count('('+s+')') for s in abbrevs};assert all(c==1 for c in counts.values()),counts
labels=br.findall('.//v:rect',N);assert len(labels)==8
assert 'ATAC_DEG' not in full and 'All motif' not in full and 'DAR-associated\ngenes'.replace('\n','') in full
assert full.count('All OCRs:')==3 and full.count('n = 89,361 OCRs | Up:')==3
p=Path('/Users/lzz/Documents/GitHub/repli-ATAC-seq/agent.md');s=p.read_text();assert 'ssh -p 10013 liaozizhuo@172.16.78.132' in s and 'DATA_INVENTORY.md' in s
out={'zip_integrity':'pass','paragraphs':343,'removed_comments':12,'remaining_comments':15,'native_editable_figure_labels':8,'unchanged_scientific_media':len(media),'citation_fields':len(br.findall('.//w:instrText',N)),'protected_ES_MS_LS_text':'unchanged','equations_tables_paragraph_properties':'unchanged','abbreviations_defined_once':counts,'sha256':hashlib.sha256(newpath.read_bytes()).hexdigest()}
(B/'verification.json').write_text(json.dumps(out,ensure_ascii=False,indent=2));print(json.dumps(out,ensure_ascii=False,indent=2))
