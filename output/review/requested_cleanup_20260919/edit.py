from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
import hashlib,json
src=Path('/Users/lzz/Desktop/Repli-ATAC-seq-20260705.docx')
base=Path('/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/review/requested_cleanup_20260919');base.mkdir(exist_ok=True)
data=src.read_bytes();(base/'before.docx').write_bytes(data)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}; w='{'+ns['w']+'}'
with ZipFile(src) as z:
 parts={x.filename:z.read(x) for x in z.infolist()};infos=z.infolist()
d=E.fromstring(parts['word/document.xml']);old=E.fromstring(parts['word/document.xml']);c=E.fromstring(parts['word/comments.xml'])
remove={'16','26','10','20','13','14'}; ps=d.xpath('//w:body/w:p',namespaces=ns)
selected=[p for p in ps if set(p.xpath('.//w:commentRangeStart/@w:id',namespaces=ns))&remove]
assert len(selected)==6
p=ps[239]; before=''.join(p.xpath('.//w:t/text()',namespaces=ns))
start=before.index(' For Figure 8D and 8E,')
# Trim only the GO detail after the preserved EndNote citation.
pos=0
for t in p.xpath('.//w:t',namespaces=ns):
 text=t.text or ''; end=pos+len(text)
 if end>start:t.text=text[:max(0,start-pos)]
 pos=end
for t in p.xpath('.//w:t',namespaces=ns):
 if t.text and ' platform with default parameters ' in t.text:t.text=t.text.replace(' platform with default parameters ',' platform ')
expected=before[:start].replace(' platform with default parameters ',' platform ')
assert ''.join(p.xpath('.//w:t/text()',namespaces=ns))==expected
for p in selected:
 for col in p.xpath('.//w:color',namespaces=ns):
  if col.get(w+'val')=='FF0000':col.set(w+'val','000000')
for tag in ['commentRangeStart','commentRangeEnd','commentReference']:
 for el in d.findall('.//'+w+tag):
  if el.get(w+'id') in remove:
   parent=el.getparent();parent.remove(el)
   if tag=='commentReference' and parent.tag==w+'r' and all(x.tag==w+'rPr' for x in parent):parent.getparent().remove(parent)
for el in list(c):
 if el.get(w+'id') in remove:c.remove(el)
parts['word/document.xml']=E.tostring(d,xml_declaration=True,encoding='UTF-8',standalone=True)
parts['word/comments.xml']=E.tostring(c,xml_declaration=True,encoding='UTF-8',standalone=True)
out=base/'Repli-ATAC-seq-20260705.docx'
with ZipFile(out,'w') as z:
 for info in infos:z.writestr(info,parts[info.filename])
# Verify exact edit scope, preserved fields, images, and remaining comments.
ops=old.xpath('//w:body/w:p',namespaces=ns)
changed=[i for i,(a,b) in enumerate(zip(ops,ps),1) if E.tostring(a)!=E.tostring(b)]
assert changed==[240,252,273,291,298,299],changed
assert old.xpath('.//w:instrText/text()',namespaces=ns)==d.xpath('.//w:instrText/text()',namespaces=ns)
assert old.xpath('.//w:fldData/text()',namespaces=ns)==d.xpath('.//w:fldData/text()',namespaces=ns)
for i in changed:
 if i!=240: assert ops[i-1].xpath('.//w:t/text()',namespaces=ns)==ps[i-1].xpath('.//w:t/text()',namespaces=ns)
ids={e.get(w+'id') for e in c}; assert len(ids)==9 and not ids&remove
for tag in ['commentRangeStart','commentRangeEnd','commentReference']:
 assert {e.get(w+'id') for e in d.findall('.//'+w+tag)}==ids
for p in selected:assert not p.xpath('.//w:color[@w:val="FF0000"]',namespaces=ns)
with ZipFile(src) as a,ZipFile(out) as b:
 assert b.testzip() is None
 assert all(a.read(n)==b.read(n) for n in a.namelist() if n not in ['word/document.xml','word/comments.xml'])
report={'source_sha256':hashlib.sha256(data).hexdigest(),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'changed_paragraphs':changed,'removed_comment_ids':sorted(remove,key=int),'remaining_comments':len(ids),'simplified_paragraph':expected,'verification':'pass'}
(base/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
print(json.dumps(report,ensure_ascii=False,indent=2))
