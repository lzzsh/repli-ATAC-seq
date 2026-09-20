from pathlib import Path
from zipfile import ZipFile
from lxml import etree as E
from difflib import SequenceMatcher
from copy import deepcopy
import hashlib,json,shutil,re
B=Path('/Users/lzz/Documents/GitHub/repli-ATAC-seq/output/review/resolve_comments_20260919');OLD=B.parent/'word_edit_20260919'
SRC=Path('/Users/lzz/Desktop/Repli-ATAC-seq-20260705.docx');OUT=B/'Repli-ATAC-seq-20260705_resolved.docx'
assert SRC.read_bytes()==(OLD/'Repli-ATAC-seq-20260705_red.docx').read_bytes(),'Original changed; do not overwrite'
shutil.copy2(SRC,B/'Repli-ATAC-seq-20260705_before_resolution.docx')
# Restore original formatting from the black version; only comment-color changes differ.
with ZipFile(OLD/'Repli-ATAC-seq-20260705_edited.docx') as z: parts={i.filename:z.read(i.filename) for i in z.infolist()};infos=z.infolist()
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main';NS={'w':W};q=lambda x:'{'+W+'}'+x
root=E.fromstring(parts['word/document.xml']);ps={i:p for i,p in enumerate(root.findall('.//w:body/w:p',NS),1)}
def txt(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
original={i:txt(p) for i,p in ps.items()};wanted=dict(original)
def sub(i,a,b):
 assert a in wanted[i],(i,a)
 wanted[i]=wanted[i].replace(a,b)
def append(i,t):wanted[i]+=' '+t
sub(219,'Raw read counts were normalized for bin length and sequencing depth using a TPM-style calculation. The normalized signals from the ES, MS, and LS fractions were then normalized to the corresponding G1-phase signals within each sample. Biological replicates from the same condition were averaged to generate condition-level replication signals for downstream analyses.',
'For bin i in library s, the normalized signal was N(i,s) = 10⁶ × [C(i,s)/L(i)] / Σj[C(j,s)/L(j)], where C denotes read count and L denotes the interval length in kilobases, calculated as (end − start + 1)/1,000 in the analysis scripts. For equal-length bins, this reduces to counts per million. Normalized libraries from the same genotype and fraction were averaged before the ES, MS, and LS signals were divided by the corresponding G1 signal. The count matrix contained two G1 and two ES libraries and one MS and one LS library for WT, oscpp8-1, oscpp11-5, and oscpp11-8; oscpp8-3 had one library per fraction.')
sub(233,'DARs were identified by comparing the OsCPP8 and OsCPP11 mutant groups with the WT samples using edgeR v3.26.8 in R with a false discovery rate (FDR) < 0.05 and |log2(fold change)| > 0.3.',
'DARs were identified by comparing the OsCPP8 and OsCPP11 mutant groups with the WT samples using edgeR v3.26.8 in R. The DAR counts and overlaps in Figure 7A and 7B used a false discovery rate (FDR) < 0.05 without an additional fold-change cutoff. A further |log2(fold change)| > 0.3 cutoff was applied to the promoter-associated DARs used in Figure 8G and 8H.')
sub(235,'Read counts were normalized for peak length and sequencing depth using a TPM-style calculation, followed by log2 transformation and quantile normalization. TF-associated OCRs in each S-phase fraction were defined by their overlap with the OsCPP8 or OsCPP11 CUT&Tag peak set. Normalized signals were averaged across replicates for each genotype at the selected replication stage.',
'Read counts at a fixed set of 89,361 OCRs were normalized using the length- and library-size formula described for Repli-seq, with OCR length replacing bin length. Signals were transformed as log2(N + 1) and quantile-normalized with limma. The fixed OsCPP8 and OsCPP11 CUT&Tag peak sets were intersected with OCRs, and the resulting associated OCRs were selected by their RT annotations for stage-specific comparisons. Signals from two ES-labeled libraries were averaged; the MS and LS comparisons each used one library per genotype.')
sub(238,'TPM normalization was performed in R based on gene-level counts from featureCounts.',
'Gene-level raw counts from featureCounts were used as the edgeR input. Low-expression genes were removed using filterByExpr, and library composition was normalized with calcNormFactors using the trimmed mean of M-values method. Genotype effects were tested using a negative-binomial generalized linear model and a likelihood-ratio test (glmFit and glmLRT), with two RNA-seq libraries per genotype. TPM values were calculated separately for expression summaries. For heatmaps, log2 counts per million (prior count = 1) were standardized by gene to obtain z-scores.')
sub(240,'using a significance cutoff of FDR < 0.05 and |log2(fold change)| > 0.5.',
'using FDR < 0.05 for the significant-gene tables and an additional |log2(fold change)| > 1 for the stringent DEG exports. The expression heatmap used FDR < 0.05 together with the explicitly selected genes shown as annotations.')
append(240,'For Figure 8D and 8E, the GO output tables used ontology-specific backgrounds of 22,598 genes for biological process, 23,764 for molecular function, and 21,660 for cellular component. Bar lengths represent −log10 of the raw P values; the exported FDR values were used to retain terms at FDR < 0.05.')
append(246,'The regression in Figure 8F used 880 promoter-associated DAR–gene records corresponding to 814 unique genes, intersected with the FDR-significant RNA table; multiple DARs assigned to the same gene were retained as separate records.')
sub(249,'DARs associated with OsCPP8 binding and located within promoter regions were first filtered by fold change threshold (|log2FC_ATAC| > 0.3), and differentially expressed genes from RNA-seq were further filtered using the stricter threshold |log2FC_RNA| > 1 for visualization. Overlapping genes between the two datasets were identified and integrated with gene annotation information.',
'The FDR-significant (FDR < 0.05) promoter-proximal, OsCPP8-associated DAR set was further filtered using |log2FC_ATAC| > 0.3, and the FDR-significant RNA table was filtered using |log2FC_RNA| > 1. Unique gene identifiers were used for the overlap in Figure 8G, yielding 1,940 DAR-associated genes, 2,250 DEGs, and 205 shared genes. The shared genes were integrated with gene annotations, and genes with usable symbols were retained for the candidate-gene map in Figure 8H.')
sub(272,'Minor-allele-frequency distributions','Minor allele frequency (MAF) distributions')
sub(272,'Variants are grouped into four frequency bins: <1%, 1–5%, 5–20%, and 20–50%.',
'Variants are grouped into four non-overlapping frequency intervals: 0 ≤ MAF < 1%, 1% ≤ MAF < 5%, 5% ≤ MAF < 20%, and 20% ≤ MAF ≤ 50%.')
sub(273,'Common-to-rare variant ratios across RT categories.',
'Common-to-rare variant ratios across RT categories, with common variants defined as MAF ≥ 5% and rare variants as MAF < 5%. The analysis used the RiceVarMap rice4k All frequency dataset containing 17,397,026 SNP and indel records from 4,726 rice accessions. Variants were assigned by genomic position to RT-annotated OCR intervals, with overlapping intervals merged within each RT category to avoid duplicate counting within that category.')
sub(291,'Each point represents one biological replicate, colored by replication stage (ES/G1, MS, LS).',
'Each point represents one sequenced library, colored according to the sample-sheet labels ES/G1, MS, or LS. Each genotype contributes two ES/G1-labeled libraries, one MS library, and one LS library.')
wanted[291]='(A) PCA of ATAC-seq signals in WT, oscpp8-1, oscpp8-3, oscpp11-5, and oscpp11-8. Each point represents a library, colored by the sample-sheet labels ES/G1, MS, and LS (two, one, and one libraries per genotype, respectively).'
wanted[293]='(C–E) ATAC-seq signals at all 89,361 OCRs in WT and oscpp8-3 for ES (C), MS (D), and LS (E). Percentages above each panel show the proportion of OCRs with higher signal in oscpp8-3 (Up) or WT (Down).'
wanted[294]=''
wanted[295]=''
append(155,'Points represent the 880 DAR–gene records (814 unique genes) included in the regression.')
append(156,'Gene identifiers were deduplicated within each set; the intersection contains 205 genes.')
# Keep all edits outside the user's protected ES/MS/LS interpretation.
assert all(wanted[i]==original[i] for i in [132,134,166])

def write_text(p,new):
 old=txt(p)
 for tag,a,b,c,d in reversed(SequenceMatcher(None,old,new,autojunk=False).get_opcodes()):
  if tag=='equal':continue
  spans=[];pos=0
  for t in p.findall('.//w:t',NS):
   end=pos+len(t.text or '');spans.append((t,pos,end));pos=end
  if a==b:
   candidates=[(t,s,e) for t,s,e in spans if s<=a<=e and e>s];assert candidates
   t,s,e=candidates[-1] if a==0 else candidates[0];off=a-s;t.text=(t.text or '')[:off]+new[c:d]+(t.text or '')[off:];t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
  else:
   affected=[(t,s,e) for t,s,e in spans if s<b and e>a];assert affected
   for j,(t,s,e) in enumerate(affected):
    val=t.text or '';t.text=val[:max(a-s,0)]+(new[c:d] if j==0 else '')+val[min(b-s,len(val)):];t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
 assert txt(p)==new
for i,v in wanted.items():
 if v!=original[i]:write_text(ps[i],v)
# Native editable Word labels correct the document figure labels while leaving
# every embedded scientific image byte unchanged.
V='urn:schemas-microsoft-com:vml'
WP='http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'
def figure_label(pnum,key,x,y,width,height,text,size=6,bold=False,center=False):
 p=ps[pnum];extent=p.find('.//{'+WP+'}extent');fw=int(extent.get('cx'))/12700;fh=int(extent.get('cy'))/12700
 # Figure 8 is centered; Figure S5 is left aligned. Column width is 453.5 pt.
 left=((453.5-fw)/2 if pnum==149 else 0)+x*fw
 run=E.SubElement(p,q('r'));pict=E.SubElement(run,q('pict'))
 rect=E.SubElement(pict,'{'+V+'}rect',id='corrected_'+key,fillcolor='#FFFFFF',stroked='f')
 rect.set('style',f'position:absolute;margin-left:{left:.3f}pt;margin-top:{y*fh:.3f}pt;width:{width*fw:.3f}pt;height:{height*fh:.3f}pt;z-index:251659264;mso-position-horizontal-relative:text;mso-position-vertical-relative:text;mso-wrap-style:none')
 box=E.SubElement(rect,'{'+V+'}textbox',inset='0,0,0,0')
 content=E.SubElement(box,q('txbxContent'));pp=E.SubElement(content,q('p'));pr=E.SubElement(pp,q('pPr'))
 E.SubElement(pr,q('spacing'),{q('before'):'0',q('after'):'0',q('line'):str(round(size*22)),q('lineRule'):'exact'})
 E.SubElement(pr,q('jc'),{q('val'):'center' if center else 'left'})
 r=E.SubElement(pp,q('r'));rp=E.SubElement(r,q('rPr'));E.SubElement(rp,q('rFonts'),{q('ascii'):'Arial',q('hAnsi'):'Arial'})
 E.SubElement(rp,q('sz'),{q('val'):str(round(size*2))});E.SubElement(rp,q('color'),{q('val'):'000000'})
 if bold:E.SubElement(rp,q('b'))
 for k,line in enumerate(text.split('\n')):
  if k:E.SubElement(r,q('br'))
  E.SubElement(r,q('t')).text=line
figure_label(149,'fig8g_rna',0.430,0.919,0.090,0.017,'RNA DEGs',4.5,True,True)
figure_label(149,'fig8g_atac',0.565,0.916,0.120,0.025,'DAR-associated\ngenes',4.5,True,True)
for key,x,phase,up,down in [('C',0.034,'ES','46.3','53.7'),('D',0.372,'MS','46.5','53.5'),('E',0.707,'LS','46.3','53.7')]:
 figure_label(289,'figs5'+key+'_title',x,0.555,0.286,0.023,f'All OCRs: WT vs oscpp8-3 [{phase}]',6.5,True)
 figure_label(289,'figs5'+key+'_n',x,0.581,min(0.320,1.008-x),0.023,f'n = 89,361 OCRs | Up: {up}% | Down: {down}%',5.5)

# Resolved: complete computational details, verified figure H (no direct/indirect labels),
# AI declaration, and previously corrected bibliography. Dataset release remains pending.
resolved={8,11,12,15,17,18,19,21,22,23,24,25}
updates={
 7:'已按代码补充长度/文库归一化公式、重复平均顺序及各阶段文库数。仍需确认现稿 Figure 3G–J 最终 WRT 输入的低信号筛选规则：连续 WRT 与用于离散分类的阈值化 _norm 表不能混用；尚未找到与最终组合图完全对应的作图脚本。',
 10:'common/rare 阈值、rice4k All 群体、17,397,026 条 SNP/indel 记录、4,726 个品种及区间内去重规则均已补入。现有文件未注明可公开引用的 RiceVarMap release 或下载日期，仍需补充版本/获取信息。',
 16:'已依据原始 GO_TCX2-3_up/down.csv 补入 Figure 8D/E 的背景基因数，并确认图中显示原始 P 值、按导出 FDR < 0.05 筛选。仍需确认 Rice Gene Index 当次分析使用的多重校正方法及注释版本；Figure 5D/E 的原始 GO 输出尚未定位。',
 20:'已按 sampleSheet1.csv 补入每组两份 ES/G1、一份 MS、一份 LS 文库。ES/G1 是原样本表的 Factor 标签；其对应 G1、early-S 或混合门仍需实际分选记录确认，不能仅凭脚本命名判断。',
 17:'Figure S5C–E 的数据与 89,361 个 OCR 已核对。嵌入 TIFF 的标题仍为 All motif、说明仍为 motif only；应在最终排版源图中改为 All OCRs、n = 89,361 OCRs。已定位单面板 PDF，但尚未找到与当前整张组合图对应的可编辑文字对象。',
 18:'Figure 8G 的 2,250 个 RNA 基因、1,940 个 DAR 关联基因及 205 个交集已从输入表重算一致；正文和图注已补口径。嵌入 TIFF 仍保留 ATAC_DEG，应在最终排版源图中改为 DAR-associated genes；本次未改动科学图像像素。',
}
comments=E.fromstring(parts['word/comments.xml'])
for c in list(comments):
 cid=int(c.get(q('id')))
 if cid in resolved:comments.remove(c)
 elif cid in updates:
  pp=c.find('w:p',NS)
  for child in list(c):c.remove(child)
  pp=E.SubElement(c,q('p'));rr=E.SubElement(pp,q('r'));tt=E.SubElement(rr,q('t'));tt.text=updates[cid]
for tag in ['commentRangeStart','commentRangeEnd','commentReference']:
 for n in list(root.findall('.//w:'+tag,NS)):
  if int(n.get(q('id'))) in resolved:
   parent=n.getparent();parent.remove(n)
   if tag=='commentReference' and parent.tag==q('r') and all(x.tag==q('rPr') for x in parent):parent.getparent().remove(parent)
# Apply red only where unresolved comments remain; all completed paragraphs retain baseline font.
remaining={int(c.get(q('id'))) for c in comments};anchors={}
for i,p in ps.items():
 ids={int(n.get(q('id'))) for n in p.findall('w:commentRangeStart',NS)}
 if ids:
  anchors[i]=sorted(ids)
  for run in p.findall('.//w:r',NS):
   if run.find('w:t',NS) is None:continue
   pr=run.find('w:rPr',NS)
   if pr is None:pr=E.Element(q('rPr'));run.insert(0,pr)
   color=pr.find('w:color',NS)
   if color is None:color=E.SubElement(pr,q('color'))
   color.attrib.clear();color.set(q('val'),'FF0000')
# Italicize newly inserted mutant names without touching field instructions or drawings.
for i in [219,233,235,238,240,246,249,272,273,291,155,156]:
 p=ps[i];text=txt(p);intervals=[m.span() for m in re.finditer(r'\b(?:oscpp(?:8|11)(?:-\d+)?|wox11)\b',text)];pos=0
 for node in list(p.findall('.//w:t',NS)):
  val=node.text or '';a=pos;b=a+len(val);pos=b;overlap=[(max(a,s)-a,min(b,e)-a) for s,e in intervals if s<b and e>a]
  if not overlap:continue
  run=node.getparent()
  if run.tag!=q('r') or len(run.findall('w:t',NS))!=1 or any(c.tag not in [q('rPr'),q('t')] for c in run):continue
  cuts=sorted({0,len(val),*[v for pair in overlap for v in pair]});parent=run.getparent();idx=parent.index(run)
  for s,e in zip(cuts,cuts[1:]):
   if s==e:continue
   clone=deepcopy(run);clone.find('w:t',NS).text=val[s:e]
   if any(s>=x and e<=y for x,y in overlap):
    pr=clone.find('w:rPr',NS)
    if pr is None:pr=E.Element(q('rPr'));clone.insert(0,pr)
    it=pr.find('w:i',NS)
    if it is None:it=E.SubElement(pr,q('i'))
    it.set(q('val'),'1')
   parent.insert(idx,clone);idx+=1
  parent.remove(run)
parts['word/document.xml']=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
parts['word/comments.xml']=E.tostring(comments,xml_declaration=True,encoding='UTF-8',standalone=True)
with ZipFile(OUT,'w') as z:
 for info in infos:z.writestr(info,parts[info.filename])
audit={'source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(OUT.read_bytes()).hexdigest(),'resolved_comment_ids':sorted(resolved),'remaining_comment_ids':sorted(remaining),'remaining_anchors':anchors,'changes':[{'paragraph':i,'before':original[i],'after':v} for i,v in wanted.items() if original[i]!=v]}
(B/'resolution_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in audit.items() if k!='changes'},ensure_ascii=False,indent=2));print('paragraphs changed',len(audit['changes']))
