import gzip, bisect
from collections import Counter
def op(p): return gzip.open(p,'rt') if str(p).endswith(('.gz','.bgz')) else open(p)
def gt(g):
    parts=g.replace('|','/').split('/')
    return None if any(x in ('.','') for x in parts) else [int(x) for x in parts]
def info_int(info,k):
    for kv in info.split(';'):
        if kv.startswith(k+'='):
            try: return int(kv.split('=')[1])
            except: return None
    return None
def parse(path):
    s=[];loci={}
    for line in op(path):
        if line.startswith('##'):continue
        if line.startswith('#CHROM'): s=line.rstrip().split('\t')[9:]; continue
        f=line.rstrip('\n').split('\t')
        ref,alt=f[3],f[4]; alleles=[ref]+([] if alt=='.' else alt.split(','))
        rlen=len(ref); per=info_int(f[7],'PERIOD'); gi=f[8].split(':').index('GT')
        rel={n:(None if gt(c.split(':')[gi]) is None else tuple(sorted(len(alleles[i])-rlen for i in gt(c.split(':')[gi])))) for n,c in zip(s,f[9:])}
        loci[(f[0],int(f[1]))]={'period':per,'filter':('PASS' if f[6] in ('.','PASS','') else f[6]),'rel':rel,'ref':ref,'reflen':rlen}
    return s,loci
_,ours=parse('/Users/jose/devel/pop_var_caller/benchmarks/ssr_tomato1/results_rerun_20260708/ours/cohort/cohort.ssr.vcf')
_,hip=parse('/Users/jose/devel/pop_var_caller/benchmarks/ssr_tomato1/results_ssr15k/hipstr/cohort.str.vcf.gz')
# catalog POS set from BED (start+1)
bed=set(); bedspan={}
for line in open('/Users/jose/devel/pop_var_caller/benchmarks/ssr_tomato1/results_ssr15k/hipstr/ssr_tomato1.hipstr_regions.bed'):
    c,st,en,per,score,name=line.rstrip('\n').split('\t')
    bed.add((c,int(st)+1)); bedspan.setdefault(c,[]).append((int(st)+1,int(en)))
for c in bedspan: bedspan[c].sort()
def bed_overlaps(c,s,e):
    v=bedspan.get(c)
    if not v: return False
    starts=[x[0] for x in v]; j=bisect.bisect_right(starts,e)-1
    while j>=0 and v[j][0]>=s-50:
        if v[j][0]<e and v[j][1]>s: return True
        j-=1
    return False
# pass index for ours
idx={}
for (c,p),r in ours.items():
    if r['filter']!='PASS' or not any(g is not None for g in r['rel'].values()): continue
    idx.setdefault(c,[]).append((p,p+r['reflen']))
for c in idx: idx[c].sort()
def ours_overlaps(c,s,e):
    v=idx.get(c)
    if not v: return False
    starts=[x[0] for x in v]; j=bisect.bisect_right(starts,e)-1
    while j>=0 and v[j][0]>=s-500:
        if v[j][0]<e and v[j][1]>s: return True
        j-=1
    return False
def hip_var(r):
    p=r['period']
    if not p: return False
    return any(g is not None and any(d!=0 and d%p==0 for d in g) for g in r['rel'].values())
drops=[]
for k,h in hip.items():
    if not hip_var(h): continue
    c,pos=k; o=ours.get(k)
    if o is not None and o['filter']=='PASS' and any(g is not None for g in o['rel'].values()): continue
    if o is not None and o['filter']!='PASS': continue
    if o is not None: continue
    if ours_overlaps(c,pos,pos+h['reflen']): continue
    drops.append((k,h))
print('genuine drops:',len(drops))
exact=sum(1 for (c,p),_ in drops if (c,p) in bed)
ov=sum(1 for (c,p),h in drops if bed_overlaps(c,p,p+h['reflen']))
print('  drop POS exactly in catalog-BED :',exact, f'({exact/len(drops):.1%})')
print('  drop tract overlaps a catalog-BED locus :',ov, f'({ov/len(drops):.1%})')
print('  NOT in catalog at all (no overlap):',len(drops)-ov, f'({(len(drops)-ov)/len(drops):.1%})')
