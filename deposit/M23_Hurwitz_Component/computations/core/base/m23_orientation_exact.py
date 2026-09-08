"""Exact orientation certificate for the M23 refinement.

The 3A witness is encoded by an explicit verified word in the standard
generators.  Every class, factorization, orbit, and transformation is then
reconstructed exactly.  Canonical sorting makes the transcript independent
of Python hash randomization and execution order.
"""
from collections import deque
import math
N=23
I=bytes(range(N+1))

def comp(p,q):
    return bytes([0]+[p[q[i]] for i in range(1,N+1)])
def inv(p):
    r=[0]*(N+1)
    for i in range(1,N+1): r[p[i]]=i
    return bytes(r)
def conj(x,g): # g^-1 x g
    return comp(comp(inv(g),x),g)
def from_cycles(cycles):
    p=list(range(N+1))
    for cyc in cycles:
        for i,x in enumerate(cyc): p[x]=cyc[(i+1)%len(cyc)]
    return bytes(p)
def ctype(p):
    seen=[False]*(N+1); lens=[]
    for i in range(1,N+1):
        if not seen[i]:
            j=i;l=0
            while not seen[j]: seen[j]=True;l+=1;j=p[j]
            if l>1:lens.append(l)
    return tuple(sorted(lens))
def cycstr(p):
    seen=[False]*(N+1); out=[]
    for i in range(1,N+1):
        if not seen[i]:
            cyc=[];j=i
            while not seen[j]: seen[j]=True;cyc.append(j);j=p[j]
            if len(cyc)>1: out.append('('+','.join(map(str,cyc))+')')
    return ''.join(out) or '()'
def order(p):
    o=1;seen=[False]*(N+1)
    for i in range(1,N+1):
        if not seen[i]:
            j=i;l=0
            while not seen[j]:seen[j]=True;l+=1;j=p[j]
            o=math.lcm(o,l)
    return o

a=from_cycles([(1,2),(3,4),(7,8),(9,10),(13,14),(15,16),(19,20),(21,22)])
b=from_cycles([(1,16,11,3),(2,9,21,12),(4,5,8,23),(6,22,14,18),(13,20),(15,17)])
bi=inv(b)
# Explicit 3A witness word in the alphabet (a,b,b^-1).
WITNESS_WORD=(1,0,0,0,2,0,0,1,0,0,0,2,1,0,0,2,0,1,2,1,0,0,0,1,2,1,1,2,0,2,2,0,0,0,0,2,2,0,0,0,1,2,1,1,2,2,1,1,1,0,0,1,1,0,0,0,1,2,1,0)
u0=I
for j in WITNESS_WORD:
    u0=comp(u0,(a,b,bi)[j])
assert ctype(u0)==(3,3,3,3,3,3)
print('u0 explicit word length',len(WITNESS_WORD),cycstr(u0))
# conjugacy class of 3A under a,b (and inverse b not necessary)
gens=[a,b]
class3={u0}
q=deque([u0])
while q:
    x=q.popleft()
    for s in gens:
        y=conj(x,s)
        if y not in class3:
            class3.add(y);q.append(y)
print('class3 size',len(class3))
# fixed h=a; factorizations u*v=h => v=u^-1 h
facts=[]; fact_index={}
for u in sorted(class3):
    v=comp(inv(u),a)
    if ctype(v)==(3,3,3,3,3,3):
        idx=len(facts);facts.append((u,v));fact_index[(u,v)]=idx
print('facts',len(facts))
assert len(facts)==2688
# conjugacy class of h and Schreier transporters x=t^-1 h t
trans={a:I}; q=deque([a])
# use a,b and inverses for better Schreier gens
actgens=[a,b,bi]
sgens=[]
while q:
    x=q.popleft(); t=trans[x]
    for s in actgens:
        y=conj(x,s)
        ts=comp(t,s)
        if y not in trans:
            trans[y]=ts;q.append(y)
        else:
            ty=trans[y]
            c=comp(ts,inv(ty))
            if c!=I:
                # assert centralizes h
                if conj(a,c)!=a:
                    raise RuntimeError('bad schreier')
                sgens.append(c)
print('class2 size',len(trans),'raw stabilizer gens',len(sgens))
# reduce stabilizer generators by unique and inverses; action on factorizations
uniq=[]; seen=set()
for c in sgens:
    if c not in seen:
        seen.add(c);seen.add(inv(c));uniq.append(c)
print('uniq stab gens',len(uniq))
# greedily keep gens that enlarge subgroup orbit of first factorization maybe enough to full C-orbits
# compute full orbits on fact set under all uniq generators; to optimize derive permutation maps only as needed
fact_lookup={u:i for i,(u,v) in enumerate(facts)}  # v determined by u
# conjugating u by c; v follows automatically
maps=[]
for k,c in enumerate(uniq):
    ci=inv(c)
    arr=[]; ok=True
    for u,v in facts:
        uu=comp(comp(ci,u),c)
        idx=fact_lookup.get(uu)
        if idx is None:
            ok=False;break
        arr.append(idx)
    if ok: maps.append(arr)
print('action maps',len(maps))
# orbit partition
orbit=[-1]*len(facts); orbits=[]
for i in range(len(facts)):
    if orbit[i]>=0: continue
    oid=len(orbits); orbit[i]=oid; qq=deque([i]); members=[]
    while qq:
        j=qq.popleft();members.append(j)
        for arr in maps:
            k=arr[j]
            if orbit[k]<0:orbit[k]=oid;qq.append(k)
    orbits.append(members)
print('orbit sizes',sorted(map(len,orbits)))
# test transformations

def oid_of_pair(u,v):
    j=fact_index.get((u,v)); return None if j is None else orbit[j]
counts={}
for name,fun in {
    'reverse_inverse':lambda u,v:(inv(v),inv(u)),
    'swap':lambda u,v:(v,u),
    'inverse_each':lambda u,v:(inv(u),inv(v)),
    'hurwitz_q':lambda u,v:(v, comp(comp(inv(v),u),v)), # (u,v)->(v,v^-1 u v), product same
    'hurwitz_q_inv':lambda u,v:(comp(comp(u,v),inv(u)),u),
}.items():
    transmat=[[0]*len(orbits) for _ in orbits]; invalid=0
    for i,(u,v) in enumerate(facts):
        uu,vv=fun(u,v)
        j=fact_index.get((uu,vv))
        if j is None: invalid+=1
        else: transmat[orbit[i]][orbit[j]]+=1
    print(name,'invalid',invalid,'matrix',transmat)
# Explicit representative each orbit and reverse_inverse destination
for oid,members in enumerate(orbits):
    u,v=facts[members[0]]
    print('ORBIT',oid,'rep u',cycstr(u),'v',cycstr(v),'uv',cycstr(comp(u,v)))
    uu,vv=inv(v),inv(u)
    print(' reverse_inverse orbit',oid_of_pair(uu,vv),'u',cycstr(uu),'v',cycstr(vv))
# enumerate subgroup generated by Schreier gens, starting with a modest subset then all uniq
sub={I}; qq=deque([I])
while qq:
    x=qq.popleft()
    for s in uniq:
        y=comp(x,s)
        if y not in sub:
            sub.add(y);qq.append(y)
print('centralizer subgroup enumerated size',len(sub))
# direct orbits using all centralizer elements by conjugating each factorization rep
for oid,members in enumerate(orbits):
    i=members[0];u,v=facts[i]
    direct=set()
    for c in sub:
        uu=conj(u,c)
        j=fact_lookup.get(uu)
        if j is not None:direct.add(j)
    print('direct orbit rep old oid',oid,'size',len(direct),'old oids',sorted(set(orbit[j] for j in direct)))
