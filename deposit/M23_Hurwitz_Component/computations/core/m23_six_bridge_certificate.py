#!/usr/bin/env python3
from collections import deque
from pathlib import Path
import importlib.util, contextlib, io, re, hashlib
from sympy.combinatorics import PermutationGroup

HERE=Path(__file__).resolve().parent
BASE=HERE/'base'

def load(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(mod)
    return mod
m=load('m23_ten_2a_explicit',BASE/'m23_ten_2a_explicit.py')
h=load('m23_huang_ten_2a',BASE/'m23_huang_ten_2a.py')
N=23

def cyc(s):
    p=list(range(N))
    for grp in re.findall(r'\(([^()]*)\)',s):
        C=[int(x)-1 for x in re.split(r'[, ]+',grp.strip()) if x]
        for a,b in zip(C,C[1:]+C[:1]): p[a]=b
    return tuple(p)

def prod(T):
    r=m.ident()
    for x in T: r=m.compose(r,x)
    return r

def hm(T,i,d):
    T=list(T); a,b=T[i],T[i+1]
    if d==1: T[i],T[i+1]=b,m.compose(m.compose(m.inv(b),a),b)
    else: T[i],T[i+1]=m.compose(m.compose(a,b),m.inv(a)),a
    return tuple(T)

def canon(T):
    sig=[sum((1<<i) for i,p in enumerate(T) if p[v]==v) for v in range(N)]
    mn=min(sig); best=None
    for root in [v for v,z in enumerate(sig) if z==mn]:
        lab=[-1]*N; order=[root]; lab[root]=0; k=0
        while k<len(order):
            u=order[k]; k+=1
            for p in T:
                w=p[u]
                if lab[w]<0: lab[w]=len(order); order.append(w)
        if len(order)!=N: continue
        key=bytes(lab[p[u]] for p in T for u in order)
        if best is None or key<best: best=key
    return best

def dec(key): return tuple(tuple(key[i*N:(i+1)*N]) for i in range(5))
def colored_move(key,g,d):
    T=dec(key)
    if g==0:
        T=hm(T,0,d); T=hm(T,0,d)
    else: T=hm(T,g,d)
    return canon(T)

def parse_colored_word(path):
    out=[]
    for t in path.read_text().split():
        if t=='q1^2': out.append((0,1))
        elif t=='q1^-2': out.append((0,-1))
        else:
            neg=t.endswith('^-1'); q=int(t[1:].split('^')[0]); out.append((q-1,-1 if neg else 1))
    return out

# Five-branch refinement of released triple.
A=(h.g1,)+tuple(cyc(x) for x in [
'(1 2 15)(3 12 7)(4 8 18)(6 23 14)(9 10 16)(11 21 17)',
'(2 17 18)(3 6 7)(5 8 22)(9 14 11)(12 15 13)(19 20 23)',
'(1 9 16)(2 18 8)(3 7 22)(4 10 11)(5 13 19)(15 17 21)',
'(1 8 6)(2 22 19)(3 11 12)(5 14 20)(7 13 21)(15 23 16)'])

u=cyc('(1 15 6)(2 18 16)(4 20 9)(5 8 22)(7 12 21)(10 14 19)')
v=cyc('(1 6 16)(2 15 18)(4 7 21)(8 22 13)(9 20 12)(10 17 19)')
h3=cyc('(1 15)(2 19)(3 8)(4 17)(6 7)(11 23)(12 21)(14 18)')
h4=cyc('(2 7)(4 13)(6 8)(9 22)(10 11)(14 17)(15 20)(16 18)')
h5=cyc('(1 20)(3 7)(4 18)(6 19)(9 22)(10 23)(13 17)(14 16)')
h6=cyc('(3 19 14)(4 21 15)(5 7 8)(10 17 11)(12 16 20)(13 22 18)')
S=(u,v,h3,h4,h5,h6)
bridge_types=[m.ctype(x) for x in S]
bridge_order=PermutationGroup([m.sp(x) for x in S]).order()
print('SIX_BRANCH_TYPES',bridge_types)
print('SIX_BRANCH_PRODUCT_ONE',prod(S)==m.ident())
print('SIX_BRANCH_GENERATED_ORDER',bridge_order)
print('UV_RECOVERS_TARGET_SEED_FIRST_ENTRY',m.compose(u,v)==m.g1)
print('H4H5_TYPE',m.ctype(m.compose(h4,h5)))

# Complete 3A class and exact 3A*3A -> fixed 2A count.
gens=[m.a,m.b,m.inv(m.a),m.inv(m.b)]
cls3={u}; q=deque([u])
while q:
    x=q.popleft()
    for s in gens:
        y=m.conj(s,x)
        if y not in cls3: cls3.add(y); q.append(y)
count=0
for x in cls3:
    y=m.compose(m.inv(x),m.g1)
    if y in cls3: count+=1
print('CLASS3_SIZE',len(cls3))
print('3A3A_TO_2A_FACTORIZATIONS',count)

# Fuse h4,h5, then q2 q1 places the unique 2A entry first.
C=(u,v,h3,m.compose(h4,h5),h6)
C=hm(C,1,1)
C=hm(C,0,1)
print('PREP_MOVES','q2 q1')
print('PREPARED_TYPES',[m.ctype(x) for x in C])

word_path=HERE/'six_bridge_colored_word.txt'
word=parse_colored_word(word_path)
sha=hashlib.sha256(word_path.read_bytes()).hexdigest()
sha_record=(HERE/'six_bridge_colored_word.sha256').read_text().split()[0]
print('COLORED_WORD_LENGTH',len(word))
print('ORDINARY_BRAID_LETTER_LENGTH',sum(2 if g==0 else 1 for g,d in word))
print('COLORED_WORD_SHA256',sha)
print('SHA256_RECORD_MATCH',sha==sha_record)
state=canon(A)
for g,d in word: state=colored_move(state,g,d)
print('COLORED_COMPONENT_ENDPOINT_MATCH',state==canon(C))

expected_types=[(3,)*6,(3,)*6,(2,)*8,(2,)*8,(2,)*8,(3,)*6]
ok=(bridge_types==expected_types and prod(S)==m.ident() and bridge_order==10200960 and
    m.compose(u,v)==m.g1 and m.ctype(m.compose(h4,h5))==(3,)*6 and len(cls3)==56672 and count==2688 and
    [m.ctype(x) for x in C]==[(2,)*8]+[(3,)*6]*4 and len(word)==14 and sha==sha_record and state==canon(C))
print('ALL_SIX_BRANCH_CERTIFICATE_CHECKS_PASS',ok)
raise SystemExit(0 if ok else 1)
