#!/usr/bin/env python3
from collections import deque
from pathlib import Path
import importlib.util, contextlib, io

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('m23_ten_2a_explicit', HERE / 'm23_ten_2a_explicit.py')
m = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(m)

N=23; I=m.ident(); a=m.a; b=m.b
act=[a,b,m.inv(b)]
cls={a}; q=deque([a])
while q:
    x=q.popleft()
    for s in act:
        y=m.conj(s,x)
        if y not in cls: cls.add(y); q.append(y)
print('CLASS2',len(cls))

facts=[]; idx={}
for x in cls:
    y=m.compose(m.inv(x),a)
    if y in cls:
        idx[x]=len(facts); facts.append((x,y))
print('FACTS',len(facts))

trans={a:I}; q=deque([a]); sg=[]
while q:
    x=q.popleft(); t=trans[x]
    for s in act:
        y=m.conj(s,x); ts=m.compose(s,t)
        if y not in trans:
            trans[y]=ts; q.append(y)
        else:
            c=m.compose(m.inv(trans[y]),ts)
            if c!=I and m.conj(c,a)==a: sg.append(c)
print('CLASS_ORBIT',len(trans),'RAW_CGEN',len(sg))

uniq=[]; seen=set()
for c in sg:
    if c not in seen:
        seen.add(c); seen.add(m.inv(c)); uniq.append(c)
print('UNIQ',len(uniq))

maps=[]
for c in uniq:
    arr=[]; ok=True
    for x,y in facts:
        xx=m.conj(c,x)
        if xx not in idx: ok=False; break
        arr.append(idx[xx])
    if ok: maps.append(arr)
print('MAPS',len(maps))

orb=[-1]*len(facts); ors=[]
for i in range(len(facts)):
    if orb[i]>=0: continue
    oid=len(ors); orb[i]=oid; qq=deque([i]); mem=[]
    while qq:
        j=qq.popleft(); mem.append(j)
        for A in maps:
            k=A[j]
            if orb[k]<0: orb[k]=oid; qq.append(k)
    ors.append(mem)
print('ORBIT_SIZES',sorted(map(len,ors)))

def point_orbits(gens):
    unseen=set(range(N)); out=[]
    while unseen:
        s=min(unseen); reached={s}; qq=deque([s])
        while qq:
            v=qq.popleft()
            for g in gens:
                w=g[v]
                if w not in reached: reached.add(w); qq.append(w)
        unseen-=reached; out.append(len(reached))
    return tuple(sorted(out,reverse=True))

for oid,mem in enumerate(ors):
    x,y=facts[mem[0]]
    print('ORBIT',oid,'SIZE',len(mem),'POINT_ORBITS',point_orbits([x,y]),'REP',m.cyclestr(x),m.cyclestr(y))

ok=(len(cls)==3795 and len(facts)==98 and sorted(map(len,ors))==[14,84])
print('ALL_2A2A_ORBIT_CHECKS_PASS',ok)
raise SystemExit(0 if ok else 1)
