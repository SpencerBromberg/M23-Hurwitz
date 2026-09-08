from collections import deque
from sympy.combinatorics import Permutation, PermutationGroup

N=23

def ident(): return tuple(range(N))
def compose(p,q): # p*q = apply q then p
    return tuple(p[q[i]] for i in range(N))
def inv(p):
    r=[0]*N
    for i,j in enumerate(p): r[j]=i
    return tuple(r)
def conj(g,p): # g p g^-1
    return compose(compose(g,p),inv(g))
def cyc(*cycles):
    p=list(range(N))
    for C in cycles:
        C=[x-1 for x in C]
        for a,b in zip(C,C[1:]+C[:1]): p[a]=b
    return tuple(p)
def cyclestr(p):
    vis=[False]*N; out=[]
    for i in range(N):
        if not vis[i] and p[i]!=i:
            C=[]; j=i
            while not vis[j]:
                vis[j]=True; C.append(j+1); j=p[j]
            out.append('('+','.join(map(str,C))+')')
        else: vis[i]=True
    return ''.join(out) or '()'
def ctype(p):
    vis=[False]*N; lens=[]
    for i in range(N):
        if not vis[i]:
            j=i;l=0
            while not vis[j]: vis[j]=True;l+=1;j=p[j]
            if l>1:lens.append(l)
    return tuple(sorted(lens))
def sp(p): return Permutation(list(p))

# Official ATLAS standard generators (natural degree-23 representation)
a=cyc((1,2),(3,4),(7,8),(9,10),(13,14),(15,16),(19,20),(21,22))
b=cyc((1,16,11,3),(2,9,21,12),(4,5,8,23),(6,22,14,18),(13,20),(15,17))

# Printed five-branch seed (2A,2A,2A,2A,3A)
g1=cyc((2,6),(4,12),(5,8),(10,17),(13,22),(14,19),(15,16),(20,21))
g2=cyc((1,15),(2,19),(3,8),(4,17),(6,7),(11,23),(12,21),(14,18))
g3=cyc((2,7),(4,13),(6,8),(9,22),(10,11),(14,17),(15,20),(16,18))
g4=cyc((1,20),(3,7),(4,18),(6,19),(9,22),(10,23),(13,17),(14,16))
g5=cyc((3,19,14),(4,21,15),(5,7,8),(10,17,11),(12,16,20),(13,22,18))
seed=[g1,g2,g3,g4,g5]

# Known 2A*2A=2A factorization from exact search.
x0=a
y0=cyc((1,10),(2,9),(3,4),(6,17),(7,8),(11,23),(13,21),(14,22))
h0=compose(x0,y0)
assert ctype(x0)==ctype(y0)==ctype(h0)==(2,)*8

# Traverse the conjugacy orbit of h0 and retain one conjugator c with c h0 c^-1 = target.
def conjugators_to_targets(targets):
    gens=[a,b,inv(a),inv(b)]
    targetset=set(targets)
    found={}
    q=deque([(h0, ident())]); seen={h0}
    while q and len(found)<len(targetset):
        h,c=q.popleft()
        if h in targetset and h not in found: found[h]=c
        for s in gens:
            hn=conj(s,h)
            if hn not in seen:
                # if h=c h0 c^-1, then s h s^-1=(s c)h0(s c)^-1
                cn=compose(s,c)
                seen.add(hn); q.append((hn,cn))
    return found, len(seen)

found, orbit_seen = conjugators_to_targets([g1,g2,g3,g4])
assert len(found)==4
pairs=[]
for gi in [g1,g2,g3,g4]:
    c=found[gi]
    x=conj(c,x0); y=conj(c,y0)
    assert compose(x,y)==gi and ctype(x)==ctype(y)==(2,)*8
    pairs.append((x,y))

# Enumerate M23 once to obtain its complete 2A class; then find a 2A*2A factorization of g5.
G=PermutationGroup([sp(a),sp(b)])
# use conjugacy orbit under generators, cheaper than group enumeration
q=deque([a]); cls2={a}; gens=[a,b,inv(a),inv(b)]
while q:
    h=q.popleft()
    for s in gens:
        hn=conj(s,h)
        if hn not in cls2:
            cls2.add(hn); q.append(hn)
assert len(cls2)==3795
zpair=None
for x in cls2:
    y=compose(inv(x),g5) # x*y=g5
    if y in cls2:
        zpair=(x,y); break
assert zpair is not None
z1,z2=zpair
assert compose(z1,z2)==g5

T=[]
for x,y in pairs: T.extend([x,y])
T.extend([z1,z2])
assert len(T)==10 and all(ctype(t)==(2,)*8 for t in T)
prod=ident()
for t in T: prod=compose(prod,t)
assert prod==ident()

H=PermutationGroup([sp(t) for t in T])
ordH=H.order()
assert ordH==10200960

# Genus of degree-23 source for ten 2A branch points.
# each 2A index is 8: 2g-2=-46+80=34, so g=18.
print('ATLAS_GROUP_ORDER',G.order())
print('2A_CLASS_SIZE',len(cls2))
print('CONJUGACY_ORBIT_VISITED',orbit_seen)
print('SEED_TYPES',[ctype(g) for g in seed])
print('SEED_PRODUCT_ONE', all(compose(compose(compose(compose(g1,g2),g3),g4),g5)[i]==i for i in range(N)))
print('TEN_2A_COUNT',len(T))
print('TEN_2A_PRODUCT_ONE',prod==ident())
print('TEN_2A_GENERATED_ORDER',ordH)
print('TEN_2A_SOURCE_GENUS',18)
for i,t in enumerate(T,1): print(f't{i} = {cyclestr(t)}')
print('COLLAPSE_PAIRS')
for i,(x,y) in enumerate(pairs,1): print(i, cyclestr(compose(x,y)), '== g'+str(i))
print('5', cyclestr(compose(z1,z2)), '== g5')
