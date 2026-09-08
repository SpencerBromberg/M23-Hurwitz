from collections import deque
from sympy.combinatorics import Permutation, PermutationGroup
N=23

def I(): return tuple(range(N))
def compose(p,q): return tuple(p[q[i]] for i in range(N))
def inv(p):
 r=[0]*N
 for i,j in enumerate(p): r[j]=i
 return tuple(r)
def conj(h,p): return compose(compose(h,p),inv(h))
def cyc(*cycles):
 p=list(range(N))
 for C in cycles:
  C=[x-1 for x in C]
  for a,b in zip(C,C[1:]+C[:1]): p[a]=b
 return tuple(p)
def ctype(p):
 vis=[False]*N;l=[]
 for i in range(N):
  if not vis[i]:
   j=i;n=0
   while not vis[j]:vis[j]=True;n+=1;j=p[j]
   if n>1:l.append(n)
 return tuple(sorted(l))
def cyclestr(p):
 vis=[False]*N; out=[]
 for i in range(N):
  if not vis[i] and p[i]!=i:
   c=[];j=i
   while not vis[j]:vis[j]=True;c.append(j+1);j=p[j]
   out.append('('+','.join(map(str,c))+')')
  else:vis[i]=True
 return ''.join(out) or '()'
def sp(p): return Permutation(list(p))
def tup(P): return tuple(P(i) for i in range(N))

a=cyc((1,2),(3,4),(7,8),(9,10),(13,14),(15,16),(19,20),(21,22))
b=cyc((1,16,11,3),(2,9,21,12),(4,5,8,23),(6,22,14,18),(13,20),(15,17))
G=PermutationGroup([sp(a),sp(b)])
# Provenance for the released three-point input: Huang--Jackson--Lee--
# Poonen--Pries--Zhang, arXiv:2608.08538, Section 3.  Their paper
# records a left-action convention; this implementation uses the
# composition convention defined above and verifies the displayed
# product-one and order identities exactly.  The 23A/23B names follow
# the paper's Galois-equivariant class labeling, rather than being
# inferred from element order alone.
g1=cyc((1,11),(2,23),(3,8),(4,16),(5,21),(7,20),(15,19),(18,22))
g2=cyc((1,2,11,10,16,9,6,3,23,19,20,14,21,17,4,8,22,5,18,15,13,7,12))
g3=cyc((1,2,3,4,10,11,12,7,19,18,8,6,9,16,17,21,22,5,14,20,13,15,23))
# unique descent involution found in centralizer of g1
h=cyc((1,5),(2,22),(3,8),(7,20),(10,17),(11,21),(12,14),(18,23))
assert compose(h,h)==I()
assert conj(h,g1)==g1
assert conj(h,inv(g2))==g3
assert conj(h,inv(g3))==g2

# enumerate 3A class by conjugacy orbit
q=deque([cyc((1,22,9),(2,11,8),(4,13,12),(5,23,17),(6,15,10),(18,19,21))])
cls={q[0]}; gens=[a,b,inv(a),inv(b)]
while q:
 x=q.popleft()
 for s in gens:
  y=conj(s,x)
  if y not in cls: cls.add(y);q.append(y)
assert len(cls)==56672
F2=[]; F3=[]
for u in cls:
 v=compose(inv(u),g2)
 if v in cls: F2.append((u,v))
 v3=compose(inv(u),g3)
 if v3 in cls: F3.append((u,v3))
print('F2',len(F2),'F3',len(F3))
assert len(F2)==len(F3)==138
S3=set(F3)
def sigma_pair(pair):
 u,v=pair
 # g2=uv; g2^-1=v^-1 u^-1; conjugate by h to g3
 return (conj(h,inv(v)),conj(h,inv(u)))
images=[sigma_pair(p) for p in F2]
print('sigma_bijection',len(set(images)),set(images)==S3)
# reverse map is same formula and h^2=1
back=[sigma_pair(p) for p in images]
print('sigma_squared_identity',all(back[i]==F2[i] for i in range(len(F2))))

# centralizer orbits of g2 on F2
C2=G.centralizer(sp(g2)); C3=G.centralizer(sp(g3))
print('C23 orders',C2.order(),C3.order())
F2idx={p:i for i,p in enumerate(F2)}
seen=set(); orbits=[]
C2els=[tup(x) for x in C2.generate_schreier_sims()]
for i,p in enumerate(F2):
 if i in seen: continue
 O=set()
 for c in C2els:
  z=(conj(c,p[0]),conj(c,p[1])); O.add(F2idx[z])
 seen|=O;orbits.append(sorted(O))
print('F2 centralizer orbit sizes',sorted(map(len,orbits)))
# orbit pairing under sigma against C3 orbits
F3idx={p:i for i,p in enumerate(F3)}
C3els=[tup(x) for x in C3.generate_schreier_sims()]
seen3=set();orbits3=[]; oid3={}
for i,p in enumerate(F3):
 if i in seen3: continue
 O=set()
 for c in C3els:
  z=(conj(c,p[0]),conj(c,p[1])); O.add(F3idx[z])
 k=len(orbits3)
 for j in O:oid3[j]=k
 seen3|=O;orbits3.append(sorted(O))
map_orb=[]
for O in orbits:
 j=F3idx[sigma_pair(F2[O[0]])]
 map_orb.append(oid3[j])
print('orbit pairing',map_orb,'num target orbits',len(orbits3))
# print one compatible pair
u,v=F2[orbits[0][0]]; U,V=sigma_pair((u,v))
print('g2 split u=',cyclestr(u));print('g2 split v=',cyclestr(v))
print('g3 conjugate split U=',cyclestr(U));print('g3 conjugate split V=',cyclestr(V))
print('checks',compose(u,v)==g2,compose(U,V)==g3)
# tail subgroup data for printed pair
K=PermutationGroup([sp(u),sp(v)])
print('tail subgroup order',K.order())
print('tail orbits',sorted([len(o) for o in K.orbits()],reverse=True))
# RH on each orbit for triple u,v,(uv)^-1
w=inv(g2)
for orb in sorted(K.orbits(), key=lambda O:(-len(O),min(O))):
 O=sorted(orb); S=set(O)
 def numcycles_restricted(p):
  seen=set();n=0
  for i in O:
   if i not in seen:
    n+=1;j=i
    while j not in seen:seen.add(j);j=p[j]
  return n
 d=len(O)
 inds=[d-numcycles_restricted(p) for p in (u,v,w)]
 genus=1-d+sum(inds)//2
 print('orbit',d,'indices',inds,'genus',genus)
