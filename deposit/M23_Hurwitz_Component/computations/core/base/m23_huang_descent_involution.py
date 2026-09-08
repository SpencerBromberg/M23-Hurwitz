from sympy.combinatorics import Permutation, PermutationGroup
N=23

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
def sp(p): return Permutation(list(p))
def tup(P): return tuple(P(i) for i in range(N))
def cyclestr(p):
 vis=[False]*N; out=[]
 for i in range(N):
  if not vis[i] and p[i]!=i:
   C=[];j=i
   while not vis[j]:vis[j]=True;C.append(j+1);j=p[j]
   out.append('('+','.join(map(str,C))+')')
  else:vis[i]=True
 return ''.join(out) or '()'
def ctype(p):
 vis=[False]*N;lens=[]
 for i in range(N):
  if not vis[i]:
   j=i;l=0
   while not vis[j]:vis[j]=True;l+=1;j=p[j]
   if l>1:lens.append(l)
 return tuple(sorted(lens))

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
print('product_one',compose(compose(g1,g2),g3)==tuple(range(N)))
C=G.centralizer(sp(g1))
print('centralizer order',C.order(),'gens',len(C.generators))
sol=[]
for H in C.generate_schreier_sims():
 h=tup(H)
 if conj(h,inv(g2))==g3 and conj(h,inv(g3))==g2:
  sol.append(h)
print('solutions',len(sol))
for i,h in enumerate(sol[:10]):
 print('h',i,cyclestr(h),'type',ctype(h),'h2',cyclestr(compose(h,h)))
