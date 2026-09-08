#!/usr/bin/env python3
from __future__ import annotations
from itertools import combinations
import sympy as sp
P=211; N=17
x=sp.symbols('x')

def pad(a,n): return a+[0]*(n-len(a))
def add(a,b):
    n=max(len(a),len(b));a=pad(a,n);b=pad(b,n);return [a[i]+b[i] for i in range(n)]
def neg(a): return [-z for z in a]
def sub(a,b): return add(a,neg(b))
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b): out[i+j]+=u*v
    return out
def scale(a,s): return [s*u for u in a]
def power(a,n):
    out=[1];base=a[:]
    while n:
        if n&1: out=mul(out,base)
        base=mul(base,base);n//=2
    return out
def deriv(a): return [i*a[i] for i in range(1,len(a))]
def shift(a,n=1): return [0]*n+a

def phi(u):
    a1,a0,b0,l1,l0,m2,m1,m0,v7,v6,v5,v4,v3,v2,v1,v0,q=u
    A=[a0,a1,1];B=[b0,1];L=[l0,l1,1];M=[m0,m1,m2,1];V=[v0,v1,v2,v3,v4,v5,v6,v7]
    core=sub(scale(shift(add(mul(deriv(L),M),scale(mul(L,deriv(M)),3))),2),mul(L,M))
    R=sub(mul(mul(A,B),core),scale(shift(mul(mul(L,M),add(scale(mul(deriv(A),B),5),mul(A,deriv(B))))),2))
    F=sub(sub(scale(mul(power(L,2),power(M,6)),4*q),shift(mul(power(A,10),power(B,2)))),mul(power(R,2),V))
    return pad(F,24)[:24]

class AD:
    __slots__=('v','g')
    def __init__(self,v,g=None): self.v=v;self.g=[0]*N if g is None else g
    def __add__(self,o):
        o=toad(o);return AD(self.v+o.v,[a+b for a,b in zip(self.g,o.g)])
    __radd__=__add__
    def __neg__(self):return AD(-self.v,[-a for a in self.g])
    def __sub__(self,o):return self+(-toad(o))
    def __rsub__(self,o):return toad(o)-self
    def __mul__(self,o):
        o=toad(o);return AD(self.v*o.v,[self.v*b+o.v*a for a,b in zip(self.g,o.g)])
    __rmul__=__mul__
def toad(z):return z if isinstance(z,AD) else AD(z)

def rref(A,p=P):
    A=[[int(v)%p for v in row] for row in A];m=len(A);n=len(A[0]);r=0;piv=[]
    for c in range(n):
        k=next((i for i in range(r,m) if A[i][c]),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r];iv=pow(A[r][c],-1,p);A[r]=[(v*iv)%p for v in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                f=A[i][c];A[i]=[(A[i][j]-f*A[r][j])%p for j in range(n)]
        piv.append(c);r+=1
        if r==m:break
    return A,piv

def solve_affine(A,b,p=P):
    m=len(A);n=len(A[0]);M=[[A[i][j]%p for j in range(n)]+[b[i]%p] for i in range(m)];r=0;piv=[]
    for c in range(n):
        k=next((i for i in range(r,m) if M[i][c]),None)
        if k is None:continue
        M[r],M[k]=M[k],M[r];iv=pow(M[r][c],-1,p);M[r]=[(v*iv)%p for v in M[r]]
        for i in range(m):
            if i!=r and M[i][c]:
                f=M[i][c];M[i]=[(M[i][j]-f*M[r][j])%p for j in range(n+1)]
        piv.append(c);r+=1
    for i in range(r,m):
        if all(M[i][j]==0 for j in range(n)) and M[i][n]:return None,None,None
    free=[j for j in range(n) if j not in piv];part=[0]*n
    for i,c in enumerate(piv):part[c]=M[i][n]
    basis=[]
    for f in free:
        z=[0]*n;z[f]=1
        for i,c in enumerate(piv):z[c]=(-M[i][f])%p
        basis.append(z)
    return part,basis,piv

# Seed construction.
Aexpr=(x-1)*(x-9);Bexpr=x-1;Lexpr=(x-9)**2;q0=95
m=sp.symbols('m');Mexpr=(x-1)*(x-9)*(x-m)
Rexpr=sp.expand(Aexpr*Bexpr*(2*x*(sp.diff(Lexpr,x)*Mexpr+3*Lexpr*sp.diff(Mexpr,x))-Lexpr*Mexpr)-2*x*Lexpr*Mexpr*(5*sp.diff(Aexpr,x)*Bexpr+Aexpr*sp.diff(Bexpr,x)))
Wexpr=sp.expand(4*q0*Lexpr**2*Mexpr**6-x*Aexpr**10*Bexpr**2)
solutions=[]
for mv in range(P):
    Rv=sp.Poly(Rexpr.subs(m,mv),x,modulus=P);Wv=sp.Poly(Wexpr.subs(m,mv),x,modulus=P)
    if Rv.is_zero:continue
    Q,rem=sp.div(Wv,Rv**2,domain=sp.GF(P))
    if rem.is_zero:solutions.append((mv,Rv,Q))
assert [s[0] for s in solutions]==[1,179]
m3=179;Rpoly=solutions[1][1];Vpoly=solutions[1][2]
Vcoeff=[int(Vpoly.nth(i))%P for i in range(8)]
u0=[201,9,210,193,81,22,111,77]+list(reversed(Vcoeff))+[q0]
assert all(v%P==0 for v in phi(u0))

uad=[]
for i,v in enumerate(u0):
    g=[0]*N;g[i]=1;uad.append(AD(v,g))
Fad=phi(uad)
J=[[int(e.g[j])%P for j in range(N)] for e in Fad]
rank=len(rref(J)[1]);assert rank==12
F0=phi(u0);rhs1=[(-(v//P))%P for v in F0]
part,K,piv=solve_affine(J,rhs1);assert part is not None and len(K)==5
JT=[list(col) for col in zip(*J)];_,Lbasis,_=solve_affine(JT,[0]*N);assert len(Lbasis)==12

def kcomb(a):return [sum(a[i]*K[i][j] for i in range(5))%P for j in range(N)]
def obs(v):return [sum(l[i]*v[i] for i in range(24))%P for l in Lbasis]
def z1_of(alpha):return [(part[j]+kcomb(alpha)[j])%P for j in range(N)]
def alpha_from_u(u4):
    a1,a2,a3,a4=u4
    return [(113-14*a1-124*a2-155*a3-91*a4)%P,a1,a2,a3,a4]
def stage2_obstruction(alpha):
    u=[u0[j]+P*z1_of(alpha)[j] for j in range(N)]
    F=phi(u);assert all(v%(P**2)==0 for v in F)
    return obs([(v//(P**2))%P for v in F])

def make_u1(alpha):
    return [u0[j]+P*z1_of(alpha)[j] for j in range(N)]

def canonical_step(u,n):
    """Given F(u)=0 mod p^n, choose free-zero digit and return u' mod p^(n+1), digit, kernel; None if insoluble."""
    F=phi(u);mod=P**n
    assert all(v%mod==0 for v in F)
    rhs=[(-(v//mod))%P for v in F]
    z,Kcur,_=solve_affine(J,rhs)
    if z is None:return None,None,None
    un=[u[j]+mod*z[j] for j in range(N)]
    assert all(v%(mod*P)==0 for v in phi(un))
    return un,z,Kcur

def canonical_to_p3(alpha):
    u1=make_u1(alpha)
    return canonical_step(u1,2)

def stage4_obstruction(u4):
    alpha=alpha_from_u(u4)
    if stage2_obstruction(alpha)!=[0]*12:
        raise ValueError('not in stage-two hyperplane')
    u1=make_u1(alpha)
    u2,z2,K2=canonical_step(u1,2)
    assert u2 is not None and K2==K
    F=phi(u2);assert all(v%(P**3)==0 for v in F)
    return obs([(v//(P**3))%P for v in F])

def lift_to_p4(u4):
    alpha=alpha_from_u(u4);u1=make_u1(alpha)
    u2,z2,K2=canonical_step(u1,2)
    if u2 is None:return None
    u3,z3,K3=canonical_step(u2,3)
    if u3 is None:return None
    return {'alpha':alpha,'u1':u1,'u2':u2,'u3':u3,'z2':z2,'z3':z3}

def next_obstruction_with_kernel(u,n,base_digit,gamma):
    """Obstruction at p^(n+1) after digit base_digit+K gamma is inserted at p^n."""
    digit=[(base_digit[j]+kcomb(gamma)[j])%P for j in range(N)]
    uu=[u[j]+(P**n)*digit[j] for j in range(N)]
    F=phi(uu);assert all(v%(P**(n+1))==0 for v in F)
    return obs([(v//(P**(n+1)))%P for v in F]),uu,digit

def lookahead(u,n):
    """u solves mod p^n. Solve current digit, then compute affine next obstruction in kernel choices.
    Returns part digit, c, A (12x5), and if solvable a gamma plus u_next solving p^(n+1)."""
    F=phi(u);assert all(v%(P**n)==0 for v in F)
    rhs=[(-(v//(P**n)))%P for v in F]
    z,Kcur,_=solve_affine(J,rhs)
    if z is None:return {'current_solvable':False}
    assert Kcur==K
    c,uu0,_=next_obstruction_with_kernel(u,n,z,[0]*5)
    cols=[]
    for i in range(5):
        g=[0]*5;g[i]=1
        vi,_,_=next_obstruction_with_kernel(u,n,z,g)
        cols.append([(vi[r]-c[r])%P for r in range(12)])
    A=[[cols[j][i] for j in range(5)] for i in range(12)]
    gamma,_,_=solve_affine(A,[(-v)%P for v in c])
    result={'current_solvable':True,'z':z,'c':c,'A':A,'A_rank':len(rref(A)[1]),'next_solvable':gamma is not None,'gamma':gamma}
    if gamma is not None:
        digit=[(z[j]+kcomb(gamma)[j])%P for j in range(N)]
        un=[u[j]+(P**n)*digit[j] for j in range(N)]
        assert all(v%(P**(n+1))==0 for v in phi(un))
        # gamma chosen ensures obstruction for following step zero, hence next current equation is solvable.
        result['digit']=digit;result['u_next']=un
    return result

def adaptive_branch(u4,max_n=6):
    """Start at u1 solving p^2. At each n choose digit whose lookahead makes p^(n+2) possible.
    Returns highest constructed precision and failure to next precision."""
    alpha=alpha_from_u(u4);u=make_u1(alpha);n=2;history=[]
    while n<=max_n:
        la=lookahead(u,n);history.append({'n':n,'A_rank':la.get('A_rank'),'c':la.get('c'),'next_solvable':la.get('next_solvable')})
        if not la.get('current_solvable',False):
            return {'alpha':alpha,'highest_precision':n,'failure':'current','history':history,'u':u}
        if not la.get('next_solvable',False):
            # current digit z creates a solution mod p^(n+1), but no kernel choice makes following obstruction vanish.
            ucurr=[u[j]+(P**n)*la['z'][j] for j in range(N)]
            assert all(v%(P**(n+1))==0 for v in phi(ucurr))
            return {'alpha':alpha,'highest_precision':n+1,'failure_precision':n+2,'failure':'lookahead','history':history,'u':ucurr,'failure_c':la['c'],'failure_A':la['A']}
        u=la['u_next'];n+=1
    return {'alpha':alpha,'highest_precision':n,'failure':None,'history':history,'u':u}
