#!/usr/bin/env python3
"""Exact order-11 boundary and clutching verification for the M23 Hurwitz component.

This script verifies:
* the order-11 power-class partition in the fixed ATLAS copy of M23;
* the 6+6 decomposition of each ordered width-11 boundary divisor;
* the component-monodromy orders for width-one and width-11 cusps;
* emptiness of the one-factor clutching bridge between those strata;
* the exact doubling from ordered width 11 to reduced width 22;
* fusion of the 11A/11B labels inside each Q'' reduction block.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, sys
from collections import Counter, deque
from pathlib import Path
from typing import Iterable, Sequence

SCRIPT=Path(__file__).resolve()
ROOT=SCRIPT.parents[1]
CODE=SCRIPT.parent
sys.path.insert(0,str(CODE))
import verify_m23_hurwitz_certificate as b
import construct_m23_reduced_jline as reduced_builder
from sympy.combinatorics import PermutationGroup
import sympy

Perm=b.Perm; Tuple4=b.Tuple4

def script_sha256(): return hashlib.sha256(SCRIPT.read_bytes()).hexdigest()

def cycles0(p: Sequence[int]) -> list[list[int]]:
    seen=[False]*len(p);out=[]
    for i in range(len(p)):
        if seen[i]:continue
        c=[];x=i
        while not seen[x]:seen[x]=True;c.append(x);x=p[x]
        out.append(c)
    return out

def point_cycles(p: Perm) -> list[list[int]]:
    seen=set();out=[]
    for i in range(1,b.N+1):
        if i in seen:continue
        c=[];x=i
        while x not in seen:seen.add(x);c.append(x);x=p[x]
        out.append(c)
    return out

def order11_conjugator_candidates(a: Perm, target: Perm) -> Iterable[Perm]:
    ca,cb=point_cycles(a),point_cycles(target)
    fa=[c[0] for c in ca if len(c)==1]; fb=[c[0] for c in cb if len(c)==1]
    aa=[c for c in ca if len(c)==11]; bb=[c for c in cb if len(c)==11]
    if len(fa)!=1 or len(fb)!=1 or len(aa)!=2 or len(bb)!=2:
        raise ValueError('expected cycle type 1*11^2')
    for sw in ((0,1),(1,0)):
        for o0 in range(11):
            for o1 in range(11):
                x=[0]*(b.N+1);x[fa[0]]=fb[0]
                for j,o in ((0,o0),(1,o1)):
                    for k in range(11):x[aa[j][k]]=bb[sw[j]][(o+k)%11]
                yield tuple(x)

def conjugator_count_in_group(a: Perm,target: Perm,G: PermutationGroup)->int:
    return sum(1 for x in order11_conjugator_candidates(a,target) if G.contains(b.to_sympy(x)))

def same_order11_class(a: Perm,target: Perm,G: PermutationGroup)->bool:
    return conjugator_count_in_group(a,target,G)>0

def build_raw_ordered_orbit(seed: Tuple4):
    _,key=b.canonical_tuple_and_key(seed)
    orbit=[seed];idx={key:0};q=deque([0])
    while q:
        i=q.popleft();t=orbit[i]
        for w in b.ALL_WORDS:
            u=b.apply_word(t,w);_,k=b.canonical_tuple_and_key(u)
            if k not in idx:
                idx[k]=len(orbit);orbit.append(u);q.append(len(orbit)-1)
    return orbit,idx

def build_operator(orbit,index,word):
    return tuple(index[b.canonical_tuple_and_key(b.apply_word(t,word))[1]] for t in orbit)

def subgroup_order(triple: Sequence[Perm])->int:
    return int(PermutationGroup([b.to_sympy(g) for g in triple]).order())

def cusp_rows(orbit,operator,prefix,G,ref11=None):
    rows=[]
    for cyc in cycles0(operator):
        rep=min(cyc);s=b.apply_word(orbit[rep],prefix);a,bb,c,d=s
        h=b.compose(a,bb);hinv=b.inverse(h)
        ho=b.permutation_order(h)
        lab=None
        if ho==11 and ref11 is not None:
            lab='11A' if same_order11_class(ref11,h,G) else '11B'
        rows.append({
            'representative_one_based':rep+1,
            'width':len(cyc),
            'node_order':ho,
            'node_label':lab,
            'left_component_group_order':subgroup_order((a,bb,hinv)),
            'right_component_group_order':subgroup_order((h,c,d)),
            'cycle_one_based':[x+1 for x in cyc],
            '_h':h,
        })
    return rows

def invop(p):
    out=[0]*len(p)
    for i,j in enumerate(p):out[j]=i
    return tuple(out)

def propagate_global_raw():
    """Rebuild the full 11,760-state braid orbit from shipped source code.

    Earlier development copies loaded cached ``reduced_state.pkl`` and
    ``reduced_ops.pkl`` files.  Those caches are unnecessary: the final
    release reconstructs the same orbit and q1,q2,q3 operators directly by
    importing ``construct_m23_reduced_jline.py``.  This keeps the certificate
    self-contained and gives the missing inputs explicit executable
    provenance.
    """
    can_orbit, index = reduced_builder.build_full_orbit()
    q1, q2, q3 = map(tuple, reduced_builder.build_q_operators(can_orbit, index))
    ops=[q1,invop(q1),q2,invop(q2),q3,invop(q3)]
    words=[b.Q1,b.Q1_INV,b.Q2,b.Q2_INV,b.Q3,b.Q3_INV]
    raw=[None]*len(can_orbit);raw[0]=b.TUPLE;dq=deque([0])
    while dq:
        i=dq.popleft();t=raw[i]
        for op,w in zip(ops,words):
            j=op[i]
            if raw[j] is None:
                raw[j]=b.apply_word(t,w);dq.append(j)
    if any(t is None for t in raw):raise RuntimeError('global raw propagation incomplete')
    return raw,q1,q2,q3

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--json',type=Path,default=ROOT/'results'/'M23_order11_boundary_certificate.json');args=parser.parse_args()
    G=PermutationGroup([b.to_sympy(b.ATLAS_A),b.to_sympy(b.ATLAS_B)])
    orbit,index=build_raw_ordered_orbit(b.TUPLE)
    B12,B13,B14=[build_operator(orbit,index,w) for w in b.FORWARD_WORDS]
    B23=build_operator(orbit,index,((1,2),))
    Q1=build_operator(orbit,index,b.Q1)
    Q1inv=invop(Q1)
    def comp(a,c): return tuple(c[a[i]] for i in range(len(a)))
    braid_conjugacy=comp(comp(Q1inv,B23),Q1)==B13

    # choose an actual order-11 node as reference
    temp=cusp_rows(orbit,B13,b.Q2,G,None)
    ref11=next(r['_h'] for r in temp if r['width']==11)
    power_counts={}
    power_labels={}
    for e in range(1,11):
        target=b.power(ref11,e);cnt=conjugator_count_in_group(ref11,target,G)
        power_counts[str(e)]=cnt
        power_labels[str(e)]='11A' if cnt else '11B'
    centralizer_order=power_counts['1']
    normalizer_order=sum(power_counts.values())

    rows13=cusp_rows(orbit,B13,b.Q2,G,ref11)
    rows14=cusp_rows(orbit,B14,b.Q3+b.Q2,G,ref11)
    def clean(rows):
        return [{k:v for k,v in r.items() if k!='_h'} for r in rows]
    def summarize(rows,width):
        sub=[r for r in rows if r['width']==width]
        return {
            'count':len(sub),
            'node_orders':dict(Counter(r['node_order'] for r in sub)),
            'node_labels':dict(Counter(r['node_label'] for r in sub if r['node_label'])),
            'component_group_order_pairs':{str(k):v for k,v in Counter((r['left_component_group_order'],r['right_component_group_order']) for r in sub).items()},
        }
    summaries={
        'B13_width1':summarize(rows13,1),'B13_width11':summarize(rows13,11),
        'B14_width1':summarize(rows14,1),'B14_width11':summarize(rows14,11),
    }
    # subgroup-order obstruction to shared three-point factor
    bridge13=set()
    for a in [r for r in rows13 if r['width']==1]:
        for z in [r for r in rows13 if r['width']==11]:
            if a['left_component_group_order']==z['left_component_group_order']:bridge13.add('left-left')
            if a['left_component_group_order']==z['right_component_group_order']:bridge13.add('left-right')
            if a['right_component_group_order']==z['left_component_group_order']:bridge13.add('right-left')
            if a['right_component_group_order']==z['right_component_group_order']:bridge13.add('right-right')
    bridge14=set()
    for a in [r for r in rows14 if r['width']==1]:
        for z in [r for r in rows14 if r['width']==11]:
            if a['left_component_group_order']==z['left_component_group_order']:bridge14.add('left-left')
            if a['left_component_group_order']==z['right_component_group_order']:bridge14.add('left-right')
            if a['right_component_group_order']==z['left_component_group_order']:bridge14.add('right-left')
            if a['right_component_group_order']==z['right_component_group_order']:bridge14.add('right-right')

    # reduced width-22 verification
    raw,q1g,q2g,q3g=propagate_global_raw()
    # Reconstruct the Q'' quotient directly from the same freshly rebuilt
    # q-operators, rather than loading a development cache.
    q_double_prime_1=reduced_builder.compose(q1g,reduced_builder.power(q3g,-1))
    shift=reduced_builder.compose(reduced_builder.compose(q1g,q2g),q3g)
    q_double_prime_2=reduced_builder.power(shift,2)
    blocks,block_of=reduced_builder.orbit_partition(
        len(raw),[q_double_prime_1,q_double_prime_2]
    )
    gamma_inf=reduced_builder.quotient_operator(q2g,blocks,block_of)
    orderings=[tuple(b.permutation_order(g) for g in t) for t in raw]
    block_orderings=[tuple(sorted(orderings[x] for x in bl)) for bl in blocks]
    color_keys=sorted(set(block_orderings));ci={k:i for i,k in enumerate(color_keys)};colors=[ci[k] for k in block_orderings]
    original=(2,2,3,5);original_color=next(i for i,k in enumerate(color_keys) if original in k)
    E22=[c for c in cycles0(gamma_inf) if len(c)==22]
    e22_label_counts=Counter();block_label_patterns=Counter();doubling_ok=True
    q2sq=tuple(q2g[q2g[i]] for i in range(len(q2g)))
    for cyc in E22:
        color_blocks=[z for z in cyc if colors[z]==original_color]
        if len(color_blocks)!=11:doubling_ok=False
        original_states=[]
        for z in color_blocks:
            ms=[x for x in blocks[z] if orderings[x]==original]
            if len(ms)!=1:doubling_ok=False;continue
            original_states.append(ms[0])
        if original_states:
            x=original_states[0];orb=[]
            while x not in orb:orb.append(x);x=q2sq[x]
            if len(orb)!=11 or set(orb)!=set(original_states):doubling_ok=False
            h=b.compose(raw[original_states[0]][1],raw[original_states[0]][2])
            lab='11A' if same_order11_class(ref11,h,G) else '11B';e22_label_counts[lab]+=1
        for z in cyc:
            labs=[]
            for x in blocks[z]:
                h=b.compose(raw[x][1],raw[x][2])
                if b.permutation_order(h)!=11:doubling_ok=False;continue
                labs.append('11A' if same_order11_class(ref11,h,G) else '11B')
            block_label_patterns[tuple(sorted(Counter(labs).items()))]+=1

    checks={
        'ordered_orbit_size_980':len(orbit)==980,
        'q1_inverse_B23_q1_equals_B13':braid_conjugacy,
        'order11_centralizer_order_11':centralizer_order==11,
        'order11_normalizer_order_55':normalizer_order==55,
        'order11_power_partition':power_labels=={'1':'11A','2':'11B','3':'11A','4':'11A','5':'11A','6':'11B','7':'11B','8':'11B','9':'11A','10':'11B'},
        'B13_width11_six_plus_six':summaries['B13_width11']['node_labels']=={'11A':6,'11B':6},
        'B14_width11_six_plus_six':summaries['B14_width11']['node_labels']=={'11A':6,'11B':6},
        'B13_width1_component_orders_12_60':summaries['B13_width1']['component_group_order_pairs']=={'(12, 60)':2},
        'B13_width11_component_orders_660_M23':summaries['B13_width11']['component_group_order_pairs']=={'(660, 10200960)':12},
        'B14_width1_component_orders_60_12':summaries['B14_width1']['component_group_order_pairs']=={'(60, 12)':2},
        'B14_width11_component_orders_M23_660':summaries['B14_width11']['component_group_order_pairs']=={'(10200960, 660)':12},
        'one_factor_bridge_empty_B13':len(bridge13)==0,
        'one_factor_bridge_empty_B14':len(bridge14)==0,
        'reduced_E22_count_12':len(E22)==12,
        'ordered_width11_to_reduced_width22_doubling':doubling_ok,
        'color_lift_E22_six_plus_six':dict(e22_label_counts)=={'11A':6,'11B':6},
        'Qprime_blocks_fuse_labels_2_plus_2':block_label_patterns==Counter({(('11A',2),('11B',2)):12*22}),
    }
    result={
        'certificate':'M23_ORDER11_BOUNDARY_AND_CLUTCHING_CERTIFICATE',
        'environment':{'python':platform.python_version(),'sympy':sympy.__version__,'script_sha256':script_sha256()},
        'order11':{'centralizer_order':centralizer_order,'normalizer_order':normalizer_order,'power_conjugator_counts':power_counts,'power_labels':power_labels,'quadratic_residues_mod_11':[1,3,4,5,9],'quadratic_nonresidues_mod_11':[2,6,7,8,10]},
        'summaries':summaries,
        'one_factor_clutching_bridge_matches':{'lambda_1':sorted(bridge13),'lambda_infinity':sorted(bridge14)},
        'reduced_width22':{'cycle_count':len(E22),'color_lift_label_counts':dict(e22_label_counts),'Qprime_block_label_patterns':{str(k):v for k,v in block_label_patterns.items()},'doubling_verified':doubling_ok,'interpretation':'the oriented order-11 label is present on an ordered color lift and is fused inside each Qprime block of the fully reduced quotient'},
        'checks':checks,
        'pass':all(checks.values()),
        'ordered_cusp_rows':{'B13':clean(rows13),'B14':clean(rows14)},
    }
    args.json.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['certificate']);print('PASS',result['pass']);print('ORDER11',result['order11']);print('SUMMARIES',summaries);print('BRIDGES',result['one_factor_clutching_bridge_matches']);print('REDUCED_WIDTH22',result['reduced_width22']);print('FAILED',[k for k,v in checks.items() if not v])
    return 0 if result['pass'] else 1

if __name__=='__main__':raise SystemExit(main())
