#!/usr/bin/env python3
"""Exact dual-pair identification for the raw representatives."""
import os, importlib.util, contextlib, io
HERE=os.path.dirname(os.path.abspath(__file__))
BASE=os.path.join(HERE,'base')
spec=importlib.util.spec_from_file_location('basecert',os.path.join(BASE,'m23_common_component_certificate.py'))
c=importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(c)
h=c.h.h
conj=lambda a: c.h.conj(h,a)
R=tuple(c.R)
Rvee=tuple(c.sig(a) for a in R)
word=c.parse_word(os.path.join(BASE,'qcompatible_interleave_word.txt'))
Uvee=Rvee
for i,d in word:
    Uvee=c.hm(Uvee,i,d)
U=tuple(c.U)
def fuse(T): return tuple(c.h.compose(T[2*i],T[2*i+1]) for i in range(5))
P=tuple(c.P); Pvee=fuse(Uvee)
checks={
 'ALL_R_ENTRIES_INVOLUTIONS': all(c.h.compose(a,a)==c.h.I() for a in R),
 'RVEE_EQUALS_H_CONJ_R': Rvee==tuple(conj(a) for a in R),
 'UVEE_EQUALS_H_CONJ_U': Uvee==tuple(conj(a) for a in U),
 'PVEE_EQUALS_H_CONJ_P': Pvee==tuple(conj(a) for a in P),
 'R_INNER_KEYS_EQUAL': c.canon(Rvee)[0]==c.canon(R)[0],
 'U_INNER_KEYS_EQUAL': c.canon(Uvee)[0]==c.canon(U)[0],
 'P_INNER_KEYS_EQUAL': c.canon(Pvee)[0]==c.canon(P)[0],
 'H_ORDER_TWO': c.h.compose(h,h)==c.h.I(),
 'H_SIGMA_H_COCYCLE': c.h.compose(h,c.sig(h))==c.h.I(),
}
for k,v in checks.items(): print(k,bool(v))
print('NO_QUADRATIC_INNER_HURWITZ_FIBER_FROM_RAW_DUALS', checks['P_INNER_KEYS_EQUAL'])
print('ALL_DUAL_COLLAPSE_CHECKS_PASS', all(checks.values()))
