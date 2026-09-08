# Ten-involution and local-factorization certificates

This directory contains the exact natural-degree-23 permutation data supporting the ten-involution refinement layer.

`m23_common_component_certificate.py` verifies all of the following directly from bundled files:

- the 29,005-move full `B_10` word from the target `2A^10` refinement;
- simultaneous `2+4+4` contraction to the released `(2A,23A,23B)` tuple;
- the branch-cycle-equivariant ten-involution refinement;
- the twelve-move interleaving and five-pair fusion;
- product one and generated group order `10,200,960` for the fused five-tuple;
- the stored eight-move `B_5` word placing that fused tuple in the target `(2A^4,3A)` component;
- the order-11 collision in the displayed representative.

`m23_2a2a_to_2a_orbits.py` independently verifies that a fixed `2A` element has 98 ordered `2A*2A` factorizations and that the centralizer orbit decomposition is `84+14`.

The remaining scripts verify the released triple, the descent involution, the `3A*3A` factorization data, and the centralizer-orbit structure. Every import resolves relative to this directory.
