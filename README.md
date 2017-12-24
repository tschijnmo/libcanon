[![Build Status](https://circleci.com/gh/tschijnmo/libcanon.svg?style=shield)](https://circleci.com/gh/tschijnmo/libcanon)

# libcanon

Libcanon is a C++ library for canonicalizing combinatorial objects with
lightening speed by full utilization of the isomorphism of the combinatorial
objects.

In addition to being a fully-templated generic library for the canonicalization
of any combinatorial objects, from its root in symbolic manipulation of
tensorial quantities, libcanon has direct support for the canonicalization of
strings and edge-labelled directed acyclic graphs (ELDAG).

Here a string over an alphabet is any sequence of elements in that alphabet.
By giving any string and a group of allowed permutations, libcanon is able to
get the lexicographically minimal form of the string by using only the allowed
permutations.


ELDAGs are special kind of DAGs where each edge comes with an label such that
the children of any node all has distinct labels.  This structure is a
specialization of common DAGs and is defined for the storage and manipulation
mathematical expressions with non-commutative operations.  For instance, to
store the mathematical expression of a simple function call, `f(1, 2)`, we can
have an ELDAG with nodes for `f`, `1`, and `2`, an edge from `f` to `1` with
label 0, and an edge from `f` to `2` with label 1.  Then libcanon is able to
give a unique canonical form for any ELDAG given the allowed permutations on
each internal node.  Canonicalization of normal directed and undirected graphs
can also be reduced to ELDAG canonicalized.  Actual direct support of plain
graph canonicalization is to be implemented.


Libcanon is developed by Jinmo Zhao and Prof Gustavo E Scuseria at Rice
University to be used in the symbolic tensorial algebra system
[drudge](https://github.com/tschijnmo/drudge).  It was supported as part of the
Center for the Computational Design of Functional Layered Materials, an Energy
Frontier Research Center funded by the U.S.  Department of Energy, Office of
Science, Basic Energy Sciences under Award DE-SC0012575.
