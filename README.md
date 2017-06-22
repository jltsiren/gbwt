# GBWT experiments

The gPBWT is an extension of the Positional Burrows-Wheeler Transform (PBWT) for graphs. Its purpose is to embed observed haplotypes in a [variation graph](https://github.com/vgteam/vg).

Haplotypes are essentially sequences of nodes in the variation graph, and the gPBWT is best seen as the multi-string BWT of the node sequences. We call this approach the Graph BWT (GBWT) to differentiate it from the earlier gPBWT.

This repository contains experiments with the construction and encoding of the GBWT.

## Better construction

* BWT for the reverse paths: build from left to right and search forward.
* Similar to RopeBWT2: load a batch of paths into memory and insert them into the existing GBWT.
* Align the paths before inserting them:
  * Rewriting a record is expensive, so we want to insert as many paths at the same time as possible.
  * Indels put paths out of sync. Different paths enter the same node at different times.
  * Construction should be faster if we resync the paths frequently.
  * Greedy alignment based on source/sink nodes of superbubbles should be enough.

## Assumptions

* We store arbitrary paths in a cyclic graph.
  * Full haplotypes in a DAG with dense node identifiers can be encoded better with the PBWT using node identifiers as positions.
  * We need a mapping from sequence id to sample id. Does vg already provide this?
* The input is an SDSL `int_vector_buffer<0>` containing sequences of node identitiers terminated by value `0`.
* How do we handle reverse complements / nodes in reverse orientation?
* Is the set of node identifiers locally dense?
  * One dense subset for the nodes in forward orientation and another dense subset for reverse orientation?
* We can support insertions and deletions if we add incoming edges (from, count) to each record.

## Record

* Incoming edges
  * Indegree
  * (from, path count) for each incoming edge
* Outgoing edges
  * Outdegree
  * (to, path rank) for each outgoing edge
* Body
  * Run-length encoding of pairs (outgoing edge rank, count)

## Encoding

* On disk, the records are stored in a single byte array.
* An index (`sd_vector`?) points to the beginning of each record.
* Runs are encoded using `Run`, while other integers are encoded using `ByteCode`.
* In-memory encoding can be the same, or we can use three `std::vector`s of pairs of integers.

## References

Richard Durbin: **Efficient haplotype matching and storage using the Positional Burrows-Wheeler Transform (PBWT)**.
Bioinformatics 30(9):1266-1272, 2014.
[DOI: 10.1093/bioinformatics/btu014](https://doi.org/10.1093/bioinformatics/btu014)

Adam M. Novak, Erik Garrison, and Benedict Paten: **A Graph Extension of the Positional Burrows-Wheeler Transform and its Applications**.
Proc. WABI 2016, Springer LNCS 9838, pp. 246-256, Aarhus, Denmark, August 22-24, 2016.
[DOI: 10.1007/978-3-319-43681-4_20](https://doi.org/10.1007/978-3-319-43681-4_20)
