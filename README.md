# GBWT

The gPBWT is an extension of the Positional Burrows-Wheeler Transform (PBWT) for graphs. Its purpose is to embed observed haplotypes in a [variation graph](https://github.com/vgteam/vg).

Haplotypes are essentially sequences of nodes in the variation graph, and the gPBWT is best seen as the multi-string BWT of the node sequences. We call this approach the Graph BWT (GBWT) to differentiate it from the earlier gPBWT.

This repository should eventually become a scalable GBWT implementation.

## Assumptions

* We store arbitrary paths in a cyclic graph.
  * Full haplotypes in a DAG with dense node identifiers can be encoded better with the PBWT using node identifiers as positions.
  * We need a mapping from sequence id to sample id in vg.
* The input is an SDSL `int_vector<0>` containing sequences of node identitiers terminated by value `0`.
  * All sequences from the input are inserted simultaneously into the existing index.
* The set of node identifiers in a chromosome is locally dense.
  * The identifiers are dense in a range *[a,b]* containing the identifiers.
  * The low-order bit tells the orientation, but we can ignore that in GBWT.
* We build BWT for the reverse sequences, as in PBWT.
  * As a result, we support forward searching instead of backward searching.
  * The construction proceeds forward in the sequences.

## Record

* The BWT has one record for each node in the alphabet.
  * We do not store the records in the maximal range *[1, offset]* of nodes not occurring in the input.
  * This makes single-chromosome indexes for multi-chromosome graphs more space-efficient.
* Incoming edges
  * (from, path count) for each incoming edge in sorted order
  * Only in the dynamic record
* Outgoing edges
  * (to, path rank) for each outgoing edge in sorted order
* Body
  * Run-length encoding of pairs (outgoing edge rank, count)

## Encoding

* On disk, the records are stored in a single byte array.
* An index (`sd_vector`) points to the beginning of each record.
* Runs are encoded using `Run`, while other integers are encoded using `ByteCode`.
* The Static in-memory encoding is the same as on disk.
* The dynamic encoding required for construction uses three `std::vector`s of pairs of integers.

## TODO

* Construction of compressed/dynamic GBWT from the other.
* Special case for merging when the node ids do not overlap.
* Query interface.
* `locate()` support for determining sequence identifiers.
* Memory mapped compressed GBWT.

## References

Richard Durbin: **Efficient haplotype matching and storage using the Positional Burrows-Wheeler Transform (PBWT)**.
Bioinformatics 30(9):1266-1272, 2014.
[DOI: 10.1093/bioinformatics/btu014](https://doi.org/10.1093/bioinformatics/btu014)

Adam M. Novak, Erik Garrison, and Benedict Paten: **A Graph Extension of the Positional Burrows-Wheeler Transform and its Applications**.
Proc. WABI 2016, Springer LNCS 9838, pp. 246-256, Aarhus, Denmark, August 22-24, 2016.
[DOI: 10.1007/978-3-319-43681-4_20](https://doi.org/10.1007/978-3-319-43681-4_20)
