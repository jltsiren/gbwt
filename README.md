# GBWT experiments

The gPBWT is an extension of the Positional Burrows-Wheeler Transform (PBWT) for graphs. Its purpose is to embed observed haplotypes in a [variation graph](https://github.com/vgteam/vg).

Haplotypes are essentially sequences of nodes in the variation graph, and the gPBWT is best seen as the multi-string BWT of the node sequences. We call this approach the Graph BWT (GBWT) to differentiate it from the earlier gPBWT.

This repository contains experiments with the construction and encoding of the GBWT.

## References

Richard Durbin: **Efficient haplotype matching and storage using the Positional Burrows-Wheeler Transform (PBWT)**.
Bioinformatics 30(9):1266-1272, 2014.
[DOI: 10.1093/bioinformatics/btu014](https://doi.org/10.1093/bioinformatics/btu014)

Adam M. Novak, Erik Garrison, and Benedict Paten: **A Graph Extension of the Positional Burrows-Wheeler Transform and its Applications**.
Proc. WABI 2016, Springer LNCS 9838, pp. 246-256, Aarhus, Denmark, August 22-24, 2016.
[DOI: 10.1007/978-3-319-43681-4_20](https://doi.org/10.1007/978-3-319-43681-4_20)
