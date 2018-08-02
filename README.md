# GBWT

Graph BWT is an independent implementation of the graph extension (gPBWT) of the positional Burrows-Wheeler transform (PBWT). Its initial purpose is to embed observed haplotypes in a [variation graph](https://github.com/vgteam/vg). Haplotypes are essentially sequences of nodes in the variation graph, and GBWT is best seen as the multi-string BWT of the node sequences.

The implementation uses [Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite) (SDSL). To compile, set `SDSL_DIR` in the Makefile to point to your SDSL directory. As the implementation uses C++11, OpenMP, and libstdc++ parallel mode, you need g++ 4.7 or newer to compile.

See [the wiki](https://github.com/jltsiren/gbwt/wiki) for further documentation.

## Citing GBWT

Jouni Sirén, Erik Garrison, Adam M. Novak, Benedict Paten, and Richard Durbin: **Haplotype-aware graph indexes**.
Proc. WABI 2018, LIPIcs 113, pp. 4:1-4:14, Helsinki, Finland, August 20-22, 2018.
DOI: [10.4230/LIPIcs.WABI.2018.4](https://doi.org/10.4230/LIPIcs.WABI.2018.4)

## Other references

Richard Durbin: **Efficient haplotype matching and storage using the Positional Burrows-Wheeler Transform (PBWT)**.
Bioinformatics 30(9):1266-1272, 2014.
DOI: [10.1093/bioinformatics/btu014](https://doi.org/10.1093/bioinformatics/btu014)

Adam M. Novak, Erik Garrison, and Benedict Paten: **A graph extension of the positional Burrows-Wheeler transform and its applications**.
Algorithms for Molecular Biology 12:18, 2017.
DOI: [10.1186/s13015-017-0109-9](https://doi.org/10.1186/s13015-017-0109-9)
