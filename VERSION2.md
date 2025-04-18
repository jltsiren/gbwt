# Plans for GBWT v2.0

Version 2.0 will be a breaking change that will simplify the GBWT and align it more closely with other data models used in (human) pangenomics.
A conversion tool will be provided for old GBWT and GBZ files.

## BWT encoding

There are currently two encodings for BWT runs, depending on local alphabet size (outdegree of the node):

* If the outdegree is at most 254, the first byte encodes the edge rank of the successor node and as much of the run length as possible. Short runs are encoded in a single byte. The remaining run length is encoded using `ByteCode`.
* Otherwise we first encode the edge rank and then the run length, both using `ByteCode`.

The threshold for the first encoding can be safely increased from 254 to 256.

It may also be a good idea to add a third encoding for runs of two consecutive values (`(ab)^k`).
That will make the endmarker record in a bidirectional GBWT more compressible and faster to use.

## Locate functionality

GBWT use a rather inefficient structure for finding the sequence identifier corresponding to a GBWT position.
The [Rust implementation](https://github.com/jltsiren/gbwt-rs) does not support the locate structure.
Version 2 will remove the structure and replace it with an optional r-index (`FastLocate`), which is larger but much faster.
Removing the locate support should also simplify the construction algorithm and make it faster.

## Metadata

The current GBWT metadata model is based on haplotype paths generated from a phased VCF file.
A path name consists of four integer components: sample, contig, haplotype / phase, and fragment / phase block.
Path names are optional in the metadata, and sample / contig identifiers may be associated with string names.

Version 2 will make path, sample, and contig names mandatory if metadata is present.
There will also be higher-level metadata for fragmented sequences consisting of multiple paths.
That metadata will, for example, tell whether the fragment field stores the rank of the path fragment or its starting offset in the full sequence.
Ideally there will be corresponding changes to the [GFA specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).

We will also need data structures that support finding paths by metadata.
The current approach of scanning the path metadata does not work well with frequency-filtered graphs that may have tens of millions of paths.

## Serialization

GBWT currently supports two serialization formats: one based on [SDSL](https://github.com/vgteam/sdsl-lite) serialization, and another on [Simple-SDS](https://github.com/jltsiren/simple-sds) serialization.
The Simple-SDS format is more versatile and better documented, and it is also supported by the Rust implementation.
It will become the only format supported in version 2.
