# Simple-SDS serialization format

GBWT version 5, Metadata version 2. Updated 2022-01-31.

Based on Simple-SDS version 0.2.0.

## Basics

This document specifies the portable simple-sds serialization format for the GBWT.
It is an alternative to the old file format based on SDSL data structures.
The format builds upon the basic data structures described in: [https://github.com/jltsiren/simple-sds/blob/main/SERIALIZATION.md](https://github.com/jltsiren/simple-sds/blob/main/SERIALIZATION.md).

The notation uses simple-sds conventions.
Array indexes start at 0.
Rank queries count the number of occurrences strictly before the query position, as in SDSL.
A rank query past the end of the sequence counts the total number of occurrences in the sequence.
Select queries use 0-based indexing (SDSL uses 1-based indexing).
A past-the-end select query should be understood as the length of the sequence, though this is usually clarified later.

### String array

A **string array** stores multiple concatenated strings in a single array.
It reduces the overhead from memory allocations and enables fast serialization and loading.

Serialization format for string arrays:

1. `index`: Starting offset of each string as a sparse vector.
2. `alphabet`: Alphabet as a vector of bytes.
3. `strings`: Concatenated strings as an integer vector.

The serialization format uses alphabet compaction.
`alphabet` and `strings` can be decompressed into a vector of bytes as `bytes[i] = alphabet[strings[i]]`.
The sequence of bytes from `bytes[index.select(i)]` (inclusive) to `bytes[index.select(i + 1)]` (exclusive) encodes the `i`th string using UTF-8.

**Note:** For the last string, the upper bound is the length of `strings`.

**Note:** The length of `alphabet` must be equal to the number of distinct values in `strings`.

**Note:** The length of `index` is not necessarily the same as the length of `strings`.
In particular, if the last string is empty, the last value in `index` is the same as the length of `strings`.

### Dictionary

A **dictionary** stores a mapping between `n` distinct strings and a semiopen interval of integers `0..n`.

Serialization format for dictionaries:

1. `strings`: Concatenated strings as a string array.
2. `sorted_ids`: Permutation of `0..n` as an integer vector.

Permutation `sorted_ids` stores the identifiers of the strings in lexicographic order.
Strings can be mapped to their identifiers using binary search in the permutation.

### Tags

**Tags** are key-value pairs that can be used for arbitrary annotations.
Both keys and values are strings.

Serialization format for tags:

1. Keys and values as a string array.

Keys are case-insensitive, and they must be distinct.
The key of tag `i` is string `2 * i` and the value is string `2 * i + 1`.

### Byte code

Some structures encode integers using a variable-length **byte code**.
The 7 low-order bits of each byte contain data, while the high-order bit is set if and only if the encoding continues in the next byte.
The encoding is little-endian.
The first byte stores bits 0 to 6 of the integer, the second byte stores bits 7 to 13, and so on.

### Run-length encoding

**Run-length encoding**, as used in the GBWT, encodes a run of `length > 0` equal integers `value` as a sequence of bytes.
The encoding depends on the size of the **local alphabet** `sigma`.
This assumes that the local alphabet consists of integers in the semiopen interval `0..sigma`.

If `sigma < 255`, we encode short runs as a single byte.
Let `threshold = floor(256 / sigma)` be a boundary value.
If `length < threshold`, the run is encoded as byte `value + sigma * (length - 1)`.
Otherwise the first byte stores `value + sigma * (threshold - 1)`.
Subsequent bytes encode the remaining run length `length - threshold` using byte code.
If `sigma >= 255`, we first encode `value` and then `length - 1` using byte code.

**Note:** The former case would also work with `sigma == 255` and `sigma == 256`.
The `sigma < 255` condition was chosen for compatibility with older GBWT indexes.

## GBWT

The **GBWT** is a run-length encoded FM-index storing integer sequences.
The integers are interpreted as node identifiers and the sequences as paths in a graph.
Each path starts from a special **endmarker** node 0, which does not exist in the graph.
**Empty paths** containing only the endmarker are supported but discouraged.

Serialization format for the GBWT:

1. GBWT header
2. Tags
3. BWT
4. Optional document array samples
5. Optional metadata

A GBWT index may be **bidirectional**.
In a bidirectional index, the forward strand of **original node** `v` is mapped to **GBWT node** `2 * v` and the reverse strand to `2 * v + 1`.
The reverse strand of the endmarker node is never used.
A bidirectional index also stores each **original path** as two **GBWT paths**.
The forward sequence for original path `i` is stored as GBWT path `2 * i`.
The reverse sequence, which visits the opposite strands of each node in reverse order, is stored as GBWT path `2 * i + 1`.

If key `source` is present in the tags, the corresponding value indicates the GBWT implementation used for serializing the structure.
The reader may use that information for determining if it can understand the serialization formats for unspecified optional structures.
The original GBWT implementation is `jltsiren/gbwt`.

### GBWT header

**GBWT header** is a 48-byte (6-element) structure with the following fields:

1. `tag`: Value `0x6B376B37` as a 32-bit integer.
2. `version`: File format version as a 32-bit integer.
3. `sequences`: Number of integer sequences / GBWT paths stored in the BWT as an element.
4. `size`: Total length of the integer sequences (including the endmarkers) as an element.
5. `offset`: Alphabet offset as an element.
6. `alphabet_size`: Size of the alphabet as an element.
7. `flags`: Binary flags as an element.

The first two fields are 32-bit unsigned little-endian integers for compatibility with the SDSL-based serialization format.
Simple-SDS serialization format requires file format version `5`.

Because the BWT stores some information for all integers in the **effective alphabet**, the GBWT uses simple alphabet compaction based on an alphabet offset:

* Node identifier `0` becomes `0` in the effective alphabet.
* Node identifiers in the closed range from `1` to `offset` are not used.
* Node identifiers `i` in the closed range from `offset + 1` to `alphabet_size - 1` become `i - offset` in the effective alphabet.

In a bidirectional GBWT, value `1` in the effective alphabet typically corresponds to the forward strand of the original node with the smallest identifier.

The following flags are supported:

* `0x0001`: The GBWT is bidirectional.
* `0x0002`: The metadata structure is present.
* `0x0004`: Simple-SDS format.

Other flag bits must not be set.
The metadata flag must be set if and only if the metadata structure is present.
If the simple-sds bit is not set, the serialized data is in the SDSL format.

### BWT

The **BWT** stores graph topology and the paths.
It is divided into records that correspond to GBWT nodes.

Serialization format for the BWT:

1. `index`: Starting offset of each record as a sparse vector.
2. `data`: Encoded records as a vector of bytes.

Let `v` be a GBWT node and `i` the corresponding value in the effective alphabet.
The corresponding record is encoded as the sequence of bytes from `data[index.select(i)]` (inclusive) to `data[index.select(i + 1)]` (exclusive).
Local alphabet size `sigma` is the number of outgoing edges from node `v`.
If the GBWT is bidirectional and original graph is a bidirected sequence graph, `sigma` is the number of edges adjacent to the exit side of the node traversal corresponding to `v`.
Edges not used on the indexed paths may be omitted from the local alphabet.

**Note:** For the last record, the upper bound is the length of `data`.

**Note:** The length of `index` must be equal to the length of `data`.

**Note:** If the BWT encodes the topology of the subgraph induced by the paths, all unused nodes and edges must be omitted.

### BWT records

The **record**  for GBWT node `v` starts with a header that encodes graph topology.
It starts with the byte code encoding of local alphabet size `sigma`.
Then, for each outgoing edge `(v, w)` in ascending order by `w`, we store a pair of byte code encoded values.
If a path ends at node `v`, the first edge is always `(v, 0)` to the endmarker.
The first value is `w - prev`, where `prev` is either the destination node of the previous edge or `0` for the first edge.
The second value is `rank(v, w)`, or the number of occurrences of `w` in the BWT before node `v`.
This can be understood as the total number of times an edge `(u, w)` is used on the indexed paths, for all `u < v`.

The header is followed by the body of the record.
Each offset `BWT(v)[j]` in the body corresponds to a visit by a path to node `v`.
The visits are sorted by the previous node `u` on the path, with ties broken by the order in node `u`.
`BWT(v)[j]` stores the next node `w` visited by the path, or the endmarker `0` if the path ends.
We encode the body by replacing each node `w` with the index of the corresponding edge in the header and by run-length encoding the result.
A path that visits offset `j` of node `v` continues to offset `rank(v, w) + BWT(v).rank(j, w)` of node `w = BWT(v)[j]`.

**Note:** If a node does not exist, the corresponding record is encoded with `sigma` set to `0`.

**Note:** The GBWT path with identifier `j` starts at offset `j` of the endmarker `0`.

**Note:** In a bidirectional GBWT, run-length encoding does not compress the body of the endmarker well.
Implementations may want to decompress the endmarker into the array of starting positions of each path.

**Note:** The GBWT is a multi-string BWT.
Each occurrence of the endmarker marking the end of a path in the GBWT has a distinct implicit character value.
Hence we cannot jump from the end of a path to its start based on the information stored in the BWT.

### Document array samples

**Document array samples** can be used for converting BWT positions to GBWT path identifiers and for jumping to a specified offset in a specified sequence.
These structures are often compilicated and both implementation-dependent and application-dependent.
Hence they are optional and not part of this specification.

## Metadata

**Metadata** stores statistics about the paths and structured names for each original path.

Serialization format for metadata:

1. Metadata header.
2. Vector of path names.
3. Sample names as a dictionary.
4. Contig names as a dictionary.

The number of path names must be either `0` or equal to the number of original paths in the GBWT index.
The number of sample/contig names must be either `0` or equal to the sample / contig count in the header.

**Note:** The absence of names is indicated by an empty structure instead of an absent optional structure.
An optional structure implies that the data is optional, while the names are the core of the metadata.

### Metadata header

**Metadata header** is a 40-byte (5-element) structure with the following fields:

1. `tag`: Value `0x6B375E7A` as a 32-bit integer.
2. `version`: File format version as a 32-bit integer.
3. `sample_count`: Number of samples as an element.
4. `haplotype_count`: Number of haplotypes as an element.
5. `contig_count`: Number of contigs as an element.
6. `flags`: Binary flags as an element.

The first two fields are 32-bit unsigned little-endian integers for compatibility with the SDSL-based serialization format.
Simple-SDS serialization format requires file format version `2`.
Contigs typically match connected components in the graph.
Haplotype count is an estimate for the total number of full-length paths in each component.

The following flags are supported:

* `0x0001`: Metadata contains path names.
* `0x0002`: Metadata contains sample names.
* `0x0004`: Metadata contains contig names.

Other flag bits must not be set.
A name flag must be set if and only if the corresponding names are present in the structure.

### Path names

**Path name** is a 16-byte (2-element) structure with the following fields:

1. `sample`: Sample identifier as a 32-bit integer.
2. `contig`: Contig identifier as a 32-bit integer.
3. `phase`: Phase / haplotype identifier as a 32-bit integer.
4. `fragment`: Fragment index as a 32-bit integer.

The fields use 32-bit unsigned little-endian integers to save space.
Each path name must be unique.
Sample / contig identifiers must be in the semiopen interval `0..sample_count` / `0..contig_count`.
If sample / contig names are present, the corresponding dictionaries can be used for mapping between identifiers and names.

The `fragment` field can be used for several purposes, including:

* Fragment index for path fragments from the same sample, contig, and phase.
* Starting offset of the fragment in the corresponding sequence.
