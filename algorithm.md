# CuttedSortedAC

Our goal is to find masked superstring $S$ of a kmer set $K$,
with aim on keeping the penalty of a result the lowest possible.
As this problem is a modification of NP-hard shortest common superstring problem,
we develop a heuristic instead of an optimal algorithm.

The penalty depends on final length of $S$ and the number of runs of `1`s in the mask.
Constants for extension and starting of a new run can be customized.

We use modified Aho-Corasick tree to store $K$ and then
greedily search it for the best overlaps according to the penalty function,
until we are able to connect all kmers from $K$ into $S$.

## The CSAC

The CSAC (Cutte sorted Aho-Corasick tree)
is a subgraph of Aho-Corasick tree built over the set of kmers of length $k$
that maintains some additional properties.
The target of the structue design is to optimize memory used to store the set,
while allowing fast searching in the next phase.

The main difference between CSAC and AC tree is
that CSAC ommits links to parent nodes (which are not used in the search),
and replaces array of child nodes with an interval.
This is possible because the nodes in CSAC are sorted accoring to their depth (descending)
and lexicographically.

As a result, each node of the tree only needs to store these information during the construction:
- depth of the node
- index of the kmer, which stores prefix of the node
    (the actual prefix can be accesed in constant time thanks to knowing the depth of the node)
- index of the next node on a failure path
    (node visited by failure link from current node in the original AC tree)
- indexes of the first child of the current node
    (index of the last child can be computed in constant time by inspecting the next node)

This allows us to use only three number data fields per node
(with size depending on the size of set $N = |K|$ and length of each kmer $k$).
As each data field has enough capacity to store the index of each node of CSAC
and there can be up to $N \times k$ nodes,
we can store both the depth and kmer index in the same data field.

The leaves of CSAC are always stored in positions $0 ... N - 1$.
Thanks to this property, we only need to store their failure links
and can reuse other data fileds for storing indexes of complements and other data.

## Building the CSAC

We build the data structure by adding layers of nodes in order of decreasing depths,
starting with the leaves of depth $k$.
These correspond to the individual kmers from set, which we sort at the beginning.

Then, in each phase, we add the next layer of nodes into the array --
we traverse the nodes of the last layer and for each group with common prefix
(that contains one to four nodes) we create new parent node and set
it's first leaf index to be the first leaf index of first node in the group.

We also maintain a sorted array of kmer suffixes of corresponding lengths.
At the beginning we just copy the initial kmers.
Then, in each phase, we resort it in linear time using stable sort
according to the nucleotide index at position corresponding to current layer depth.

For each failure path starting at leaf of the tree
(and following the failure links up all the way to the root),
we maintain the last node on that path we have found so far.

If a node we are adding (representing the prefix of at least one kmer from $K$)
lies on a failure path from at least one leaf (is a suffix of kmer from $K$),
we update failure link of a last node on this failure path to point to the current node.
Then we update the last node of the failure path, also to current node.
We also merge failure paths as soon as they represent the same suffix.

This way, we are able to create the whole tree with a subset of failure links
that lie on at least one failure path from a leaf.
We don't need any other failure links, though.
All of the internal nodes also keep an interval of leaves that the form a prefix of --
we only need to store the index of the first leaf,
as index of the last leaf can be always inferred from the first-child-index of the next node.
None of the nodes need to keep index of their parent or children.

As a last step, we traverse all the leaves and for each kmer we find it's reverse complement
using binary search and store it in the corresponding leaf node.

## Searching for the masked superstring

We use modification of global greedy strategy with three parameters:
- penalty for superstring extension by base -- extension penalty
- penalty for creating a new run of ones in the mask (for inserting a block of zeros) --
run penalty
- precision -- an inverse binary logarithm of the fraction of *unexteded* leaves
    remaining at which the algorithm should stop

Again, we divide the algorithm into phases, where each phase has it's own maximal penalty,
which increases by extension penalty per phase.
We then try to connect each *unextended* leaf to an *unused* leaf, in such a way
that the leaf being *extended* has the best possible overlap with the leaf being *used*
(the fist one being a prefix, the second one being a suffix).

We use DFS search limited to steps that don't exceed the maximal penalty
and stop if there are no remaining steps to take or we find at least one *unused* leaf.
In the end we either *extend* the original leaf or leave it to be *extended* in the next phases.
We stop the whole algorithm when we run out of phases or there are less
*unextended* leaves than were defined by precision parameter.
Running out of phases happens when the penalty threshold allows us
to *extend* any leaf by visiting the root and choosing any of remaining leaves --
with zero overlap, which means the same penalty can also be achieved
by simply joining any of resulting chains together.
This happens when the penalty threshold reaches a value of
$k \times extension~penalty + run~penalty$.

From each node visited in the search (which can also be a leaf),
we try to proceed by it's failure link
and use any of the leaves in the leaf interval of the corresponding failure node.
The penalty of such step depends on the depth difference of the failure node
from the current one (with possible application of new run penalty).
If there is no leaf that can be used, we try to continue the search from all the leaves
of the failure node and the failure node itself.

We always find the best possible extension within a given maximal penalty,
as we increase the penalty threshold exactly by penalty extension in each step, and
we always try to search from the current node (by it's failure link) first --
if the optimal extension would need to start a new run
(with a possibility to have penalty which is not a multiple of extension penalty),
that would be tried as the first one.

To speed the search up,
for each internal node, we keep the index of first unused leaf of that node
and increase it if we try to connect the leaf but are not able to.
This can also happen when connecting the next leaf would create a cycle,
so the leaves skipped this way might be marked as used even if they were not.
Due to this, we reset each counter after it reaches the end of it's range.

To keep track of possible cycles, we use union-find structure with
path update. It also keeps track of beginnings of chains of connected leaves
(thanks to the asymmetrical way of connecting),
which is used when printing the resulting superstring.

To further increase the speed (at cost of increasing overall penalty),
we introduce a maximal remaining penalty counter for each leaf.
It kepps track of the maximal remaining penalty with which a leaf was visited while searching.
The leaf then cannot be visited with the same or lower remaining penalty and is skipped instead.
*TODO remove: Preventing the leaf from being visited with lower reamining penalty that it keeps
comes from an observation:
since we have visited the leaf with higher remaining penalty,
it must have been impossible to visit it with a lower remaining penalty ()...
well, technically, i guess this doesn't hold because we use it even if we don't extend
the origin leaf through this one, and still, if it made sense,
there's no point in having equality there, but well, that way it somehow works.*

To keep track of an optimal way of extending a leaf
(which can use a number of other leaves present), we use an array of backtrack indexes
and other array called backtracks.
If the way of extension contains at least one other leaf,
we append all those leaves in an reverse order into the backtracks array,
starting with the leaf to be reached (the *used* one).
We then mark the start of this chain (it's last index) in the backtrack indexes
for the leaf being *extended*. To be able to use this approach, we
also need to keep track of previous visited leaf for each leaf we visit in the algorithm
(those marks never need to be cleared as they are only used when valid).

In the end, when a sufficient amount of leaves has been *extended*,
we search through them and when we find the beginning of a chain of extending leaves
(marked by a leaf that points to itself in union-find),
we print all the leaves in the chain.
For each leaf that has a backtrack index defined,
we print the whole path to the next leaf as stored in backtracks
(we stop when we reach the next leaf in backtracks).
