---
title : "State-of-the-art review"
description : "Another project reports"
date : 2022-05-20
author : "Hugo Abreu, Fanny Terrier, Hugo Thomas"
---



This report presents a review of a paper written by {{< cite apers_quantum_2019
>}} in which a new quantum algorithm for graph sparsification relying on nearly
linear classical algorithms is introduced, leading to quantum speedups for
several problems such as extremal cuts and Laplacian solving.


# Introduction

$$\textit{" Graphs are nice [...], but sparse graphs are nicer. "}$$

Graphs are a very common data structure in many areas of computer
science, such as optimization and networks. Many practical problems can
indeed be reduced to graph problems, and as such are of interest to
computer scientists. Recent works, such as that by Chen *et al.*,
yielded a near-linear time classical algorithm for the exact
maximum-cost flow problem {{< cite chen_maximum_2022 >}}. It is nevertheless
possible to get an even better speedup by considering approximate
algorithms. The paper contribution is the creation of a quantum
algorithm for $\varepsilon$-spectral sparsification of graphs in time
$\tilde{O}(\frac{\sqrt{nm}}{\varepsilon})$, proving by the way the lower
bound of their algorithm. Taking into account the algorithm of Chen *et
al.*, it results in an algorithm generalizable to most graph problems.

## Graphs

Let $G = (V, E, \omega)$ be a weighted graph, where $V$ is a set of
vertices, $E$ a set of edges, and
$\omega : V \times V \rightarrow \mathbb{R}$ a weight function, with
$|V| = n$ and $|E| = m \leq \binom{n}{2}$. $G$ is said to be undirected
if for all $i,j \in V$ then $(i,j) \in E$ implies $(j,i) \in E$.

[]{#fig:weighted-graph-example label="fig:weighted-graph-example"}

The access to the graph are done via the *adjacency list*.

## Graph Laplacian

Let $G = (V, E, \omega)$ be a weighted graph, the Laplacian of $G$ is an
$n
\times n$ matrix defined as $$L_G = D - A$$ where $D$ is the degree
matrix and $A$ is the adjacency matrix, defined such that
$(D)\_{ii} = \sum_j \omega(i,j)$ and $(A)_{ij} = \omega(i,j)$. The graph
shown in Fig.
[\[fig:weighted-graph-example\]](#fig:weighted-graph-example){reference-type="ref"
reference="fig:weighted-graph-example"} has the following adjacency and
degree matrices $$A = \begin{pmatrix}
        0 & 1 & 2 & 1 \\\\\\
        1 & 0 & 2 & 2 \\\\\\
        2 & 2 & 0 & 0 \\\\\\
        1 & 2 & 0 & 0
    \end{pmatrix} \text{,}
    \quad
    D = \begin{pmatrix}
        4 & 0 & 0 & 0 \\\\\\
        0 & 5 & 0 & 0 \\\\\\
        0 & 0 & 4 & 0 \\\\\\
        0 & 0 & 0 & 3
    \end{pmatrix} \text{,}$$ which yield the Laplacian
$$L = \begin{pmatrix}
        4 & -1 & -2 & -1 \\\\\\
        -1 & 5 & -2 & -2 \\\\\\
        -2 & -2 & 4 & 0 \\\\\\
        -1 & -2 & 0 & 3
    \end{pmatrix} \text{.}$$

An interesting property of the Laplacian that arises from those
definitions is that $L_G$ is invertible if and only if the graph $G$ is
connected. Equivalently, the Laplacian of a graph can be expressed in
terms of a weighted sum of its edges Laplacian:
$$L_Q = \sum_{(i,j)\in E} \omega(i,j)L_{(i,j)}$$ where $L_{(i,j)}$
denotes the Laplacian of the edge $(i,j)$, defined as
$$L_{(i,j)} = (\boldsymbol e_i - \boldsymbol e_j)(\boldsymbol e_i - \boldsymbol e_j)^T$$
$\boldsymbol e_i$ is a unit vector with a 1 in coordinate $i$ and zeros
everywhere else. One can remark that $L_{(i,j)}$ is a very sparse
matrix, with only 4 nonzero entries.

The Laplacian is a positive semidefinite matrix i.e. the eigenvalues of
the Laplacian are non-negative.

The pseudo inverse of a Laplacian $L$ denoted $L^+$, is such that
$LL^+L = L$ and $L^+LL^+ = L^+$.

## Quadratic forms of a Laplacian

The quadratic form of a Laplacian has a number of nice properties, and
can be used to calculate quantities associated to the graph. All
quadratic forms of a Laplacian can be expressed, by linearity of the
sum, in terms of a weighted sum of its edges Laplacian:
$$\begin{aligned}
        \chi^T L_G \chi
            & = \sum_{(i,j)\in E} \omega(i,j)\chi^T L_{(i,j)}\chi\\\\\\
            & = \sum_{(i,j)\in E} \omega(i,j)(\chi(i) - \chi(j))^2
    \end{aligned}$$

An interesting example, showing how quadratic forms underlie graphs
properties, is that if $\chi_s$ is an indicator vector on
$S \subseteq V$, the quadratic form $\chi_s L_G \chi_s$ is equal to the
value of the cut $(S, S^c)$.

## Spectral Sparsification

Spectral sparsification of graphs aims to reduce the number of edges,
while keeping an approximation of interesting quantities i.e.,
approximately preserving all quadratic forms.

#### Definition
##### $\varepsilon$-sparsifier
H is an $\varepsilon$-sparsifier of G if and only if for all $\chi \in \mathbb{R}^n$n
the following holds:
$$\chi^T L_H \chi= (1 \pm \varepsilon)\chi^T L_G \chi \ .$$

Using the pseudo-inverse of the Laplacian, this definition can be
equivalently formulated as
$\chi L_H^+ = (1 \pm O(\varepsilon)) \chi L_G^+\chi$. It is also
possible to define an $\varepsilon$-spectral sparsifier taking into
account the positive semidefinite property of the Laplacian, such that
$(1-\varepsilon) L_G \preccurlyeq
L_h \preccurlyeq (1+\varepsilon) L_G$, where $\preccurlyeq$ denotes the
partial ordering on symmetric matrices. The three above definitions are
equivalent and one should use one or the other depending on the context.

#### Theorem
##### Graph Sparsifier
Every graph $G$ has an
$\varepsilon$-spectral sparsifier $H$ with a number of edges in
$\tilde{O}(\frac{n}{\varepsilon^2})$. Moreover, $H$ can be found in time
$\tilde{O}(m)$.

One should note that this is relevant only when
$\varepsilon\leq \sqrt\frac{n}{m}$

The existence of such $\varepsilon$-spectral sparsifier was proved by {{< cite
spielman_spectral_2011 >}}. Additional work of {{< cite
batson_twice-ramanujan_2012 >}} reduced the lower bound on the number of edges
in the sparsifier to $O(\frac{n}{\varepsilon^2})$.

# Classical Sparsification Algorithm {#seq:classical_sparsification}

The classical algorithm for graph sparsification is based on edge
sampling, where each edge is added to the sparsifier according to a
fixed probability distribution.

In order to be sure $L_H$ effectively approximates $L_G$, the choice of each
$p_e$ cannot be done at random. A nearly-linear classical sparsification
algorithm was introduced by {{< cite spielman_graph_2011 >}}, by approximating
effective resistance between any two edges in the graph efficiently, and thus
introducing a way to correctly sampling the edges.

## Effective resistance

In the case of unweighted graphs, the effective resistance of an edge
can be related to the connectivity of the graph: edges belonging to
strong components have a low effective resistance, and vertex cut (whose
removal renders G disconnected) tends to have a high effective
resistance.

The effective resistance of the bold edge is roughly $r_e = 1$, while
the effective resistance of the other is $r_e \in O(\frac{1}{n})$.
Spielman and Srivastava associated then the probability for an edge to
be kept while constructing the sparsifier proportionally to
$r_e \omega_e$.

In a graph $G=(V,E)$ the effective resistance $r_e$ of an edge
$e=(u,v)$, with $u,v$ two nodes, can be expressed with the quadratic
form $$\label{eq:effective-resistance}
   R_{u,v} = (\chi_u - \chi_v)^TL_G^+(\chi_u - \chi_v) \text{ ,}$$ where
$\chi_i$ is the $i^{th}$ vector of the canonical basis.
{{< cite spielman_graph_2011 >}}

## Graph spanner

The *distance* between two nodes $u$ and $v$ with respect to $G$ is
defined as:
$$\delta_G(u,v) = \min_{u\rightarrow v} \sum_{i} \omega(p_{i}, p_{i+1})^{-1},$$
which is consistent with the previous definition of effective resistance
of an edge, where $\{ u\rightarrow v\}$ is the set of all paths from $u$
to $v$ in $G$, each element $u\rightarrow v = \{p_0, \cdots, p_k\}$ is a
set of vertices of $V$. A spanner $H$ of a graph $G$ is a subgraph of
$G$ with fewer edges, where a trade-off is made between the number of
edges and the stretching of distances.

#### Definition
##### Graph spanner
An $(\alpha,\beta)$-spanner of the graph $G=(V,E)$ is a subgraph $H = (V,E_H)$
with $E_H \subseteq E$, such that $\forall u, v \in V$, $$\delta_G(u,v) \leq
\delta_H(u,v) \leq \alpha\delta_G(u,v) + \beta.$$

This definition holds for weighted graphs, in which case the weight of
the kept edges stay unchanged. In the following only multiplicative
spanners are considered, i.e. $\beta = 0$ and $\alpha = 2 \log n$,
namely, $2\log
n$-spanners. Furthermore, key objects of the algorithm for graph
$\varepsilon$-spectral sparsification described below are $r$-packings
spanners.

#### Definition
Let $G$ be a graph, an $r$-packings spanner of $G$ is an ordered set $H=(H_1,
\cdots, H_r)$ of $r$ edge-disjoint subgraphs of $G$ such that $H_i$ is a spanner
for the graph $G - \bigcup_{j=1}^{i-1} H_j$.

Koutis and Xu proposed the following algorithm {{< cite koutis_simple_2016 >}},
using the effective resistance of each edge as exhibited by Spielman and
Srivastava to construct $t$-packings spanner of the input graph:

<pre class="algorithm" style="display:none">
\begin{algorithm}
    \caption{\textbf{ClassicallySparsify}($G,\epsilon$)}
    \begin{algorithmic}
    \STATE \textit{construct an $O(\frac{\log n}{\epsilon ^2})$-packing spanner $H$ of $G$}
    \STATE $\tilde{G} \leftarrow H $
    \ForAll{$e \notin H$}
    \STATE \textit{w.p. $\frac{1}{4}$ add e to $\tilde{G}$, with weight $4\omega_e$}
    \EndFor
    \RETURN $\tilde{G}$
    \end{algorithmic}
\end{algorithm}
</pre>
and provided the following theorem:


#### Theorem {#thmClassicalSparsifier}
##### Classical sparsifier
The output $\tilde{G}$ of **ClassicallySparsify** on inputs $G$ and
$\varepsilon$ satisfies with probability $1 - \frac{1}{n^2}$ $$(1 - \varepsilon)
L_G \preccurlyeq L_{\tilde{G}} \preccurlyeq (1 + \varepsilon) L_G$$ Moreover,
the expected number of edges in $\tilde{G}$ is at most
$\tilde{O}(\frac{n}{\varepsilon^2} + \frac{m}{2})$.

The proof of the theorem of [classical sparsification]({{< relref
 "#thmClassicalSparsifier" >}}) ensures that the output of
**ClassicallySparsify** is an $\varepsilon$-spectral sparsifier. A single
iteration of the above procedure divides the number of edges in the output graph
by roughly two. Hence, repeating $t \in O(\log \frac{m}{n})$ times
**ClassicallySparsify($G, \varepsilon'$)** with $\varepsilon' \in
O(\frac{\varepsilon}{t})$ results in an $\varepsilon$-spectral sparsifier with
$\tilde{O}(\frac{n}{\varepsilon^2})$ edges.

Complexity-wise, the execution time of the provided algorithm is mostly
dominated by the construction of the
$\tilde{O}(\frac{1}{\varepsilon^2})$ spanners, each of which requires
time $\tilde{O}(m)$, giving a total time complexity of
$\tilde{O}(\frac{m}{\varepsilon^2})$.

# Quantum speed-up for graph sparsification

Apers and de Wolf propose a quantum analog to the sparsificationalgorithm described in . They build on results from classical and
quantum algorithms, in particular the classical algorithm for
sparsification by Koutis and Xu ([Algorithm
[\[alg:classical-sparsify\]](#alg:classical-sparsify){reference-type="ref"
reference="alg:classical-sparsify"}](#alg:classical-sparsify)), the
spanner algorithm by {{< cite Thorup_Zwick_2005 >}}, the
quantum algorithm for single-source shortest-path trees by
{{< cite Durr_Heiligman_2006 >}}, and an efficient
$k$-independent hash function by
{{< cite christiani_independence_2015 >}}.

## Quantum spanner algorithm

The quantum spanner algorithm proposed by {{< cite apers_quantum_2019 >}} is heavily
inspired by the best classical introduced by
{{< cite Thorup_Zwick_2005 >}}: as such, the classical algorithm will be
introduced before the quantum one.

### Classical spanner algorithm

In order to efficiently construct a graph spanner, {{< cite Thorup_Zwick_2005 >}} designed a classical algorithm based on
shortest-path trees.

#### Definition
##### Shortest-path tree
Inside a graph $G=(V,E)$, a shortest path tree $\mathcal T$, rooted at a node
$v_0$ and covering a subset $S\subseteq V$ of vertices, is a subgraph of $G$
such that for all nodes $v_S \in S$, the distance between $v_0$ and $v_S$ is the
same as in the original graph $G$, and is minimal in $G$.


Their algorithm constructs a $(2k-1)$-spanner $H$ of $G$ with
$O(k n^{1+1/k})$ edges, for some $k\in \mathbb{N}$. To do so, a family
$\{A_0,\dots, A_k\}$ of node subsets is generated at random such that
$A_0=V$, $A_k = \emptyset$ and for all $i < k$, $$\label{eq:ai-spanner}
   A_i = \{v \in A_{i-1}\ w.p.\ n^{-1/k}\} \text{ ,}$$

i.e., $A_i$ contains each edge of the previous subset with probability
$n^{-1/k}$. At each iteration $i\leq k$, for all nodes
$v \notin A_i \cap A_{i-1}$, a shortest path tree $\mathcal T(v)$
spanning the ensemble of nodes $V'$ is built from node $v$. $V'$ is
defined such that for all $v' \in V'$, the distance between $v$ and $v'$
is smaller than the distance between $v'$ and all the nodes that belongs
to $A_i$. The resulting spanner is the union of all the shortest path
trees created thereby.

### Quantum spanner algorithm

The runtime of Thorup and Zwick's algorithm is dominated by the
construction of the shortest path trees. A quantum algorithm speeding up
this construction exists {{< cite Durr_Heiligman_2006 >}}, and is strongly
inspired by Dijkstra's algorithm. In the latter, a tree $\mathcal{T}$
rooted at a node $v_0$ is recursively grown by adding the cheapest
border[^1] edge, i.e the edge $(i,j)$ such that
$$cost(i,j)= \min\{ cost(u,v) | u\in\mathcal{T},v\notin \mathcal{T}\} \ ,$$
where $cost(i,j) = \delta(v_0,i) + \frac{1}{w(i,j)}$. The quantum time
improvement arises from a speedup for the selection of the cheapest
border edge. The quantum routine called in the quantum shortest path
tree algorithm is the *minimum finding* quantum algorithm
$\text{MINFIND}(d,f,g)$, which takes as inputs

-   a *value* function $f : [N] \rightarrow \mathbb{R} \cup \{\infty\}$

-   a *type* function $g : [N] \rightarrow \mathbb{N}$

-   an integer $d \leq \dfrac N2$

and outputs a subset $I \subseteq [N]$ of size $|I| = \min\{d,M\}$,
where $M = |Im(g)|$, such that every distinct elements of $I$ have a
different type, i.e. for all $i,j \in I$ $$g(i) \neq g(j)\ ,$$ and for
$j\notin I$ and $i \in I$, having $f(j) < f(i)$ implies that there exists
an $i' \in I$ so that $$f(i')\leq f(i) \text{ and } g(i')= g(j)\ ,$$
i.e., $j$ and $i'$ have the same type.

Let $P_L$ be a subset of nodes and $E(P_L)$ the set of edges such that
$\forall (u,v) \in E(P_L)$, $u \in P_L \text{ or } v \in P_L$. In
[Algorithm [\[alg:quantum_spt\]](#alg:quantum_spt){reference-type="ref"
reference="alg:quantum_spt"}](#alg:quantum_spt), the functions $f$ and
$g$ are both defined on $E(P_L)$, in such a way that $g((u,v)) = v$ and
$$f((u,v)) = \begin{cases}
        cost(u,v)=dist(u) + \frac{1}{w(u,v)} & \text{if } u\in P_L,\, v\notin T \\\\\\
        \infty & \text{otherwise.}
    \end{cases}$$

In other words, we are looking for a subset of the border edges of the
set of nodes $P_L$ that contains at most one edge for each node in
$P_L$, and if several edges are possible, the least costly is kept. A
brief explanation follows.


<div class="algo">
$T =( V_T=\{v_0\}, E_T =\emptyset)$ $P_1 = \{v_0\}$ and $L=1$ <br>

set $\text{dist}(v_0) = 0$ and $\forall u\in V, u \neq v_0$,
$\text{dist}(u) = \infty$. <br>

$B_L = \textsc{MINFIND}(|P_L|,f,g)$ <br>
Let  $(u,v) \in B_1 \cup \dots \cup B_L$ have minimal $\text{cost}(u,v)$ with
$v\notin P_1 \cup \dots \cup P_L$. <br>
$\mathcal{T}$
$V_T \leftarrow V_T \cup \{v\}$ , $E_T \leftarrow E_T \cup \{(u,v)\}$
$\text{dist}(v) = \text{dist}(u) + 1/w(u,v)$ $P_{L+1} \leftarrow \{v\}$
, $L\leftarrow L+1$
<p> Name of the algorithm </p>
</div>


In [Algorithm
[\[alg:quantum_spt\]](#alg:quantum_spt){reference-type="ref"
reference="alg:quantum_spt"}](#alg:quantum_spt), a set of $L$ partitions
$\{P_l\}_{l=1}^L$ of the vertices covered by the shortest path tree
$\mathcal T$ is generated, and the algorithm stops only when
$\mathcal T$ covers the connected component of $v_0$.

Step 1 initializes the distances, as does Dijkstra's algorithm. In Step
4, a set $B_L$ containing the $|P_L|$ cheapest border edges with
disjoint target vertices is generated by the quantum routine
$\text{MINFIND}(|P_L|,f,g)$. Step 10 updates the distance of the
selected vertex, in a same manner as in Dijkstra's. After all the merges
of Step 13, the $P_k$ are sets of vertices of the growing tree, so that
$|P_k| = 2^{L-k}$. This ensures that since $B_k$ contains $|P_k|$ edges,
then at least one of these edges has its target outside of $\cup_{j=1}^L
P_k$ , implying in Step 5 at least one border edge exists, and is
effectively selected, thus the correctness of the algorithm (see
{{< cite apers_quantum_2019 >}} Appendix A, Proposition 5}). As a side note, at
each step $V_{\mathcal{T}}$ contains the growing tree.

#### Theorem {#thm:comp-spt}
In the worst case, [Algorithm
[\[alg:quantum_spt\]](#alg:quantum_spt){reference-type="ref"
reference="alg:quantum_spt"}](#alg:quantum_spt) returns a shortest path tree
covering the graph $G=(V,E)$ in time $\tilde{O}(\sqrt{mn})$.*

More precisely, the running time depends on the size of the connected
component in which the starting node $v_0$ is. Taking into account , one
can conclude on the overall time complexity of [Algorithm
[\[alg:quantum-spanner\]](#alg:quantum-spanner){reference-type="ref"
reference="alg:quantum-spanner"}](#alg:quantum-spanner).

#### Theorem {#thmqSpanner}
There exists a quantum algorithm that outputs in time
$\tilde{O}(k n^{1/k} \sqrt{mn})$ with high probability a
$(2k-1)$-spanner of $G$ with an expected number of edges
$O(k n^{1+1/k})$. *
:::

::: algorithm
[]{#alg:quantum-spanner label="alg:quantum-spanner"}

::: algorithmic
$A_0 = V$ ans $A_k = \emptyset$ *$H$ is initially an empty graph* *set
$A_i$ such as defined in  *
$\mathcal T \leftarrow \mathbf{QuantumSPT}(G, v)$
$H \leftarrow H \cup \mathcal T$ H
:::
:::

To conclude, setting $k= \log n + 1/2$, one can construct
$(2 \log n)$-spanners of an input graph with $n$ nodes and $m$ edges in
time $\tilde{O}(\sqrt{mn})$.

## Implicit construction of the graph though a *string*

In order to stay within a sublinear runtime, one cannot use an explicit
representation as used by Koutis and Xu in the algorithm of
[classical spartification]({{< relref
 "#thmClassicalSparsifier" >}}): indeed,
after a single iteration, the outputted graph could have up to
$\frac m2\in O(n^2)$ edges (see e.g., ).

Apers and de Wolf address this issue by constructing a random string
$r \in \{0,1\}^{m}$ encoding the discarded edges at some iteration with
0-valued bits, and later implicitly setting the corresponding weights in
the graph to 0, as shown in [Algorithm
[\[alg:quantum-sparsify\]](#alg:quantum-sparsify){reference-type="ref"
reference="alg:quantum-sparsify"}](#alg:quantum-sparsify). This enables
the construction of a spanner in the remaining graph. One can then use a
Grover search to the undiscarded $\tilde{O}(n/\varepsilon^{2})$ edges,
whose union forms the spectral sparsifier. In addition, it is possible
to further improve the classical complexity, since a $k$-query quantum
algorithm cannot distinguish a $2k$-wise independent strings from a
uniformly random one {{< cite zhandry_secure_2015 >}}.

At first, a family of independent random bit-strings
$$r_i \in \{0,1\}^{m} \text{, } \; i \in \big[\log \frac mn\big] \text{ ,}$$
is considered, such that all bits are *independent* and *equal to 1 with
probability 1/4*.

Thus the graph is represented throughout the execution with a bit-string
$r$, where each bit $b_e$ is sampled only when edge $e$ is queried.

However, thanks to the result of {{< cite zhandry_secure_2015 >}}, it is possible to discard
the random strings by considering $k$-independent hash functions, whose
[definition]({{< ref "#defkIndependentHashing" >}}) is recalled below. Hence,
such a structure allows to query for a bit-string element in time $\tilde{O}(1)$
without even having to store the bit-string, but still being able to retrieve
it. It is important to stress that it is a purely classical result.


#### Definition {#defkIndependentHashing}
##### $k$-independent hashing functions
Let $\mathcal U$ be the set of keys. A family $\mathcal{H} = \{ h : \mathcal
U \rightarrow [m]\}$ is said to be $k$-independent if for all keys $x_1,
\cdots, x_k$ in $\mathcal U$ pairwise distinct and for all values $v_1, \cdots,
v_k$ in $[m]$,
$$\big| \{ h \in \mathcal H ; \quad  h(x_1)=v_1, \cdots, h(x_k)=v_k \} \big| =
\frac{|\mathcal H |}{m^k} \ ,$$
in other words, by providing $\mathcal H$ with the uniform probability, for any
$h\in \mathcal H$,
$$\mathbb{P}\big(h(x_1)=v_1, \cdots,  h(x_k)=v_k \big) = \frac{1}{m^k} \text{
.}$$

A quantum oracle that keeps track of the weight updates is easily
constructed. Considering the $i^{th}$ iteration; given an edge $e$, let
$k$ denote the number of spanners in which $e$ appears before this
iteration.

If $k=0$, the weight of the edge $e$ is re-weighted as follows:
$$\omega_e' =
    \begin{cases}
    4^i\omega_e & \text{ if } \bigvee\limits_{l=1}^i r_l(e) = 1 \\\\\\
    0 & \text{ otherwise.}
    \end{cases}$$ Otherwise, in the case where $k>0$, the weight of the
edge $e$ is re-weighted in a different manner, so that $$\omega_e' =
    \begin{cases}
    4^{i-k}\omega_e & \text{if } \bigvee\limits_{l=1}^{j+1} r_l(e) = 1 \\\\\\
    0 & \text{otherwise,}
    \end{cases}$$ where $\vee$ is the logical disjunction.

<!-- ::: algorithm
::: algorithmic
$\forall e, \; w_e' = w_e$ and $l=\lceil \log\frac mn \rceil$
$\forall i\in[\log(m/n)], \; r_i \in \{0,1\}^m$, *create $H_i$, union of
an $O(\frac{\log^2 n}{\varepsilon^2})$-packing of spanners of
$G' = (V,E,w')$* $w_e' \leftarrow 4\,w_e'$ *use repeated Grover search
to find $\tilde{E} = \{ e \in E | w_e^{'} > 0\}$ the edges of
$\tilde{G}$* $\tilde{G}$
:::
:::
 -->



Intuitively, the unions of the
$O(\frac{\log^2 n}{\varepsilon^2})$-packing of spanners select the
*most* important edges of the graph, and the conditional reweighting
(Steps 4,5) is a way of keeping a fraction of the remaining edges in
order to *spectrally* preserve the graph (i.e., asserts that in the end
it effectively $(1+\varepsilon)$-approximates the input graph). In each
iteration, the remaining graph is classically sparsified using
[Algorithm
[\[alg:classical-sparsify\]](#alg:classical-sparsify){reference-type="ref"
reference="alg:classical-sparsify"}](#alg:classical-sparsify). The
sparsified graph is the one induced by the vertices of the initial graph
and the edges whose weight $w_e^{'}$ is greater than $0$.

#### Proposition
The probability that all $\log \frac mn$ iterations succeed
is $1 - O(\frac{\log n}{n^2})$.



*Proof* : Let $p_s$ be the probability of success and $p_f$ be the probability of
failure. If $p_s = 1-\frac{1}{n^2}$ then $p_e = \frac{1}{n^2}$, since
$\log \frac{m}{n}$ are done, the global probability of failure $P_f$ is
the sum of each $p_f$, such that
$$P_f = \frac{1}{n^2} \times \log \frac{m}{n}.$$ Since $m$ is the number
of edges of the input graph, $$m \leq \binom{n}{2}\in O(n^2),$$ thus
$$\log \frac mn \in O\big(\log \frac{n^2}{n}\big) = O(\log n)\text{ ,}$$
hence $$p_e = O\big(\frac{\log n}{n^2}\big)\text{ ,}$$ the result
follows.
:::

The overall time complexity of the algorithm depends on whether the
$(r_i)_i$ are represented with a random string. By [Definition
[4](#def:packing-spanner){reference-type="ref"
reference="def:packing-spanner"}](#def:packing-spanner), the set of
spanners is assumed ordered, allowing to binary search through the set
in time $\tilde{O}(1)$, and there is $O(i)$ calls to the aforementioned
oracle. The algorithm requires $O(\log n)$ qubits, which is the number
of qubits needed for the quantum spanner algorithm and the repeated
Grover search. In addition a QRAM [^2] of $\tilde{O}(n/\varepsilon^2)$
bits is required since the classical space complexity is dominated by
the output size, i.e. the size of the graph.

It is possible to simulate the random strings in [Algorithm
[\[alg:quantum-sparsify\]](#alg:quantum-sparsify){reference-type="ref"
reference="alg:quantum-sparsify"}](#alg:quantum-sparsify) with
$k$-independent hash functions, and hence improve the classical space
complexity from $\tilde{O}(n/\varepsilon^2)$ to
$\tilde{O}(\sqrt{mn}/\varepsilon^2)$.

Considering the efficiently computable search function
$f : E \rightarrow \{0,1\}$ such that $$f(e) = \begin{cases}
       1  & \text{if } w'_e >0 \\\\\\
       0 & \text{otherwise}
       \end{cases}
    \text{,}$$ Grover's algorithm finds a single edge in time
$\tilde{O}(\sqrt{\frac{m}{n/\varepsilon^2}})$. Therefore, retrieving
$\tilde{O}(n/\varepsilon^2)$ edges belonging to $\tilde{G}$ takes time
$\tilde{O}(\frac{\sqrt{m\,n}}{\varepsilon})$.\
As stated in , one can construct a $(\log^2 n / \varepsilon^2)$-spanner
in time $\tilde{O}(\sqrt{mn}/\varepsilon^2)$.

The overall time complexity is the sum of the runtimes needed to
simulate the random string, to construct a spanner and for the repeated
Grover search. Therefore, the total runtime is
$$2 \tilde{O}(\sqrt{mn}/\varepsilon^2) + \tilde{O}(\sqrt{mn}/\varepsilon) = \tilde{O}(\sqrt{mn}/\varepsilon^2) \ .$$

#### Theorem {#th:qu-spectral-sparsification .theorem}
##### Quantum Spectral Sparsification
The algorithm **QuantumSparsify**$(G,\varepsilon)$ returns with probability
$1-O(\log n/n^2)$ an $\varepsilon$-spectral sparsifier of $G$ with
$\tilde{O}(n/\varepsilon^2)$ edges, in time
$\tilde{O}(\sqrt{m\,n}/\varepsilon^2)$ and using a QRAM of
$\tilde{O}(\sqrt{m\,n}/\varepsilon^2)$ bits.

## Time improvement of quantum sparsification

### Approximate resistance oracle and spectral sparsification

From the result of Spielman and Srivastava {{< cite spielman_graph_2011 >}}, one can
compute a matrix $Z$ such that for all pairs $(s,t)$ of edges in $G$,
$$\label{eq:matrix-z} (1-\varepsilon) R_{s,t} \leq || Z\cdot(\chi_s-\chi_t)^2 ||
\leq (1+\varepsilon) R_{s,t}$$
in time $\tilde{O}(m/\varepsilon^2)$. $Z$ is defined as $Z =
QW^{\frac{1}{2}}BL^+$, where $L=B^TWB$ with $B$ the incidence matrix and $W$ a
diagonal matrix such that $(W)\_{ii} = \omega_{e_i}$, and $Q$ a random $\pm
1/\sqrt k$ matrix (i.e., independent Bernoulli entries).  Consequently, thanks
to , the matrix $Z$ helps $\varepsilon$-approximate the effective resistance
between any edge $e = (s,t)$ of the initial graph.

The proof of the existence of such a $Z$ matrix allows one to
efficiently create an oracle for the quantum algorithm.

#### Theorem {#th:app-resistance-oracle .theorem}
##### Sparsification with approximate resistances {{< cite spielman_graph_2011 >}}
Let $R_e/2 \leq \tilde{R_e} \leq 2R_e$ be a rough approximation of $R_e$, for
each $e\in E$ and $p_e = min(1,Cw_e\tilde{R_e}\log(n)/\varepsilon^2)$. Then,
with probability $1-1/n$, an $\varepsilon$-spectral sparsifier $H$ with
$O(n\log(n)/\varepsilon^2)$ edges can be obtained by keeping every edge $e$
independently with probability $p_e$ and rescaling its weight with $1/p_e$.

allows one to efficiently define the $\{p_e\}$ according the effective
resistance approximations $\{\tilde {R_e}\}$.

Since $H$ is an $\varepsilon$-spectral sparsifier of $G$, we have that
for all edge $e$ of $H$,
$$(1-1/\varepsilon) R^G_e \leq R^H_e \leq (1+1/\varepsilon) R^G_e \ ,$$
where $R^G_e$ and $R^H_e$ are effective resistances in $G$ and $H$,
respectively. From , the effective resistances $R_e$ can be approximated
with the matrix $Z$ in such a way that for an edge $e = (s,t)$, the
approximated resistance is
$$\tilde{R_e} = || Z \cdot (\chi_s-\chi_t)^2 || \text{ .}$$ The
probability $p_e$ of keeping an edge $e$ is taken to be proportional to
$\tilde{R_e}$, since an edge will be more important if it belongs to a
weak component, i.e. if it has a high effective resistance. Thanks to
the result of {{< cite bollobas_modern_1998 >}} Theorem 25,
$$\sum_e w_e R_e = n-1$$ for connected graphs of order $n$ [^3], and
thus, one has that $\sum_e w_e \tilde{R_e} = O(n)$. Since $\sum_e p_e$
represents the number of edges in the sparsifier, if one wants to end up
with $O(n \log n /\varepsilon^2)$ edges in the resulting graph, $p_e$
should be taken proportional to
$w_e \tilde{R_e} \log (n) / \varepsilon^2$.

In order to keep a satisfying approximation of the weights $\{w_e\}$ in
the sparsifier, we want to keep unchanged the expectation value of the
weight of each edge. Hence[^4], every weight $w_e$ is re-scaled by
$1/p_e$ i.e., $\tilde{w_e} = \frac{w_e}{p_e}$.

### Edge sampling

Classically, [Algorithm
[\[alg:classical-edge-sampling\]](#alg:classical-edge-sampling){reference-type="ref"
reference="alg:classical-edge-sampling"}](#alg:classical-edge-sampling)
shows how one can sample a subset of edges that contains every edge $e$
independently with probability $p_e$, in time
$\tilde{O}(m + \sum_e p_e)$.

::: algorithm
::: algorithmic
$S = \emptyset$ *approximate $\{p_e\}_{e\in E}$ using  * *add edge $e$
to $S$ with probability $p_e$* S
:::
:::

A quantum algorithm could sample a subset of edges more efficiently. We
assume we have access to a random string $r \in \{0,1\}^{\tilde{O}(m)}$
through the hash function $h_r: E\times [0,1] \rightarrow \{0,1\}$, such
that for all $e\in E$, $h_r(e,p_e) = 1$ with probability $p_e$ and
$h_r(e,p_e)=0$ otherwise. From this, it is possible to construct the
following oracle
$$O_s : \ket{e} \ket{v} \ket{w} \longmapsto \ket{e} \ket{v \oplus p_e} \ket{w \oplus h_r(e,p_e)} \ .$$

Due to the fact that the expected number of edges $e$ for which
$h_r(e,p_e)=1$ is $\sum_e p_e$, a repeated Grover search finds the
desired edges in time $\tilde{O}(\sqrt{m \sum_e p_e})$.

### Refined quantum sparsification algorithm

The runtime of [Algorithm
[\[alg:quantum-sparsify\]](#alg:quantum-sparsify){reference-type="ref"
reference="alg:quantum-sparsify"}](#alg:quantum-sparsify) can be
improved to $\tilde{O}(\sqrt{mn}/\varepsilon)$ by creating a first
\"rough\" $\varepsilon$-sparsifier $H$, estimating the effective
resistances of $G$ from $H$ using Laplacian solving, and then using
quantum sampling in order to sample a subset containing
$\tilde{O}(n/\varepsilon^2)$ edges.

::: algorithm
::: algorithmic
*use [Algorithm
[\[alg:quantum-sparsify\]](#alg:quantum-sparsify){reference-type="ref"
reference="alg:quantum-sparsify"}](#alg:quantum-sparsify) to construct a
$(1/100)$-spectral sparsifier $H$ of $G$* *create a
$(1/100)$-approximate resistance oracle of $H$ using , yielding
estimations $\tilde{R_e}$* *use quantum sampling to sample a subset of
the edges, keeping every edge with probability $p_e =
    min(1,Cw_e\tilde{R_e}\log(n)/\varepsilon^2)$*
:::
:::

The Step 1 of [Algorithm
[\[alg:quantum-sparsify2\]](#alg:quantum-sparsify2){reference-type="ref"
reference="alg:quantum-sparsify2"}](#alg:quantum-sparsify2) requires for
$\tilde{O}(\sqrt{mn})$ to construct the $1/100$-spectral sparsifier $H$,
in which each edge $e$ is such that its effective resistance $R_e^H$
satisfies $$(1-1/100) R^G_e \leq R^H_e \leq (1+1/100) R^G_e\ .$$
According to , there exists an oracle to derive approximated resistances
$\{\tilde{R_e}\}$ in Step 2 such that, for all edges $e=(s,t)$,
$$(1-1/100) R^H_e \leq \tilde{R^H_e} \leq (1+1/100) R^H_e \ ,$$ where
$\tilde{R^H_e} = || Z\cdot(\chi_s-\chi_t)^2 ||$. One can then deduce
that $$(1-1/100)^2 R^G_e \leq \tilde{R^H_e} \leq (1+1/100)^2 R^G_e \ .$$

Supposing that each edge $e$ is kept with probability
$$p_e = \min(1,C w_e \tilde{R^H_e} \log(n)/\varepsilon^2) \ ,$$ an
$\varepsilon$-spectral sparsifier can be constructed with
$O(n \log(n)/\varepsilon^2)$ edges according to . The approximate oracle
needed for this step requires time $\tilde{O}(n)$ to be constructed.

The quantum routine of Step 3 takes time
$\tilde{O}(\sqrt{m \sum_e p_e})$ where $$\begin{aligned}
    \sum_e p_e
        & \leq \frac{C \log(n)}{\varepsilon^2} \sum_e w_e \tilde{R^H_e} \\\\\\
        & \leq \dfrac{(1+1/100)^2 C \log(n)}{\varepsilon^2} \sum_e w_e R^G_e \ .
\end{aligned}$$ As stated in {{< cite bollobas_modern_1998 >}}, one always has
that for a connected graph of order $n$, $\sum_e w_e R^G_e = n-1$.
Therefore, we can conclude that
$\sum_e p_e \in \tilde{O}(n/\varepsilon^2)$ which implies that the total
runtime of the quantum sampling routine is
$\tilde{O}(\sqrt{mn}/\varepsilon)$.

One can notice that Step 2 only succeeds with probability $1-\frac 1n$
as claimed by . According to that, one can abort the algorithm as soon
as the runtime exceeds $\tilde{O}(\sqrt{mn}/\varepsilon)$ and start
again, yielding a runtime of
$\tilde{O}(2 \sqrt{mn}/\varepsilon) = \tilde{O}(\sqrt{mn}/\varepsilon)$
in the worst case.\
The total runtime of [Algorithm
[\[alg:quantum-sparsify2\]](#alg:quantum-sparsify2){reference-type="ref"
reference="alg:quantum-sparsify2"}](#alg:quantum-sparsify2) is the sum
of the runtimes of the three steps and it is therefore
$\tilde{O}(\sqrt{mn})$.

::: {#th:raf-spectral-sparisifaction .theorem}
#### Theorem
##### Quantum Spectral Sparsification
**RefinedQuantumSparsify**($G$,$\varepsilon)$ returns with high
probability an $\varepsilon$-spectral sparsifier $H$ with
$\tilde{O}(n/\varepsilon^2)$ edges, and has runtime
$\tilde{O}(\sqrt{mn}/\varepsilon)$. The algorithm uses $O(\log\, n)$
qubits and a QRAM of $\tilde{O}(\sqrt{mn}/\varepsilon)$.

Having arrived to , the main result of the paper was made explicit.
{{< cite apers_quantum_2019 >}}'s algorithm thus implies a quantum speedup for
solving Laplacian systems and for approximating a range of cut problems
such as min cut and sparsest cut.

This result can probably be combined with recent classical results such
as {{< cite chen_maximum_2022 >}} to yield even faster algorithms. Stay tuned\...

# QRAM Model {#ap:qram}

To achieve the speed-up promised by the quantum algorithms presented
hereby, we assume the existence of a quantum device able to run quantum
subroutines on at most $O(\log N)$ qubits, where $N$ is the size of the
problem or the input.\
Besides, we assume an access to a Quantum Random Access Memory (QRAM)
which is, as its classical analog, composed of an *input* register, a
*memory* array and an *output* register. The main variations are that
the input and output registers are composed of qubits rather than bits.
Thus, the quantum computer can address memory in superposition meaning
that a superposition of inputs returns a superposition of outputs, so
that one can design the following quantum unitary
$$\sum_j \lambda_j \ket{j}_{in} \ket{0}_{out}  \ \xrightarrow{QRAM \ access} \ \sum_j \lambda_j  \ket{j}_{in} \ket{v_j}_{out}  \ \ ,$$
where $in$, $out$ represent respectively the *input* and the *output*
registers and $v_j$ the value contained in the $j-th$ register. Hence, a
reading operation corresponds to a quantum query to the classical bits
stored in the memory array, whereas the operation of writing a bit in
the QRAM stays classical.\
Within this computational model, the complexity of an algorithm can have
several definitions. One can consider either the *time complexity*,
which counts the number of elementary gates (classical and quantum), of
quantum queries to the input and of QRAM operations, or the *query
complexity* which only counts the number of quantum queries to the
input. As an example of actual QRAM, a quantum optical implementation is
presented in {{< cite "QRAM" >}}.




[^1]: A border edge $(i,j)$ of $\mathcal T$ is so that $i\in \mathcal T$
    and $j \notin \mathcal T$.

[^2]: See for details about the QRAM model used herein.

[^3]: i.e. with $n$ edges

[^4]: Let $\tilde{W}$ be the random variable such that
    $P(\tilde{W}= \tilde{w_e}) = p_e$ and $P(\tilde{W}= {0}) = 1 - p_e$.
    Then the expectation value of $\tilde{W}$ is
    $\mathbb{E}(\tilde{W}) = p_e \tilde{w_e} + (1-p_e) \times 0  = p_e \tilde{w_e}$.


# Bibliography

{{< bibliography cited >}}

[*get the original pdf*](./resources/graph-sparsification.pdf)