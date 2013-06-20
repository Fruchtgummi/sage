r"""
Hypergraphs

This module consists in a very basic implementation of :class:`Hypergraph`,
whose only current purpose is to provide method to visualize them. This is
done at the moment through `\LaTeX` and TikZ, and can be obtained from Sage
through the ``view`` command::

    sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
    Hypergraph on 6 vertices containing 4 sets
    sage: view(H) # not tested

Note that hypergraphs are very hard to visualize, and unless it is very small
(`\leq 10` sets) or has a very specific structure (like the following one),
Sage's drawing will only bring more confusion::

    sage: g = graphs.Grid2dGraph(5,5)
    sage: sets = Set(map(Set,list(g.subgraph_search_iterator(graphs.CycleGraph(4)))))
    sage: H = Hypergraph(sets)
    sage: view(H) # not tested

.. SEEALSO::

    :class:`Hypergraph` for information on the `\LaTeX` output

Classes and methods
-------------------
"""
from sage.misc.latex import latex
from sage.sets.set import Set

class Hypergraph:
    r"""
    A hypergraph.

    A *hypergraph* `H = (V, E)` is a set of vertices `V` and a collection `E` of
    sets of vertices called *hyperedges* or edges. In particular `E \subseteq
    \mathcal{P}(V)`. If all (hyper)edges contain exactly 2 vertices, then `H` is
    a graph in the usual sense.

    .. rubric:: Latex output

    The `\LaTeX` for a hypergraph `H` is consists of the vertices set and a
    set of closed curves. The set of vertices in each closed curve represents a
    hyperedge of `H`. A vertex which is encircled by a curve but is not
    located on its boundary is **NOT** included in the corresponding set.

    The colors are picked for readability and have no other meaning.

    INPUT:

    - ``sets`` -- A list of hyperedges

    EXAMPLES::

        sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{6}]); H
        Hypergraph on 6 vertices containing 4 sets

    REFERENCES:

    - :wikipedia:`Hypergraph`
    """
    def __init__(self, sets):
        r"""
        Constructor

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Hypergraph on 6 vertices containing 4 sets
        """
        from sage.sets.set import Set
        self._sets = map(Set, sets)
        self._domain = set([])
        for s in sets:
            for i in s:
                self._domain.add(i)

    def __repr__(self):
        r"""
        Short description of ``self``.

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Hypergraph on 6 vertices containing 4 sets
        """
        return ("Hypergraph on "+str(len(self.domain()))+" "
                "vertices containing "+str(len(self._sets))+" sets")

    def edge_coloring(self):
        r"""
        Compute a proper edge-coloring.

        A proper edge-coloring is an assignment of colors to the sets of the
        hypergraph such that two sets with non-empty intersection receive
        different colors. The coloring returned minimizes the number of colors.

        OUTPUT:

        A partition of the sets into color classes.

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Hypergraph on 6 vertices containing 4 sets
            sage: C = H.edge_coloring()
            sage: C # random
            [[{3, 4, 5}], [{4, 5, 6}, {1, 2, 3}], [{2, 3, 4}]]
            sage: Set(sum(C,[])) == Set(H._sets)
            True
        """
        from sage.graphs.graph import Graph
        g = Graph([self._sets,lambda x,y : len(x&y)],loops = False)
        return g.coloring(algorithm="MILP")

    def _spring_layout(self):
        r"""
        Return a spring layout for the vertices.

        The layout is computed by creating a graph `G` on the vertices *and*
        sets of the hypergraph. Each set is then made adjacent in `G` with all
        vertices it contains before a spring layout is computed for this
        graph. The position of the vertices in the hypergraph is the position of
        the same vertices in the graph's layout.

        .. NOTE::

            This method also returns the position of the "fake" vertices,
            i.e. those representing the sets.

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Hypergraph on 6 vertices containing 4 sets
            sage: L = H._spring_layout()
            sage: L # random
            {1: (0.238, -0.926),
             2: (0.672, -0.518),
             3: (0.449, -0.225),
             4: (0.782, 0.225),
             5: (0.558, 0.518),
             6: (0.992, 0.926),
             {3, 4, 5}: (0.504, 0.173),
             {2, 3, 4}: (0.727, -0.173),
             {4, 5, 6}: (0.838, 0.617),
             {1, 2, 3}: (0.393, -0.617)}
            sage: all(v in L for v in H.domain())
            True
            sage: all(v in L for v in H._sets)
            True
        """
        from sage.graphs.graph import Graph

        g = Graph()
        for s in self._sets:
            for x in s:
                g.add_edge(s,x)

        _ = g.plot(iterations = 50000,save_pos=True)

        # The values are rounded as TikZ does not like accuracy.
        return {k:(round(x,3),round(y,3)) for k,(x,y) in g.get_pos().items()}

    def domain(self):
        r"""
        Return the set of vertices.

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Hypergraph on 6 vertices containing 4 sets
            sage: H.domain()
            set([1, 2, 3, 4, 5, 6])
        """
        return self._domain.copy()

    def _latex_(self):
        r"""
        Return a TikZ representation of the hypergraph.

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Hypergraph on 6 vertices containing 4 sets
            sage: view(H) # not tested

        With sets of size 4::

            sage: g = graphs.Grid2dGraph(5,5)
            sage: C4 = graphs.CycleGraph(4)
            sage: sets = Set(map(Set,list(g.subgraph_search_iterator(C4))))
            sage: H = Hypergraph(sets)
            sage: view(H) # not tested

        """
        from sage.rings.integer import Integer
        from sage.functions.trig import arctan2

        from sage.misc.misc import warn
        warn("\nThe hypergraph is drawn as a set of closed curves. The curve "
             "representing a set S go **THROUGH** the vertices contained "
             "in S.\n A vertex which is encircled by a curve but is not located "
             "on its boundary is **NOT** included in the corresponding set.\n"
             "\n"
             "The colors are picked for readability and have no other meaning.")

        latex.add_package_to_preamble_if_available("tikz")

        if not latex.has_file("tikz.sty"):
            raise RuntimeError("You must have TikZ installed in order "
                               "to draw a hypergraph.")

        domain = self.domain()
        pos = self._spring_layout()
        tex = "\\begin{tikzpicture}[scale=3]\n"

        colors = ["black", "red", "green", "blue", "cyan", "magenta", "yellow","pink","brown"]
        colored_sets = [(s,i) for i,S in enumerate(self.edge_coloring()) for s in S]

        # Prints each set with its color
        for s,i in colored_sets:
            current_color = colors[i%len(colors)]

            if len(s) == 2:
                s = list(s)
                tex += ("\\draw[color="+str(current_color)+","+
                        "line width=.1cm,opacity = .6] "+
                        str(pos[s[0]])+" -- "+str(pos[s[1]])+";\n")
                continue

            tex += ("\\draw[color="+str(current_color)+","
                    "line width=.1cm,opacity = .6,"
                    "line cap=round,"
                    "line join=round]"
                    "plot [smooth cycle,tension=1] coordinates {")

            # Reorders the vertices of s according to their angle with the
            # "center", i.e. the vertex representing the set s
            cx,cy = pos[s]
            s = map(lambda x:pos[x],s)
            s = sorted(s, key = lambda (x,y) : arctan2(x-cx,y-cy))

            for x in s:
                tex += str(x)+" "
            tex += "};\n"

        # Prints each vertex
        for v in domain:
            tex += "\\draw node[fill,circle,scale=.5,label={90:$"+latex(v)+"$}] at "+str(pos[v])+" {};\n"

        tex += "\\end{tikzpicture}"
        return tex

