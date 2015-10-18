class DeBruijnGraph:
    @staticmethod
    def chop(st, k):
        for i in xrange(0, len(st)-(k-1)):
            yield (st[i:i+k], st[i:i+k-1], st[i+1:i+k])

    class Node:
        def __init__(self, km1mer):
            self.km1mer = km1mer
            self.nin = 0
            self.nout = 0

        def isSemiBalanced(self):
            return abs(self.nin - self.nout) == 1

        def isBalanced(self):
            return self.nin == self.nout

        def __hash__(self):
            return hash(self.km1mer)

        def __str__(self):
            return self.km1mer

    def __init__(self, strIter, k, circularize=False):
        self.G = {}
        self.nodes = {}
        for st in strIter:
            if circularize:
                st += st[:k-1]
            for kmer, km1L, km1R in self.chop(st, k):
                nodeL, nodeR = None, None
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = self.Node(km1L)
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = self.Node(km1R)
                nodeL.nout += 1
                nodeR.nin += 1
                self.G.setdefault(nodeL, []).append(nodeR)
        self.nsemi, self.nbal, self.nneither = 0, 0, 0
        self.head, self.tail = None, None
        for node in self.nodes.itervalues():
            if node.isBalanced():
                self.nbal += 1
            elif node.isSemiBalanced():
                if node.nin == node.nout + 1:
                    self.tail = node
                if node.nin == node.nout - 1:
                    self.head = node
                self.nsemi += 1
            else:
                self.nneither += 1

    def nnodes(self):
        return len(self.nodes)

    def nedges(self):
        return len(self.G)

    def hasEulerianWalk(self):
        return self.nneither == 0 and self.nsemi == 2

    def hasEulerianCycle(self):
        return self.nneither == 0 and self.nsemi == 0

    def isEulerian(self):
        return self.hasEulerianWalk() or self.hasEulerianCycle()

    def eulerianWalkOrCycle(self):
        assert self.isEulerian()
        g = self.G
        if self.hasEulerianWalk():
            g = g.copy()
            assert self.head is not None
            assert self.tail is not None
            g.setdefault(self.tail, []).append(self.head)
        tour = []
        src = g.iterkeys().next()
        def __visit(n):
            while len(g[n]) > 0:
                dst = g[n].pop()
                __visit(dst)
            tour.append(n)

        __visit(src)
        tour = tour[::-1][:-1]

        if self.hasEulerianWalk():
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]
        return map(str, tour)
