
# ---- KDTree taken from python cookbook (Eppstein) -------------------------

def dist2(p,q):
    """Squared distance between p and q."""
    d = 0
    for i in range(len(p)):
        d += (p[i]-q[i])**2
    return d

def dist(p,q):
    """Distance between p and q."""
    d = 0
    for pp, qq in zip( p, q ):
        d += ( pp - qq )**2
    return sqrt( d )

class kdtree:
    def __init__(self,dim=2,index=0):
        self.dim = dim
        self.index = index
        self.split = None

    def addPoint(self,p):
        """Include another point in the kD-tree."""
        if self.split is None:
            self.split = p
            self.left = kdtree(self.dim, (self.index + 1) % self.dim)
            self.right = kdtree(self.dim, (self.index + 1) % self.dim)
        elif self.split[self.index] < p[self.index]:
            self.left.addPoint(p)
        else:
            self.right.addPoint(p)

    def nearestNeighbor(self,q,maxdist2):
        """Find pair (d,p) where p is nearest neighbor and d is squared
        distance to p. Returned distance must be within maxdist2; if
        not, no point itself is returned.
        """
        solution = (maxdist2+1,None)
        if self.split is not None:
            solution = min(solution, (dist2(self.split,q),self.split))
            d2split = (self.split[self.index] - q[self.index])**2
            if self.split[self.index] < q[self.index]:
                solution = min(solution,
                    self.left.nearestNeighbor(q,solution[0]))
                if d2split < solution[0]:
                    solution = min(solution,
                        self.right.nearestNeighbor(q,solution[0]))
            else:
                solution = min(solution,
                    self.right.nearestNeighbor(q,solution[0]))
                if d2split < solution[0]:
                    solution = min(solution,
                        self.left.nearestNeighbor(q,solution[0]))
        return solution

# ---------------------------------------------------------------------------
