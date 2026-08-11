from sortedcontainers import SortedDict

class AVLTree(SortedDict):
    """Minimal bintrees.AVLTree replacement backed by sortedcontainers.SortedDict."""

    def floor_item(self, key):
        idx = self.bisect_right(key) - 1
        if idx < 0:
            raise KeyError(key)
        k = self.keys()[idx]
        return k, self[k]
