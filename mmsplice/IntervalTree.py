"""
Accept GenomicInterval object, with attributes chrom, start, end, strand (optional).
In case of unstranded data, interval.strand == '.'
"""

import math
import random
import warnings

class Interval(object):
    ''' Compatible with pybedtools Interval
    '''

    def __init__(self,
                 chrom=None,
                 start=None,
                 end=None,
                 strand="*",
                 name=None):
        self.chrom = chrom
        assert start <= end, "start position smaller or equal to the end position"
        self.start = start
        self.end = end
        self.name = name
        self.strand = strand

    def __repr__(self):
        return "{0}, {1}:{2}-{3},{4}".format(self.name, self.chrom, self.start, self.end, self.strand)

class IntervalTree(object):
    def __init__(self):
        self.chroms = {}

    def insert(self, interval):
        # This interval is the interval to construct the tree, e.g. gtf annotations
        if interval.chrom in self.chroms:
            self.chroms[interval.chrom] = self.chroms[interval.chrom].insert(interval)
        else:
            self.chroms[interval.chrom] = IntervalNode(interval)

    def intersect(self, interval, report_func, ignore_strand=False):
        # This interval from the query
        if interval.chrom in self.chroms:
            self.chroms[interval.chrom].intersect(interval, report_func, ignore_strand=ignore_strand)
        else:
            warnings.warn("Interval chromosome name no match", UserWarning)

    def traverse(self, func):
        for item in self.chroms.itervalues():
            item.traverse(func)

class IntervalNode(object):
    def __init__(self, interval):
        self.priority = math.ceil((-1.0 / math.log(.5)) * math.log(-1.0 / (random.uniform(0, 1) - 1)))
        self.interval = interval
        self.start = interval.start
        self.end = interval.end
        self.strand = interval.strand
        self.maxend = self.end
        self.minend = self.end
        self.left = None
        self.right = None

    def insert(self, interval):
        root = self
        if interval.start > self.start:
            # insert to right tree
            if self.right:
                self.right = self.right.insert(interval)
            else:
                self.right = IntervalNode(interval)
            # rebalance tree
            if self.priority < self.right.priority:
                root = self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left = self.left.insert(interval)
            else:
                self.left = IntervalNode(interval)
            # rebalance tree
            if self.priority < self.left.priority:
                root = self.rotateright()
        if root.right and root.left:
            root.maxend = max(root.end, root.right.maxend, root.left.maxend)
            root.minend = min(root.end, root.right.minend, root.left.minend)
        elif root.right:
            root.maxend = max(root.end, root.right.maxend)
            root.minend = min(root.end, root.right.minend)
        elif root.left:
            root.maxend = max(root.end, root.left.maxend)
            root.minend = min(root.end, root.left.minend)
        return root

    def rotateright(self):
        root = self.left
        self.left = self.left.right
        root.right = self
        if self.right and self.left:
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend)
        return root

    def rotateleft(self):
        root = self.right
        self.right = self.right.left
        root.left = self
        if self.right and self.left:
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend)
        return root

    def intersect(self, interval, report_func, ignore_strand=False):
        if interval.strand == '*' or ignore_strand:  # unstranded data, not going to compare strand
            if interval.start <= self.end and interval.end >= self.start:
                report_func(self)
        else:
            if interval.start <= self.end and interval.end >= self.start and interval.strand == self.strand:
                report_func(self)
        if self.left and interval.start <= self.left.maxend:
            self.left.intersect(interval, report_func, ignore_strand=ignore_strand)
        if self.right and interval.end >= self.start:
            self.right.intersect(interval, report_func, ignore_strand=ignore_strand)

    def traverse(self, func):
        if self.left:
            self.left.traverse(func)
        func(self)
        if self.right:
            self.right.traverse(func)
