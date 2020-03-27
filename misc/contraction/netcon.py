"""Get optimal contraction sequence using netcon algorithm

Reference:
    R. N. C. Pfeifer, et al.: Phys. Rev. E 90, 033315 (2014)
"""

# tdt.py, netcon.py and config.py are originally from https://github.com/smorita/Tensordot
# They are redistributed under the MIT license (see LICENSE in this directory)

__author__ = "Satoshi MORITA <morita@issp.u-tokyo.ac.jp>"
__date__ = "24 March 2016"

import sys
import logging
import time
import config
import itertools


class TensorFrame:
    """Tensor class for netcon.

    Attributes:
        rpn: contraction sequence with reverse polish notation.
        bits: bits representation of contracted tensors.
        bonds: list of uncontracted bonds.
        is_new: a flag.
    """

    def __init__(self,rpn=[],bits=0,bonds=[],cost=0.0,is_new=True):
        self.rpn = rpn[:]
        self.bits = bits
        self.bonds = bonds
        self.cost = cost
        self.is_new = is_new

    def __repr__(self):
        return "TensorFrame({0}, bonds={1}, cost={2:.6e}, bits={3}, is_new={4})".format(
            self.rpn, self.bonds, self.cost, self.bits, self.is_new)

    def __str__(self):
        return "{0} : bonds={1} cost={2:.6e} bits={3} new={4}".format(
            self.rpn, self.bonds, self.cost, self.bits, self.is_new)


class NetconOptimizer:
    def __init__(self, prime_tensors, bond_dims):
        self.prime_tensors = prime_tensors
        self.BOND_DIMS = bond_dims[:]

    def optimize(self):
        """Find optimal contraction sequence.

        Args:
            tn: TensorNetwork in tdt.py
            bond_dims: List of bond dimensions.

        Return:
            rpn: Optimal contraction sequence with reverse polish notation.
            cost: Total contraction cost.
        """
        tensordict_of_size = self.init_tensordict_of_size()

        n = len(self.prime_tensors)
        xi_min = float(min(self.BOND_DIMS))
        mu_cap = 1.0
        prev_mu_cap = 0.0 #>=0

        while len(tensordict_of_size[-1])<1:
            logging.info("netcon: searching with mu_cap={0:.6e}".format(mu_cap))
            next_mu_cap = sys.float_info.max
            for c in range(2,n+1):
                for d1 in range(1,c//2+1):
                    d2 = c-d1
                    t1_t2_iterator = itertools.combinations(tensordict_of_size[d1].values(), 2) if d1==d2 else itertools.product(tensordict_of_size[d1].values(), tensordict_of_size[d2].values())
                    for t1, t2 in t1_t2_iterator:
                        if self.are_overlap(t1,t2): continue
                        if self.are_direct_product(t1,t2): continue

                        cost = self.get_contracting_cost(t1,t2)
                        bits = t1.bits ^ t2.bits

                        if next_mu_cap <= cost:
                            pass
                        elif mu_cap < cost:
                            next_mu_cap = cost
                        elif t1.is_new or t2.is_new or prev_mu_cap < cost:
                            t_old = tensordict_of_size[c].get(bits)
                            if t_old is None or cost < t_old.cost:
                                tensordict_of_size[c][bits] = self.contract(t1,t2)
            prev_mu_cap = mu_cap
            mu_cap = max(next_mu_cap, mu_cap*xi_min)
            for s in tensordict_of_size:
                for t in s.values(): t.is_new = False

            logging.debug("netcon: tensor_num=" +  str([ len(s) for s in tensordict_of_size]))

        t_final = tensordict_of_size[-1][(1<<n)-1]
        return t_final.rpn, t_final.cost


    def init_tensordict_of_size(self):
        """tensordict_of_size[k][bits] == calculated lowest-cost tensor which is contraction of k+1 prime tensors and whose bits == bits"""
        tensordict_of_size = [{} for size in range(len(self.prime_tensors)+1)]
        for t in self.prime_tensors:
            rpn = t.name
            bits = 0
            for i in rpn:
                if i>=0: bits += (1<<i)
            bonds = frozenset(t.bonds)
            cost = 0.0
            tensordict_of_size[1].update({bits:TensorFrame(rpn,bits,bonds,cost)})
        return tensordict_of_size


    def get_contracting_cost(self,t1,t2):
        """Get the cost of contraction of two tensors."""
        cost = 1.0
        for b in (t1.bonds | t2.bonds):
            cost *= self.BOND_DIMS[b]
        cost += t1.cost + t2.cost
        return cost


    def contract(self,t1,t2):
        """Return a contracted tensor"""
        assert (not self.are_direct_product(t1,t2))
        rpn = t1.rpn + t2.rpn + [-1]
        bits = t1.bits ^ t2.bits # XOR
        bonds = frozenset(t1.bonds ^ t2.bonds)
        cost = self.get_contracting_cost(t1,t2)
        return TensorFrame(rpn,bits,bonds,cost)


    def are_direct_product(self,t1,t2):
        """Check if two tensors are disjoint."""
        return (t1.bonds).isdisjoint(t2.bonds)


    def are_overlap(self,t1,t2):
        """Check if two tensors have the same basic tensor."""
        return (t1.bits & t2.bits)>0


    def print_tset(self,tensors_of_size):
        """Print tensors_of_size. (for debug)"""
        for level in range(len(tensors_of_size)):
            for i,t in enumerate(tensors_of_size[level]):
                print(level,i,t)

