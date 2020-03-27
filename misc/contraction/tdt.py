#!/usr/bin/env python

# tdt.py, netcon.py and config.py are originally from https://github.com/smorita/Tensordot
# They are redistributed under the MIT license (see LICENSE in this directory)

import sys
import argparse
import logging
import time
import config
import netcon

# Global variables
TENSOR_NAMES = []
TENSOR_MATH_NAMES = []
BOND_NAMES = []
BOND_DIMS = []
VECTORS = []
FINAL_ORDER = None

class Tensor:
    def __init__(self,name=None,bonds=[]):
        if name==None:
            self.name = []
        elif isinstance(name, list):
            self.name = name[:]
        else:
            self.name = [name]
        self.bonds = bonds[:]

    def __repr__(self):
        return "Tensor(" + str(self.name) + ", " + str(self.bonds) +")"

    def __str__(self):
        return str(self.name) + ", " + str(self.bonds)


class Bond:
    def __init__(self,t0=-1,t1=-1):
        self.t0 = t0
        self.t1 = t1

    def __str__(self):
        return "({0},{1})".format(self.t0,self.t1)

    def isFree(self):
        return (self.t0 < 0 or self.t1 < 0)

    def connect(self,tensor_index):
        assert self.isFree(), "edge already connected to two tensors"
        if self.t0<0:
            self.t0 = tensor_index
        else:
            assert not self.t0==tensor_index, "edge connects to the same tensor"
            self.t1 = tensor_index

    def has(self,tensor_index):
        return (self.t0==tensor_index or self.t1==tensor_index)


class TensorNetwork:
    def __init__(self):
        self.tensors = []
        self.bonds = []
        self.total_memory = 0.0
        self.max_memory = 0.0
        self.cpu_cost = 0.0

    def __str__(self):
        s = ""
        for i,t in enumerate(self.tensors):
            s += "tensor {0} : {1}\n".format(i,t)
        for i,b in enumerate(self.bonds):
            s += "bond {0} : {1}, {2} {3}\n".format(i,BOND_NAMES[i],b,BOND_DIMS[i])
        s += "memory : {0}\n".format(self.total_memory)
        s += "cpu : {0}\n".format(self.cpu_cost)
        return s


    def clone(self):
        tn = TensorNetwork()
        tn.total_memory = self.total_memory
        tn.max_memory = self.max_memory
        tn.cpu_cost = self.cpu_cost
        tn.bonds = [ Bond(b.t0,b.t1) for b in self.bonds ]
        tn.tensors = [ Tensor(t.name,t.bonds) for t in self.tensors ]
        return tn


    def output_log(self,prefix=""):
        if not prefix=="": prefix += " "
        for i,t in enumerate(self.tensors):
            logging.info(prefix + "tensor{0} : {1} {2}".format(i,TENSOR_NAMES[i],t.bonds))
        for i,b in enumerate(self.bonds):
            logging.info(prefix + "bond{0} : {1} {2} {3}".format(i,BOND_NAMES[i],b,BOND_DIMS[i]))


    def add_tensor(self, t_name, b_names):
        t_index = len(self.tensors)
        b_indexs = []
        for b in b_names:
            if b not in BOND_NAMES:
                self.bonds.append(Bond())
                BOND_NAMES.append(b)
                BOND_DIMS.append(config.DEFAULT_BOND_DIM)

            i = BOND_NAMES.index(b)
            self.bonds[i].connect(t_index)
            b_indexs.append(i)

        TENSOR_NAMES.append(t_name)
        self.tensors.append(Tensor(t_index,b_indexs))


    def find_bonds(self, tensor_a, tensor_b):
        bonds_a = self.tensors[tensor_a].bonds
        bonds_b = self.tensors[tensor_b].bonds
        contract = [ b for b in bonds_a if b in bonds_b]
        replaced_a = [ b for b in bonds_a if b not in bonds_b ]
        replaced_b = [ b for b in bonds_b if b not in bonds_a ]
        return contract, replaced_a, replaced_b


    def contract(self, t0, t1, bc, br0, br1):
        tn = self.clone()

        # create the contracted tensor
        t_new = tn.tensors[t0]
        ## change names of tensors using Reverse Polish Notation
        t_new.name = self.tensors[t0].name+self.tensors[t1].name+[-1]
        ## remove contracted bonds
        for b in bc: t_new.bonds.remove(b)
        ## add bonds from deleted tensor
        for b in br1: t_new.bonds.append(b)

        # clear the removed tensor
        tn.tensors[t1] = Tensor()

        # update bonds
        bonds = tn.bonds
        ## remove contracted bonds from the bond list
        for b in bc: bonds[b].t0 = bonds[b].t1 = -1
        ## change bond connections
        old_idx = t1
        new_idx = t0
        for b in br1:
            if bonds[b].t0==old_idx: bonds[b].t0=new_idx
            elif bonds[b].t1==old_idx: bonds[b].t1=new_idx

        return tn


def get_memory(tn_orig,rpn):
    """Caluculate memory cost for contractions from Reverse Polish Notation"""
    tn = tn_orig.clone()
    cost = []
    for item in rpn:
        if item==-1:
            c1 = cost.pop()
            c0 = cost.pop()

            index1 = c1[0]
            index0 = c0[0]

            t0 = tn.tensors[index0]
            t1 = tn.tensors[index1]

            bc, br0, br1 = tn.find_bonds(index0, index1)

            mem_start = c0[2] + c1[2]
            mem_end = 1.0
            for b in br0 + br1: mem_end *= BOND_DIMS[b]
            mem_req = max(c0[1]+c1[2], c0[1]+c1[3], c0[2]+c1[1], c0[3]+c1[1], mem_end+c0[3]+c1[3])

            tn = tn.contract(index0, index1, bc, br0, br1)

            cost.append( (index0, mem_req, mem_start, mem_end) )

        else:
            t = tn.tensors[item]
            val = 1.0
            for b in t.bonds: val *= BOND_DIMS[b]
            cost.append( (item, val, val, val) ) # (index, mem_req, mem_start, mem_end)

    return cost[0][1]


def get_math(rpn):
    """Generate mathematical formula from Reverse Polish Notation"""
    stack = []
    for c in rpn:
        if c==-1:
            t1 = stack.pop()
            t0 = stack.pop()
            new_name = "("+t0+"*"+t1+")"
            stack.append( new_name )

        else:
            stack.append(TENSOR_MATH_NAMES[c])
    return stack[0]


def get_script(tn_orig,rpn):
    """Generate tensordot script from Reverse Polish Notation"""
    tn = tn_orig.clone()
    index = []
    name = []
    for c in rpn:
        if c==-1:
            index1 = index.pop()
            index0 = index.pop()
            name1 = name.pop()
            name0 = name.pop()

            t0 = tn.tensors[index0]
            t1 = tn.tensors[index1]

            bc, br0, br1 = tn.find_bonds(index0, index1)

            axes0 = [ t0.bonds.index(b) for b in bc]
            axes1 = [ t1.bonds.index(b) for b in bc]

            tn = tn.contract(index0, index1, bc, br0, br1)

            trace = (len(br0)==0 and len(br1)==0)
            new_name = tensordot_script(name0,name1,axes0,axes1,trace)

            index.append(index0)
            name.append(new_name)

        else:
            index.append(c)
            name.append([TENSOR_NAMES[c]])

    bond_order = tn.tensors[index.pop()].bonds

    return name.pop(), bond_order


def tensordot_script(name0,name1,axes0,axes1,trace=False):
    if config.STYLE == "numpy":
        func_name = config.NUMPY+".tensordot"
        axes = "(" + str(axes0) + ", " + str(axes1) + ")"
    elif config.STYLE == "mptensor":
        func_name = "trace" if trace else "tensordot"
        str_axes0 = str(tuple(axes0)) if len(axes0)>1 else "("+str(axes0[0])+")"
        str_axes1 = str(tuple(axes1)) if len(axes1)>1 else "("+str(axes1[0])+")"
        axes = "Axes" + str_axes0 + ", " + "Axes" + str_axes1

    script = []
    script.append( func_name + "(" )
    for l in name0: script.append(config.INDENT + l)
    script[-1] += ", " + name1[0]
    for i in range(1,len(name1)): script.append(config.INDENT + name1[i])
    script[-1] += ", " + axes
    script.append( ")" )

    return script


def transpose_script(name,axes):
    if config.STYLE == "numpy":
        func_name = config.NUMPY+".transpose"
        axes = str(axes)
    elif config.STYLE == "mptensor":
        func_name = "transpose"
        str_axes = str(tuple(axes)) if len(axes)>1 else "("+str(axes[0])+")"
        axes = "Axes" + str_axes

    script = []
    script.append( func_name + "(" )
    for l in name: script.append(config.INDENT + l)
    script[-1] += ", " + axes
    script.append( ")" )

    return script


def multiply_vector_script(name,vec_list,rank):
    if config.STYLE == "numpy":
        newaxis = ","+config.NUMPY+".newaxis"
        script = "("+name
        for axis,vec_name in vec_list:
            if axis==rank-1:
                script += "*"+vec_name
            else:
                script += "*"+vec_name+"[:"+newaxis*(rank-axis-1)+"]"
        script += ")"
    if config.STYLE == "mptensor":
        arg = []
        for axis,vec_name in vec_list:
            arg.append(vec_name)
            arg.append(str(axis))
        script = name + ".multiply_vector(" + ",".join(arg)+ ")"

    math = "("+name
    for axis,vec_name in vec_list: math += "*"+vec_name
    math += ")"
    return script, math


def add_transpose(tn,script,bond_order):
    if FINAL_ORDER == None: return script, bond_order

    f_order = [ BOND_NAMES.index(b) for b in FINAL_ORDER ]

    if not sorted(f_order)==sorted(bond_order):
        logging.warning("The final bond order is invalid. It is ignored.")
        return script,bond_order
    elif f_order == bond_order:
        logging.info("The final bond order was requested, but Transpose is not necessary.")
        return script, bond_order

    axes = [ bond_order.index(b) for b in f_order ]
    return transpose_script(script,axes), f_order


def add_multiply_vector(tn):
    """Change names of tensors by vector multiplications"""
    if len(VECTORS)==0: return

    mod_list = [[] for _ in TENSOR_NAMES]
    for v_name, b_name in VECTORS:
        assert (b_name in BOND_NAMES), "Vector ({0}) is multiplied to a non-existent bond.".format(v_name)
        b_index = BOND_NAMES.index(b_name)
        bond = tn.bonds[b_index]
        t0, t1 = bond.t0, bond.t1
        # find a smaller tensor
        if t0>-1 and t1>-1:
            mem0 = mem1 = 1.0
            for b in tn.tensors[t0].bonds: mem0 *= BOND_DIMS[b]
            for b in tn.tensors[t1].bonds: mem1 *= BOND_DIMS[b]
            t = t0 if mem0<mem1 else t1
        else:
            t = max(t0, t1)
        axis = tn.tensors[t].bonds.index(b_index)
        mod_list[t].append((axis,v_name))
        logging.debug("vector: "+v_name+" on bond"+str(b_index)+" -> tensor"+str(t))

    for i,l in enumerate(mod_list):
        if len(l)==0: continue
        rank = len(tn.tensors[i].bonds)
        new_name, new_math = multiply_vector_script(TENSOR_NAMES[i],sorted(l),rank)
        logging.info("vector: "+TENSOR_NAMES[i]+" -> "+new_math)
        TENSOR_NAMES[i] = new_name
        TENSOR_MATH_NAMES[i] = new_math


def read_file(infile, tn):
    """Read input file"""
    global FINAL_ORDER

    for line in infile:
        data = line.split()
        if data==[]: continue

        command = data[0].lower()
        if command=="style":
            set_style(data[1].lower())
        elif command=="numpy":
            config.NUMPY = data[1]
        elif command=="indent":
            config.INDENT = " " * int(data[1])
        elif command=="default_dimension":
            # Should be set the top of input file.
            config.DEFAULT_BOND_DIM = int(data[1])
        elif command=="debug" or command=="verbose":
            config.LOGGING_LEVEL = logging.DEBUG

        elif command=="tensor":
            tn.add_tensor(data[1], data[2:])
        elif command=="bond":
            for b in data[1:-1]: set_bond_dim(b, int(data[-1]))
        elif command=="bond_dim":
            for b in data[2:]: set_bond_dim(b, int(data[1]))
        elif command=="order":
            FINAL_ORDER = data[1:]
        elif command=="vector":
            VECTORS.append((data[1], data[2]))
    infile.close()


def set_bond_dim(bond_name, dim):
    BOND_DIMS[ BOND_NAMES.index(bond_name) ] = dim


def set_style(style):
    if style=="numpy":
        config.STYLE = "numpy"
        config.COMMENT_PREFIX = "#"
    elif style=="mptensor":
        config.STYLE = "mptensor"
        config.COMMENT_PREFIX = "//"


def check_bond_order(tn):
    return FINAL_ORDER == None or \
        frozenset(FINAL_ORDER) == frozenset( BOND_NAMES[i] for i,b in enumerate(tn.bonds) if b.isFree() )


def check_vector():
    for v in VECTORS:
        if v[1] not in BOND_NAMES: return False
    return True


def output_result(outfile,script,math_script,cpu,mem,bond_order,input_file):
    final_bonds = "(" + ", ".join([BOND_NAMES[b] for b in bond_order]) + ")"
    BR = "\n"
    SP = " "
    output = [config.COMMENT_PREFIX*30,
              config.COMMENT_PREFIX + SP + input_file,
              config.COMMENT_PREFIX*30,
              config.COMMENT_PREFIX + SP + math_script,
              config.COMMENT_PREFIX + SP + "cpu_cost= {0:g}  memory= {1:g}".format(cpu, mem),
              config.COMMENT_PREFIX + SP + "final_bond_order " + final_bonds,
              config.COMMENT_PREFIX*30]
    output += script
    outfile.write(BR.join(output) + BR)


def parse_args():
    parser = argparse.ArgumentParser(description="Code generator for tensor contruction")
    parser.add_argument('-s', metavar='style', dest='style',
                        type=str, default=None,
                        choices=['numpy', 'mptensor'],
                        help='set output style ("numpy" or "mptensor")')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                        help='verbose mode')
    parser.add_argument('-o', metavar='outfile', dest='outfile',
                        type=argparse.FileType('w'), default=sys.stdout,
                        help='write the result to outfile')
    parser.add_argument('infile',
                        type=argparse.FileType('r'),
                        help='tensor-network definition file')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    tn = TensorNetwork()

    # Read input file
    read_file(args.infile, tn)

    # Overwrite by command-line option
    set_style(args.style)
    if args.verbose:
        config.LOGGING_LEVEL = logging.DEBUG

    assert len(tn.tensors)>0, "No tensor."
    assert len(tn.bonds)>0, "No bond."
    assert check_bond_order(tn), "Final bond order is invalid."
    assert check_vector(), "Vectors will be put on non-existent bond."
    logging.basicConfig(format="%(levelname)s:%(message)s", level=config.LOGGING_LEVEL)

    tn.output_log("input")
    rpn, cpu = netcon.NetconOptimizer(tn.tensors, BOND_DIMS).optimize()
    mem = get_memory(tn, rpn)

    TENSOR_MATH_NAMES = TENSOR_NAMES[:]
    add_multiply_vector(tn)
    script, bond_order = get_script(tn, rpn)
    script, bond_order = add_transpose(tn, script, bond_order)

    output_result(args.outfile,
                  script,get_math(rpn),cpu,mem,bond_order,
                  args.infile.name)
