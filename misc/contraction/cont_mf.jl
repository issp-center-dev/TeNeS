# /* TeNeS - Massively parallel tensor network solver /
# / Copyright (C) 2019- The University of Tokyo */
#
# /* This program is free software: you can redistribute it and/or modify /
# / it under the terms of the GNU General Public License as published by /
# / the Free Software Foundation, either version 3 of the License, or /
# / (at your option) any later version. */
#
# /* This program is distributed in the hope that it will be useful, /
# / but WITHOUT ANY WARRANTY; without even the implied warranty of /
# / MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
# / GNU General Public License for more details. */
#
# /* You should have received a copy of the GNU General Public License /
# / along with this program. If not, see http://www.gnu.org/licenses/. */

using DocOpt

mutable struct Tensor
    name :: String
    bonds :: Vector{String}
end

mutable struct Tensors
    centers :: Vector{Tensor}
    Ccenters :: Vector{Tensor}
    ops :: Vector{Tensor}
    nrow :: Int
    ncol :: Int
    function Tensors(nrow::Integer, ncol::Integer; orig::Integer=0, pass_as_vector::Bool=false)
        centers = gen_centers(nrow, ncol; orig, pass_as_vector)
        Ccenters = gen_centers(nrow, ncol; conj=true, orig, pass_as_vector)
        ops = gen_ops(nrow, ncol, orig=orig, pass_as_vector=pass_as_vector)
        tns = new(centers, Ccenters, ops, nrow, ncol)
        contract!(tns)
        check_bonds(tns)
        return tns
    end
end

function gen_centers(nrow, ncol; conj::Bool=false, orig::Integer=0, pass_as_vector::Bool=false)
    centers = Tensor[]
    for irow in 1:nrow
        irow = orig+irow-1
        for icol in 1:ncol
            icol = orig+icol-1
            bond_base = "Tn_$(irow)_$(icol)"
            if conj
                name = pass_as_vector ? "conj(*(Tn[$(irow)][$(icol)]))" : "conj(Tn_$(irow)_$(icol))"
                bond_prefix = "C"
            else
                name = pass_as_vector ? "*(Tn[$(irow)][$(icol)])" : "Tn_$(irow)_$(icol)"
                bond_prefix = ""
            end
            bonds =  ["$(bond_prefix)$(bond_base)_$dir" for dir in ("l", "t", "r", "b")]
            push!(bonds, "$(bond_prefix)p_$(irow)_$(icol)")
            push!(centers, Tensor(name, bonds))
        end
    end
    return centers
end

function gen_ops(nrow, ncol; orig=0, pass_as_vector::Bool=false)
    ops = Tensor[]
    for irow in 1:nrow
        irow = orig + irow -1
        for icol in 1:ncol
            icol = orig + icol -1
            name = pass_as_vector ? "*(op[$(irow)][$(icol)])" : "op_$(irow)_$(icol)"
            push!(ops, Tensor(name, ["p_$(irow)_$icol", "Cp_$(irow)_$icol"] ))
        end
    end
    return ops
end

function contract!(tensors::Tensors)
    nrow = tensors.nrow
    ncol = tensors.ncol
    centers = tensors.centers
    Ccenters = tensors.Ccenters
    ops = tensors.ops

    function coord2index(irow, icol)
        return icol + ncol*(irow-1)
    end

    # contract centers
    for icol in 1:ncol
        for irow in 1:nrow
            index = coord2index(irow, icol)
            if icol == 1
                Ccenters[index].bonds[1] = centers[index].bonds[1]
            else
                centers[index].bonds[1] = centers[coord2index(irow, icol-1)].bonds[3]
                Ccenters[index].bonds[1] = Ccenters[coord2index(irow, icol-1)].bonds[3]
            end
            if icol == ncol
                Ccenters[index].bonds[3] = centers[index].bonds[3]
            end

            if irow == 1
                Ccenters[index].bonds[2] = centers[index].bonds[2]
            else
                centers[index].bonds[2] = centers[coord2index(irow-1, icol)].bonds[4]
                Ccenters[index].bonds[2] = Ccenters[coord2index(irow-1, icol)].bonds[4]
            end
            if irow == nrow
                Ccenters[index].bonds[4] = centers[index].bonds[4]
            end
        end
    end
end

function phys_bonds(tensors::Tensors)
    ret = Set{String}()
    for op in tensors.ops
        push!(ret, op.bonds[1])
        push!(ret, op.bonds[2])
    end
    return ret
end

function virtual_bonds(tensors::Tensors)
    ret = Set{String}()
    for c in tensors.centers
        for b in c.bonds
            push!(ret, b)
        end
    end
    for c in tensors.Ccenters
        for b in c.bonds
            push!(ret, b)
        end
    end
    return ret
end

function check_bonds(tensors::Tensors)
    counts = Dict{String, Int}()
    for T in Iterators.flatten((tensors.centers, tensors.Ccenters, tensors.ops))
        for b in T.bonds
            c = get(counts, b, 0)
            if c>1
                push!(fails, b)
            end
            counts[b] = c+1
        end
    end
    OK = true
    for (b,c) in counts
        if c != 2
            println("FAILED: bond $(b) appears $(c) time(s)")
            OK = false
        end
    end
    if !OK
        dump_tensors(stdout, tensors)
        error("some bonds have trouble")
    end
end


import Base: print
function print(io::IO, tensor::Tensor)
    print(io, "tensor ", tensor.name, " ")
    for b in tensor.bonds
        print(io, " ", b)
    end
end

function dump(tensors::Tensors; args...)
    dump(stdout, tensors; args...)
end

function dump_tensors(io::IO, tensors::Tensors)
    for t in tensors.centers
        println(io, t)
    end
    println(io)

    for t in tensors.Ccenters
        println(io, t)
    end
    println(io)

    for t in tensors.ops
        println(io, t)
    end
    println(io)
end

function dump(io::IO, tensors::Tensors; pdim::Integer=2, vdim::Integer=4, lang::AbstractString="cpp")
    dump_tensors(io, tensors)

    print(io, "bond_dim ", pdim)
    for b in phys_bonds(tensors)
        print(io, " ", b)
    end
    println(io)

    print(io, "bond_dim ", vdim)
    for b in virtual_bonds(tensors)
        print(io, " ", b)
    end
    println(io)
    println(io)

    if lang == "cpp"
        indent = 2
        style = "mptensor"
    else
        indent = 4
        style = "numpy"
    end

    println(io, "style ", style)
    println(io, "indent ", indent)
end

function dump_fn_declare(io::IO, tensors::Tensors; lang::AbstractString="cpp", pass_as_vector::Bool=false)

    iter = Iterators.flatten((tensors.centers, tensors.ops))

    if lang=="cpp"
        println(io, "template <class tensor>")
        println(io, "typename tensor::value_type")
        println(io, "Contract_MF_$(tensors.nrow)x$(tensors.ncol)(")
        if pass_as_vector
            println(io, "  const std::vector<std::vector<const tensor*>> &Tn,")
            println(io, "  const std::vector<std::vector<const tensor*>> &op")
        else
            args = String[]
            for T in iter
                push!(args, "  const tensor &$(T.name)")
            end
            println(io, join(args, ",\n"))
        end
        println(io, ")")
    else
        println(io, "def Contract_MF_$(tensors.nrow)x$(tensors.ncol)(")
        if pass_as_vector
            println(io, "    Tn: List[List[np.ndarray]],")
            println(io, "    op: List[List[np.ndarray]]")
        else
            args = String[]
            for T in iter
                push!(args, "    $(T.name): np.ndarray")
            end
            println(io, join(args, ",\n"))
        end
        println(io, ") -> np.ndarray:")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    doc = """

    Usage:
        cont.jl [--output=<output>] [--pdim=<pdim>] [--vdim=<vdim>] [--lang=<lang>] [--tdt_path=<tdt_path>] [--pass_as_vector] <nrow> <ncol>

    Options:
        --output=<output>         filename to be saved [default: STDOUT]
        --pdim=<pdim>             dim. of physical bonds [default: 2]
        --vdim=<vdim>             dim. of virtual bonds [default: 4]
        --lang=<lang>             cpp or python [default: cpp]
        --tdt_path=<tdt_path>     path to tdt.py (if empty, don't run tdt.py) [default: ""]
        --pass_as_vector          generate API as a function of vectors
    """

    ARGS_ORIG = ARGS[:]
    args = docopt(doc, version=v"1.0.0")
    output = args["--output"]
    io = output=="STDOUT" ? stdout : open(output, "w")
    pdim = parse(Int, args["--pdim"])
    vdim = parse(Int, args["--vdim"])
    lang = args["--lang"]
    tdt_path = args["--tdt_path"]
    pass_as_vector = args["--pass_as_vector"]

    tensors = Tensors(parse(Int, args["<nrow>"]), parse(Int, args["<ncol>"]), orig=0, pass_as_vector=pass_as_vector)

    println(io, """
            # This file is generated by $(PROGRAM_FILE) $(join(ARGS_ORIG, " "))
            """)

    dump(io, tensors, pdim=pdim, vdim=vdim, lang=lang)
    println(io)
    println(io)
    dump_fn_declare(io, tensors, lang=lang, pass_as_vector=pass_as_vector)

    if io != stdout && length(tdt_path) > 0
        close(io)
        res = read(`python $(tdt_path) $(output)`, String)
        open(output, "a") do io
            if lang == "cpp"
                println(io, "{")
                indent = "  "
                write_return = false

                if pass_as_vector
                    println(io, "#ifndef NDEBUG")
                    println(io, "  const size_t nrow = Tn.size();")
                    println(io, "  const size_t ncol = Tn[0].size();")
                    println(io, "  for(const auto& r: Tn) assert(r.size() == ncol);")
                    println(io, "  assert(op.size() == nrow);")
                    println(io, "  for(const auto& r: op) assert(r.size() == ncol);")
                    println(io, "#endif")
                end
            else
                indent = "    "

                if pass_as_vector
                    println(io, "    nrow = len(Tn)")
                    println(io, "    ncol = len(Tn[0])")
                    println(io, "    for r in Tn:")
                    println(io, "        assert(len(r) == ncol)")
                    println(io, "    assert(len(op) == nrow")
                    println(io, "    for r in op:")
                    println(io, "        assert(len(r) == ncol)")
                end
            end
            for line in split(res, "\n")
                if length(line) == 0
                    continue
                end
                if lang == "cpp" && !write_return
                    if ! startswith(line, "//")
                        write_return = true
                        println(io, indent, "return")
                    end
                end
                println(io, indent, line)
            end
            if lang == "cpp"
                println(io, indent, ";")
                println(io, "}")
            end
        end
    end
end
