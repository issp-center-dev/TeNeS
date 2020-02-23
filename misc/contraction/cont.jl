using DocOpt

mutable struct Tensor
    name :: String
    bonds :: Vector{String}
end

mutable struct Tensors
    corners :: Vector{Tensor}
    edges :: Vector{Tensor}
    centers :: Vector{Tensor}
    Ccenters :: Vector{Tensor}
    ops :: Vector{Tensor}
    nrow :: Int
    ncol :: Int
    function Tensors(nrow::Integer, ncol::Integer; orig::Integer=0) 
        corners = gen_corners()
        edges = gen_edges(nrow, ncol, orig=orig)
        centers = gen_centers(nrow, ncol, orig=orig)
        Ccenters = gen_centers(nrow, ncol, conj=true, orig=orig)
        ops = gen_ops(nrow, ncol, orig=orig)
        tns = new(corners, edges, centers, Ccenters, ops, nrow, ncol)
        contract!(tns)
        check_bonds(tns)
        return tns
    end
end

function gen_corners()
    corners = Tensor[]
    push!(corners, Tensor("C_tl", ["C_tl_b", "C_tl_r"]))
    push!(corners, Tensor("C_tr", ["C_tr_l", "C_tr_b"]))
    push!(corners, Tensor("C_br", ["C_br_t", "C_br_l"]))
    push!(corners, Tensor("C_bl", ["C_bl_r", "C_bl_t"]))
    return corners
end

function gen_edges(nrow, ncol; orig=0)
    edges = Tensor[]
    for icol in 1:ncol
        i = orig+icol-1
        name = "eT_t$(i)"
        push!(edges, Tensor(name, ["$(name)_$dir" for dir in ("l", "r", "b", "Cb")]))
    end
    for irow in 1:nrow
        i = orig+irow-1
        name = "eT_r$(i)"
        push!(edges, Tensor(name, ["$(name)_$dir" for dir in ("t", "b", "l", "Cl")]))
    end
    for icol in reverse(1:ncol)
        i = orig+icol-1
        name = "eT_b$(i)"
        push!(edges, Tensor(name, ["$(name)_$dir" for dir in ("r", "l", "t", "Ct")]))
    end
    for irow in reverse(1:nrow)
        i = orig+irow-1
        name = "eT_l$(i)"
        push!(edges, Tensor(name, ["$(name)_$dir" for dir in ("b", "t", "r", "Cr")]))
    end
    return edges
end

function gen_centers(nrow, ncol; conj::Bool=false, orig::Integer=0)
    centers = Tensor[]
    for irow in 1:nrow
        irow = orig+irow-1
        for icol in 1:ncol
            icol = orig+icol-1
            bond_base = "Tn_$(irow)_$(icol)"
            if conj
                name = "conj(Tn_$(irow)_$(icol))"
                bond_prefix = "C"
            else
                name = "Tn_$(irow)_$(icol)"
                bond_prefix = ""
            end
            bonds =  ["$(bond_prefix)$(bond_base)_$dir" for dir in ("l", "t", "r", "b")]
            push!(bonds, "$(bond_prefix)p_$(irow)_$(icol)")
            push!(centers, Tensor(name, bonds))
        end
    end
    return centers
end

function gen_ops(nrow, ncol; orig=0)
    ops = Tensor[]
    for irow in 1:nrow
        irow = orig + irow -1
        for icol in 1:ncol
            icol = orig + icol -1
            push!(ops, Tensor("op_$(irow)_$(icol)", ["p_$(irow)_$icol", "Cp_$(irow)_$icol"] ))
        end
    end
    return ops
end

function contract!(tensors::Tensors)
    nrow = tensors.nrow
    ncol = tensors.ncol
    corners = tensors.corners
    edges = tensors.edges
    centers = tensors.centers
    Ccenters = tensors.Ccenters
    ops = tensors.ops

    function coord2index(irow, icol)
        return icol + ncol*(irow-1)
    end

    # contract corners
    corners[1].bonds[1] = edges[end].bonds[2]
    corners[1].bonds[2] = edges[1].bonds[1]

    corners[2].bonds[1] = edges[ncol  ].bonds[2]
    corners[2].bonds[2] = edges[ncol+1].bonds[1]

    corners[3].bonds[1] = edges[ncol+nrow  ].bonds[2]
    corners[3].bonds[2] = edges[ncol+nrow+1].bonds[1]

    corners[4].bonds[1] = edges[2ncol+nrow  ].bonds[2]
    corners[4].bonds[2] = edges[2ncol+nrow+1].bonds[1]

    # contract top/bottom edges
    for icol in 1:ncol
        ie = icol
        ie2 = icol + ncol+nrow
        if icol < ncol
            edges[ie].bonds[2] = edges[ie+1].bonds[1]
            edges[ie2].bonds[2] = edges[ie2+1].bonds[1]
        end
        edges[ie].bonds[3] =  centers[ie].bonds[2]
        edges[ie].bonds[4] = Ccenters[ie].bonds[2]
        edges[ie2].bonds[3] =  centers[coord2index(nrow, ncol-icol+1)].bonds[4]
        edges[ie2].bonds[4] = Ccenters[coord2index(nrow, ncol-icol+1)].bonds[4]
    end

    # contract right/left edges
    for irow in 1:nrow
        ie = irow+ncol
        ie2 = ie + ncol+nrow
        if irow < nrow
            edges[ie].bonds[2] = edges[ie+1].bonds[1]
            edges[ie2].bonds[2] = edges[ie2+1].bonds[1]
        end
        edges[ie].bonds[3] =  centers[coord2index(irow, ncol)].bonds[3]
        edges[ie].bonds[4] = Ccenters[coord2index(irow, ncol)].bonds[3]
        edges[ie2].bonds[3] =  centers[coord2index(nrow-irow+1, 1)].bonds[1]
        edges[ie2].bonds[4] = Ccenters[coord2index(nrow-irow+1, 1)].bonds[1]
    end

    # contract centers
    for icol in 1:ncol
        for irow in 1:nrow
            index = coord2index(irow, icol)
            if icol > 1
                centers[index].bonds[1] = centers[coord2index(irow, icol-1)].bonds[3]
            end
            if irow > 1
                centers[index].bonds[2] = centers[coord2index(irow-1, icol)].bonds[4]
            end
            if icol < ncol
                centers[index].bonds[3] = centers[coord2index(irow, icol+1)].bonds[1]
            end
            if irow < nrow
                centers[index].bonds[4] = centers[coord2index(irow+1, icol)].bonds[2]
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

function ctm_bonds(tensors::Tensors)
    ret = Set{String}()
    for e in tensors.edges
        push!(ret, e.bonds[1])
        push!(ret, e.bonds[2])
    end
    return ret
end

function check_bonds(tensors::Tensors)
    fails = Set{String}()
    counts = Dict{String, Int}()
    for T in Iterators.flatten((tensors.corners, tensors.edges, tensors.centers, tensors.Ccenters, tensors.ops))
        for b in T.bonds
            c = get(counts, b, 0)
            if c>1
                push!(fails, b)
            end
            counts[b] = c+1
        end
    end
    if length(fails) > 0
        println("FAILED: following bonds have more than 2 endpoints")
        for b in fails
            println(b)
        end
        error()
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

function dump(io::IO, tensors::Tensors; pdim::Integer=2, vdim::Integer=4, Cdim::Integer=vdim^2, lang::AbstractString="cpp")
    for t in tensors.corners
        println(io, t)
    end
    println(io)

    for t in tensors.edges
        println(io, t)
    end
    println(io)

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

    print(io, "bond_dim ", Cdim)
    for b in ctm_bonds(tensors)
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

function dump_fn_declare(io::IO, tensors::Tensors; lang::AbstractString="cpp")

    iter = Iterators.flatten((tensors.corners, tensors.edges, tensors.centers, tensors.ops))

    if lang=="cpp"
        println(io, "template <class tensor>")
        println(io, "typename tensor::value_type")
        println(io, "Contract_$(tensors.nrow)x$(tensors.ncol)(")
        args = String[]
        for T in iter
            push!(args, "  const tensor &$(T.name)")
        end
        println(io, join(args, ",\n"))
        println(io, ")")
    else
        println(io, "def Contract_$(tensors.nrow)x$(tensors.ncol)(")
        args = String[]
        for T in iter
            push!(args, "    $(T.name): np.ndarray")
        end
        println(io, join(args, ",\n"))
        println(io, ") -> np.ndarray:")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    doc = """

    Usage:
        cont.jl [--output=<output>] [--pdim=<pdim>] [--vdim=<vdim>] [--cdim=<cdim>] [--lang=<lang>] [--tdt_path=<tdt_path>] <nrow> <ncol>

    Options:
      --output=<output>         filename to be saved [default: STDOUT]
      --pdim=<pdim>             dim. of physical bonds [default: 2]
      --vdim=<vdim>             dim. of virtual bonds [default: 4]
      --cdim=<cdim>             dim. of ctm bonds [default: vdim^2]
      --lang=<lang>             cpp or python [default: cpp]
      --tdt_path=<tdt_path>     path to tdt.py (if empty, don't run tdt.py) [default: ""]

    """

    ARGS_ORIG = ARGS[:]
    args = docopt(doc, version=v"1.0.0")
    output = args["--output"]
    io = output=="STDOUT" ? stdout : open(output, "w")
    pdim = parse(Int, args["--pdim"])
    vdim = parse(Int, args["--vdim"])
    cdim = args["--cdim"] == "vdim^2" ? vdim^2 : parse(Int, args["--cdim"])
    lang = args["--lang"]
    tdt_path = args["--tdt_path"]

    tensors = Tensors(parse(Int, args["<nrow>"]), parse(Int, args["<ncol>"]), orig=0)

    println(io, """
            # This file is generated by $(PROGRAM_FILE) $(join(ARGS_ORIG, " "))
            """)

    dump(io, tensors, pdim=pdim, vdim=vdim, Cdim=cdim, lang=lang)
    println(io)
    println(io)
    dump_fn_declare(io, tensors, lang=lang)

    if io != stdout && length(tdt_path) > 0
        close(io)
        res = read(`python $(tdt_path) $(output)`, String)
        open(output, "a") do io
            if lang == "cpp"
                println(io, "{")
                indent = "  "
                write_return = false
            else
                indent = "    "
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
