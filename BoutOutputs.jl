module BoutOutputs

export openboutdata
export boutcollect

using EndpointRanges
using Glob
using Lazy: @forward
using NCDatasets: NCDataset, dimnames, Attributes

"""
struct containing a set of BOUT++ dump files, opened as NCDataset objects
"""
struct BoutFiles
    data_files::Array{NCDataset,1}
    nt::Int
    nx::Int
    ny::Int
    nz::Int
    mxg::Int
    myg::Int
    is_doublenull::Bool
    ny_inner::Int
    nxpe::Int
    nype::Int
    mxsub::Int
    mysub::Int
    yproc_upper_target::Union{Nothing, Int}
end

struct BoutArray{T, N} <: AbstractArray{T, N}
    data::Array{T, N}
    attributes::Attributes
    dimensions::Tuple
end

# Define AbstractArray interface methods for BoutArray
@forward BoutArray.data (Base.size, Base.getindex, Base.setindex!,
                         Base.IndexStyle, Base.iterate, Base.length,
                         Base.similar, Base.axes)

IntOrRange = Union{Int, AbstractRange, EndpointRanges.EndpointRange}

"""
    openboutdata()

Opens a set of NetCDF BOUT++ output files for reading

Parameters
----------
datadir : String, default "."
    The directory to open files from
prefix : String, default "BOUT.dmp"
    Prefix for file names. E.g. by default BOUT++ output files are named like
    "BOUT.dmp.<i>.nc" where <i> is an integer.
"""
function openboutdata(datadir::String="."; prefix::String="BOUT.dmp")
    f = NCDataset(joinpath(datadir, prefix * ".0.nc"))

    # If BOUT_VERSION is earlier than 3.5, nz is defined differently
    bout_version = f["BOUT_VERSION"].var[:]
    if bout_version < 3.5
        error("BOUT++ version " * string(bout_version) * " is too early to be "
              * "supported by bout_outputs.jl")
    end

    nxpe = f["NXPE"].var[:]
    nype = f["NYPE"].var[:]

    npe = nxpe * nype
    if datadir[1] == '/'
        # Glob.jl doesn't support absolute paths (yet) because of Windows
        # weirdness. Passing "/" as second argument does work, not sure why!
        nfiles = length(glob(joinpath(datadir[2:end], prefix * ".*.nc"), "/"))
    else
        nfiles = length(glob(joinpath(datadir, file_pattern)))
    end
    if npe < nfiles
        @warn ("More data files found ($(nfiles)) than expected ($(npe)), "
               * "ignoring extra files")
    end

    nt = "t_array" ∈ f ? f.dim["t"] : 1
    nx = f["nx"].var[:]
    mxg = f["MXG"].var[:]
    myg = f["MYG"].var[:]
    jyseps2_1 = f["jyseps2_1"].var[:]
    jyseps1_2 = f["jyseps1_2"].var[:]
    is_doublenull = jyseps2_1 != jyseps1_2
    ny = f["ny"].var[:] + 2*myg
    if is_doublenull
        ny += 2*myg
    end
    ny_inner = f["ny_inner"].var[:]
    nz = f["nz"].var[:]
    mxsub = f["MXSUB"].var[:]
    mysub = f["MYSUB"].var[:]
    yproc_upper_target = is_doublenull ? ny_inner ÷ mysub - 1 : nothing
    if is_doublenull && ny_inner % mysub != 0
        error("trying to keep upper boundary cells but mysub=$(mysub) does "
              * "not divide ny_inner=$(ny_inner)")
    end

    file_names = ("$(prefix).$(i-1).nc" for i ∈ 2:npe)
    data_files = Array{NCDataset,1}(undef, npe)
    data_files[1] = f
    for (i, name) ∈ zip(2:nfiles, file_names)
        data_files[i] = NCDataset(joinpath(datadir, name))
    end

    return BoutFiles(data_files, nt, nx, ny, nz, mxg, myg, is_doublenull,
                     ny_inner, nxpe, nype, mxsub, mysub, yproc_upper_target)
end

"""
    boutcollect(outputs, "variable")

Read a variable from output files opened with open_bout_data()

When passing ranges to the `tind`, `xind`, `yind` or `zind` dimensions, it
might be useful to use `ibegin` and `iend` which can be used to set start and
stop values as the beginning and end of the array without requiring specific
values. `ibegin` and `iend` are re-exported by bout_post from the
EndpointRanges package.

Parameters
----------
outputs : BoutFiles
    A BoutFiles object returned by openboutdata()
varname : String
    Name of the variable to read
tind : Range, optional
    Range of indices to select in the t-dimension, defaults to all indices
xind : Range, optional
    Range of indices to select in the x-dimension, defaults to all indices
yind : Range, optional
    Range of indices to select in the y-dimension, defaults to all indices
zind : Range, optional
    Range of indices to select in the z-dimension, defaults to all indices
info : Bool, default false
    If set to true, print extra information when reading data
"""
function boutcollect(outputs::BoutFiles, varname::String;
        tind::IntOrRange=ibegin:iend,
        xind::IntOrRange=ibegin:iend,
        yind::IntOrRange=ibegin:iend,
        zind::IntOrRange=ibegin:iend,
        info::Bool=false)

    if !(varname ∈ outputs.data_files[1])
        error(varname * " is not present in dump files")
    end

    nt = outputs.nt
    nx = outputs.nx
    ny = outputs.ny
    nz = outputs.nz
    mxg = outputs.mxg
    myg = outputs.myg
    is_doublenull = outputs.is_doublenull
    ny_inner = outputs.ny_inner
    nxpe = outputs.nxpe
    nype = outputs.nype
    mxsub = outputs.mxsub
    mysub = outputs.mysub
    yproc_upper_target = outputs.yproc_upper_target

    if length(outputs.data_files) == 1
        var = outputs.data_files[1][varname]
        dimensions = dimnames(var)

        if dimensions == ()
            # scalar
            ranges = ()
        elseif dimensions == ("t")
            # time-dependent scalar
            ranges = (tind)
        elseif dimensions == ("x", "y")
            # Field2D
            ranges = (xind, yind)
        elseif dimensions == ("x", "z")
            # FieldPerp
            ranges = (xind, zind)
        elseif dimensions == ("t", "x", "y")
            # Evolving Field2D
            ranges = (tind, xind, yind)
        elseif dimensions == ("t", "x", "z")
            # Evolving FieldPerp
            ranges = (tind, xind, zind)
        elseif dimensions == ("x", "y", "z")
            # Field3D
            ranges = (xind, yind, zind)
        elseif dimensions == ("t", "x", "y", "z")
            # Evolving Field3D
            ranges = (tind, xind, yind, zind)
        else
            ranges = nothing
        end

        return BoutArray(var.var[ranges], var.attrib, dimensions)
    end

    data_files = outputs.data_files
    nfiles = length(data_files)
    f = data_files[1]
    dimensions = dimnames(f[varname])
    ndims = length(dimensions)
    var_attributes = f[varname].attrib

    if !any(dim ∈ dimensions for dim ∈ ("x", "y", "z"))
        # No spatial dependence, so just read from first file
        if "t" ∈ dimensions
            if dimensions != ("t",)
                error(varname * " has a 't' dimension, but it is not the only "
                      * "dimension in non-spatial dimensions="
                      * string(dimensions))
            end
            return BoutArray(f[varname][tind], var_attributes, dimensions)
        else
            # No space or time dimensions, so no slicing
            data = f[varname][:]
        end
        return BoutArray(data, var_attributes, dimensions)
    end

    function converttonicerange(s, n)
        if hasproperty(s, :start)
            if s.start == ibegin
                start = 1
            else
                start = s.start
            end
        else
            start = 1
        end
        if hasproperty(s, :stop)
            if s.stop == iend
                stop = n
            else
                stop = s.stop
            end
        else
            stop = n
        end
        if hasproperty(s, :step)
            step = s.step
        else
            step = 1
        end

        size = ceil(Int, (stop + 1 - start) / step)

        return start:step:stop, size
    end

    tind, tsize = converttonicerange(tind, nt)
    xind, xsize = converttonicerange(xind, nx)
    yind, ysize = converttonicerange(yind, ny)
    zind, zsize = converttonicerange(zind, nz)

    # Create a list with the size of each dimension
    sizes = Dict("t"=>tsize, "x"=>xsize, "y"=>ysize, "z"=>zsize)
    result_dims = Tuple(sizes[d] for d ∈ dimensions)
    ndims = length(result_dims)

    if dimensions ∈ (("z", "x", "t"), ("z", "x"))
        is_fieldperp = true
        yindex_global = nothing
        # The pe_yind that this FieldPerp is going to be read from
        fieldperp_yproc = nothing
    else
        is_fieldperp = false
    end

    function getlocalxyinds(i)
        # Get X and Y processor indices
        pe_yind = i ÷ nxpe
        pe_xind = i % nxpe

        inrange = true

        # Get local ranges
        xstart = xind.start - pe_xind * mxsub
        xstop = xind.stop - pe_xind * mxsub

        # Check inner x boundary
        if pe_xind == 0
            # Keeping inner boundary
            if xstop < 1
                inrange = false
            end
            if xstart < 1
                xstart = 1
            end
        else
            if xstop < mxg + 1
                inrange = false
            end
            if xstart < mxg + 1
                xstart = mxg + 1
            end
        end

        # Outer x boundary
        if pe_xind == nxpe - 1
            # Keeping outer boundary
            if xstart > mxsub + 2*mxg
                inrange = false
            end
            if xstop > mxsub + 2 * mxg
                xstop = mxsub + 2 * mxg
            end
        end

        # Get local ranges
        ystart = yind.start - pe_yind * mysub
        ystop = yind.stop - pe_yind * mysub

        # Check lower y boundary
        if pe_yind == 0
            # Keeping lower boundary
            if ystop < 1
                inrange = false
            end
            if ystart < 1
                ystart = 1
            end
        else
            if ystop < myg + 1
                inrange = false
            end
            if ystart < myg + 1
                ystart = myg + 1
            end
        end
        # and lower y boundary at upper target
        if yproc_upper_target != nothing && pe_yind - 1 == yproc_upper_target
            ystart = ystart - myg
        end

        # Upper y boundary
        if pe_yind == nype - 1
            # Keeping upper boundary
            if ystart > mysub + 2 * myg
                inrange = false
            end
            if ystop > mysub + 2 * myg
                ystop = mysub + 2 * myg
            end
        else
            if ystart > mysub + myg
                inrange = false
            end
            if ystop > mysub + myg
                ystop = mysub + myg
            end
            # upper y boundary at upper target
            if yproc_upper_target != nothing && pe_yind == yproc_upper_target
                ystop = ystop + myg
            end
        end

        return xstart, xstop, ystart, ystop, inrange
    end

    function getlocalslices(i)
        xstart, xstop, ystart, ystop, _ = getlocalxyinds(i)
        ranges = Dict("t" => tind, "x" => xstart:xstop, "y" => ystart:ystop,
                      "z" => zind)
        return collect(ranges[d] for d ∈ dimensions)
    end

    function getglobalxyinds(i)
        # Get X and Y processor indices
        pe_yind = i ÷ nxpe
        pe_xind = i % nxpe

        xstart, xstop, ystart, ystop, _ = getlocalxyinds(i)

        xgstart = xstart + pe_xind * mxsub
        xgstop = xstop + pe_xind * mxsub

        ygstart = ystart + pe_yind * mysub
        ygstop = ystop + pe_yind * mysub
        if yproc_upper_target != nothing && pe_yind > yproc_upper_target
            ygstart = ygstart + 2 * myg
            ygstop = ygstop + 2 * myg
        end

        return xgstart, xgstop, ygstart, ygstop
    end

    function getglobalslices(i)
        xgstart, xgstop, ygstart, ygstop = getglobalxyinds(i)
        ranges = Dict("t" => tind, "x" => xgstart:xgstop,
                      "y" => ygstart:ygstop, "z" => zind)
        return collect(ranges[d] for d ∈ dimensions)
    end

    function getdata(i, f)
        # Get X and Y processor indices
        pe_yind = i ÷ nxpe
        pe_xind = i % nxpe

        # get indices
        xstart, xstop, ystart, ystop, inrange = getlocalxyinds(i)
        xgstart, xgstop, ygstart, ygstop = getglobalxyinds(i)

        if !inrange
            return nothing # Don't need this file
        end

        if info
            println("Reading from processor $i : [$xstart:$xstop, "
                    * "$ystart:$ystop] -> [$xgstart:$xgstop, "
                    * "$ygstart:$ygstop]")
        end
        if is_fieldperp
            temp_yindex = var_attributes["yindex_global"]
            if temp_yindex < 0
                # No data for FieldPerp on this processor
                return nothing
            end
        end

        return f[varname].var[getlocalslices(i)...]
    end

    function checkfieldperp!(i, f, yindex_global, var_attributes,
                             fieldperp_yproc)
        # FieldPerp should only be defined on processors which contain its
        # yindex_global
        pe_yind = i ÷ nxpe

        f_attributes = f[varname].attribs
        temp_yindex = f_attributes["yindex_global"]

        if temp_yindex >= 0
            if yindex_global == nothing
                yindex_global = temp_yindex

                # we have found a file with containing the FieldPerp, get the
                # attributes from here
                var_attributes = f_attributes
            end

            if temp_yindex != yindex_global
                error("Multiple global yindices $temp_yindex and "
                      * "$yindex_global found for FieldPerp $varname")
            end

            # Check we only read from one pe_yind
            if !(fieldperp_yproc == nothing || fieldperp_yproc == pe_yind)
                error("FieldPerp $varname found at multiple yproc indices "
                      * "$fieldperp_yproc and $pe_yind")
            end
            fieldperp_yproc = pe_yind
        else
            error("get_data() returned data, but FieldPerp $varname should "
                  * "not have data")
        end
    end

    data = zeros(Float64, result_dims...)

    for (j, f) ∈ enumerate(data_files)
        # Processor indices are 0-based not 1-based
        i = j - 1

        result = getdata(i, f)
        if result != nothing
            if is_fieldperp
                checkfieldperp!(i, f, yindex_global, var_attributes,
                                fieldperp_yproc)
            end

            data[getglobalslices(i)...] = result
        end
    end

    # if a step was requested in x or y, need to apply it here
    if xind.step != 1 || yind.step != 1
        ranges = Dict("t" => ibegin:iend, "x" => ibegin:xind.step:iend,
                      "y" => ibegin:yind.step:iend, "z" => ibegin:iend)
        data = data[(ranges[d] for d ∈ dimensions)]
    end

    return BoutArray(data, var_attributes, dimensions)
end

end
