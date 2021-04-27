module SquashOutput

export squashoutput

using NCDatasets
using NCDatasets: Attributes

using BoutOutputs
using BoutOutputs: BoutArray, BoutScalar

"""
Combine BOUT++ output into a single file

Parameters
----------
datadir : String, default "."
    Directory to read from.
outputname : String, default "boutdata.nc"
    Name of the single output file, which will be written in datadir.
prefix : String, default "BOUT.dmp"
    Prefix for file names. E.g. by default BOUT++ output files are named like
    "BOUT.dmp.<i>.nc" where <i> is an integer.
complevel : Int, default 0
    Compression level. 0 is no compression. 1 should be fastest compression,
    and 9 should yield highest compression level.
singleprecision : Bool, default false
    Convert floats to single precision
kwargs
    Extra keyword arguments are passed to boutcollect
"""
function squashoutput(datadir::String="."; outputname::String="boutdata.nc",
                      prefix::String="BOUT.dmp", complevel::Int=0,
                      singleprecision=false, kwargs...)
    boutfiles = openboutdata(datadir, prefix=prefix)

    # Create output file, copying file-level attributes from root file
    outputfile = createoutputfile(joinpath(datadir, outputname), boutfiles.nx,
                                  boutfiles.ny, boutfiles.nz,
                                  boutfiles.data_files[1].attrib)

    for varname ∈ keys(boutfiles.data_files[1])
        println(varname)
        data = boutcollect(boutfiles, varname, kwargs...)
        addvariable(outputfile, varname, data, complevel, singleprecision)
    end
end

function createoutputfile(path::String, nx::Int, ny::Int, nz::Int,
                          attributes::Attributes)
    if isfile(path)
        error("$path already exsists")
    end

    file = NCDataset(path, "c")

    # Use length=Inf to make time dimension 'unlimited'
    defDim(file, "t", Inf)
    defDim(file, "x", nx)
    defDim(file, "y", ny)
    defDim(file, "z", nz)

    file.attrib = attributes

    return file
end

function addvariable(file::NCDataset, name::String, data::BoutArray,
                     complevel::Int, singleprecision::Bool)
    outtype = typeof(data.data[1])
    if singleprecision && outtype == Float64
        outtype = Float32
    end
    v = defVar(file, name, outtype, data.dimensions, attrib=data.attributes,
               deflatelevel=complevel)

    v[:] = data.data
end

function addvariable(file::NCDataset, name::String, data::BoutScalar,
                     complevel::Int, singleprecision::Bool)
    outtype = typeof(data.data)
    if singleprecision && outtype == Float64
        outtype = Float32
    end
    v = defVar(file, name, outtype, (), attrib=data.attributes)

    v[:] = data.data
end

end
