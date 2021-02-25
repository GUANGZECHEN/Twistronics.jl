__precompile__()

module Twistronics

using LinearAlgebra
using DelimitedFiles
using Statistics
using Dates

using Random
using SparseArrays
using Arpack

using PyCall
const plt=PyNULL()
const colormap=PyNULL()
function __init__()
    copy!(plt,pyimport_conda("matplotlib.pyplot","matplotlib"))
    copy!(colormap,pyimport_conda("matplotlib.colors","matplotlib"))
end

include("Geometry_twisted.jl")
include("tight-binding_twisted.jl")
include("quantities_twisted.jl")
include("valley_operator_twisted.jl")

end
