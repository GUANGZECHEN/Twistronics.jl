using LinearAlgebra
using DelimitedFiles
using Statistics
using Dates

using Random
using SparseArrays
using Arpack

using PyCall
const plt=pyimport_conda("matplotlib.pyplot","matplotlib")
const colormap=pyimport_conda("matplotlib.colors","matplotlib")


include("Geometry_twisted.jl")
include("tight-binding_twisted.jl")
include("quantities_twisted.jl")
include("valley_operator_twisted.jl")

