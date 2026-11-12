# Packages

## Data management
using BioSequences
using FASTX
using CSV
using DataFrames
using XLSX
using FastaIO
using UnPack

## Solvers
using DifferentialEquations
using SciMLBase: successful_retcode
using Statistics
using DelimitedFiles
using LinearAlgebra
using Printf

## Analysis
using Plots
using LaTeXStrings  
using Colors
using Dates
using JLD2
using Measures

# Code optimization
using JET
using BenchmarkTools
using Base.Threads
println("🔧 Julia started with $(Threads.nthreads()) threads")

# Files 
include("datatypes.jl")
include("functions.jl")
include("ode_diagnostics.jl")
include("rhs_metabolites.jl")
println("  Leu-RS: competitive (LEU_RS_COMPETITIVE=", LEU_RS_COMPETITIVE, ")")
include("rhs_sequence.jl")
include("parameters.jl")
include("stall_freeze.jl")
include("reduced_cell_load.jl")

