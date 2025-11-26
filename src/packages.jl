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
using Statistics
using DelimitedFiles
using LinearAlgebra

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
include("rhs_metabolites.jl")
include("rhs_sequence.jl")
include("parameters.jl")

