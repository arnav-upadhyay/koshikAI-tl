# Packages

## Data management
using BioSequences
using FASTX
using CSV
using DataFrames
using XLSX
using FastaIO

## Solvers
using DifferentialEquations
using Statistics
using DelimitedFiles
using LinearAlgebra

## Analysis
using Plots
using LaTeXStrings  
using Colors
using JET
using Dates

# Files 

include("datatypes.jl")
include("efficiency_factors.jl")
include("functions.jl")
include("rhs_metabolites.jl")
include("rhs_sequence.jl")
include("ODE.jl")
include("parameters.jl")

