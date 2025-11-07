function ODE!(du, u, p, t)

    fill!(du, 0.0)  # Reset to zero, just like initializing with zeros

    rhs_metabolites!(du, u, p)
    
    rhs_sequence!(du, u, p)

    nothing 
end 