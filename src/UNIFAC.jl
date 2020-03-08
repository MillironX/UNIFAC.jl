# Thomas A. Christensen II
# Spring 2020
# UNIFAC Solver

module UNIFAC

export unifac

"""
    unifac(ν, x, T)

Calculates the activity coefficients of the chemical species described by the
UNIFAC groups in `ν`, at the liquid mole fractions `x` and temperature `T` in
Kelvins.

`ν` is constructed as a 56 x n `Array{Int}`, such that `ν[i,j]` is the number
of times subgroup `i` appears in molecule `j`.

# Example

```julia-repl
julia> ν = zeros(56, 2);
julia> ν[1, 1] = 2;
julia> ν[2, 1] = 1;
julia> ν[33, 1] = 1;
julia> ν[1, 2] = 2;
julia> ν[2, 2] = 5;
julia> unifac(ν, [0.4 0.6], 308.15)
1×2 Array{Float64,2}:
 1.133  1.047
```
"""
function unifac(ν, x, T)
# unifac takes three arguments:
# ν (that's the Greek lowercase nu) is a 56 x n Integer array where n is the
#   number of species in the system. The number at ν[k,i] indicates how many
#   instances of subgroup k are present in a single molecule of species i
# x is a 1 x n Float64 array that contains the liquid species mole fractions
# T is the temperature of the system in Kelvins

# Most subgroup parameters and interactions are taken from Smith, Van Ness, &
# Abbott, _Introduction to Chemical Engineering Thermodynamics_, 7th Ed. pp.
# 791-7. Some subgroup parameters are taken from Dr. Aston's given table, and
# some of the interactions for those species are taken from
# http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html

    # Declare the UNIFAC subgroup parameters
    R = Array{Float64}(undef, 56, 1)
    R[1]  = 0.9011
    R[2]  = 0.6744
    R[3]  = 0.4469
    R[4]  = 0.2195
    R[10] = 0.5313
    R[12] = 1.2663
    R[13] = 1.0396
    R[15] = 1.0000
    R[17] = 0.9200
    R[19] = 1.6724
    R[20] = 1.4457
    R[25] = 1.1450
    R[26] = 0.9183
    R[27] = 0.6908
    R[32] = 1.4337
    R[33] = 1.2070
    R[34] = 0.9795
    R[41] = 1.8701
    R[42] = 1.6434
    R[56] = 1.7818

    Q = Array{Float64}(undef, 56, 1)
    Q[1]  = 0.848
    Q[2]  = 0.540
    Q[3]  = 0.228
    Q[4]  = 0.000
    Q[10] = 0.400
    Q[12] = 0.968
    Q[13] = 0.660
    Q[15] = 1.200
    Q[17] = 1.400
    Q[19] = 1.488
    Q[20] = 1.180
    Q[25] = 1.088
    Q[26] = 0.780
    Q[27] = 0.468
    Q[32] = 1.244
    Q[33] = 0.936
    Q[34] = 0.624
    Q[41] = 1.724
    Q[42] = 1.416
    Q[56] = 1.560

    # Declare the subgroup interaction parameters
    a = zeros(56, 56)
    a[1:4, 10]    .= 61.13
    a[1:4, 12:13] .= 76.50
    a[1:4, 15]    .= 986.50
    a[1:4, 17]    .= 1318.00
    a[1:4, 19:20] .= 476.40
    a[1:4, 25:27] .= 251.50
    a[1:4, 32:34] .= 255.70
    a[1:4, 41:42] .= 597.00
    a[1:4, 56]    .= 661.50

    a[10, 1:4]   .= -11.12
    a[10, 12:13] .= 167.00
    a[10, 15]     = 636.10
    a[10, 17]     = 903.80
    a[10, 19:20] .= 25.77
    a[10, 25:27] .= 32.14
    a[10, 32:34] .= 122.80
    a[10, 41:42] .= 212.50
    a[10, 56]     = 168.00

    a[12:13, 1:4]   .= -69.70
    a[12:13, 10]    .= -146.80
    a[12:13, 15]    .= 803.20
    a[12:13, 17]    .= 5695.00
    a[12:13, 19:20] .= -52.10
    a[12:13, 25:27] .= 213.10
    a[12:13, 32:34] .= -49.29
    a[12:13, 41:42] .= 6096.00
    a[12:13, 56]    .= 3629.0

    a[15, 1:4]   .= 156.40
    a[15, 10]     = 89.60
    a[15, 12:13] .= 25.82
    a[15, 17]     = 353.50
    a[15, 19:20] .= 84.00
    a[15, 25:27] .= 28.06
    a[15, 32:34] .= 42.70
    a[15, 41:42] .= 6.712
    a[15, 56]     = 256.5

    a[17, 1:4]   .= 300.00
    a[17, 10]     = 362.30
    a[17, 12:13] .= 377.60
    a[17, 15]     = -229.10
    a[17, 19:20] .= -195.40
    a[17, 25:27] .= 540.50
    a[17, 32:34] .= 168.00
    a[17, 41:42] .= 112.60
    a[17, 56]     = 220.60

    a[19:20, 1:4]   .= 26.79
    a[19:20, 10]    .= 140.10
    a[19:20, 12:13] .= 365.80
    a[19:20, 15]    .= 164.50
    a[19:20, 17]    .= 472.50
    a[19:20, 25:27] .= -103.60
    a[19:20, 32:34] .= -174.20
    a[19:20, 41:42] .= 481.70
    a[19:20, 56]    .= 137.50

    a[25:27, 1:4]   .= 83.36
    a[25:27, 10]    .= 52.13
    a[25:27, 12:13] .= 65.69
    a[25:27, 15]    .= 237.70
    a[25:27, 17]    .= -314.70
    a[25:27, 19:20] .= 191.10
    a[25:27, 32:34] .= 251.50
    a[25:27, 41:42] .= -18.51
    a[25:27, 56]    .= 95.108

    a[32:34, 1:4]   .= 65.33
    a[32:34, 10]    .= -22.31
    a[32:34, 12:13] .= 223.00
    a[32:34, 15]    .= -150.00
    a[32:34, 17]    .= -448.20
    a[32:34, 19:20] .= 394.60
    a[32:34, 25:27] .= -56.08
    a[32:34, 41:42] .= 147.10

    a[41:42, 1:4]   .= 24.82
    a[41:42, 10]    .= -22.97
    a[41:42, 12:13] .= -138.40
    a[41:42, 15]    .= 185.40
    a[41:42, 17]    .= 242.80
    a[41:42, 19:20] .= -287.50
    a[41:42, 25:27] .= 38.81
    a[41:42, 32:34] .= -108.50
    a[41:41, 56]    .= -0.5150

    a[56, 1:4]   .= -32.690
    a[56, 10]    .= 10.380
    a[56, 12:13] .= -97.050
    a[56, 15]    .= 261.60
    a[56, 17]    .= 417.90
    a[56, 19:20] .= -142.60
    a[56, 25:27] .= -94.490
    a[56, 41:42] .= 0.28270

    # Calculate r
    r = zeros(1, size(ν, 2))
    for i in 1:size(ν, 2)
        r[i] = sum(ν[:, i] .* R)
    end

    # Calculate q
    q = zeros(1, size(ν, 2))
    for i in 1:size(ν, 2)
        q[i] = sum(ν[:, i] .* Q)
    end

    # Calculate e
    e = zeros(56, size(ν, 2))
    for i in 1:size(ν,2 )
        e[:, i] = (ν[:, i] .* Q) ./ q[i]
    end

    # Calculate tau
    τ = exp.(-a ./ T)

    # Calculate beta
    β = zeros(56, size(ν, 2))
    for i in 1:size(ν, 2)
        for k in 1:56
            β[k, i] = sum(e[:, i] .* τ[:, k])
        end
    end

    # Calculate Theta
    Θ = zeros(56, 1)
    for k in 1:56
        Θ[k] = sum(x .* q .* e[k, :]') / sum(x .* q)
    end

    # Calculate s
    s = zeros(56, 1)
    for k in 1:56
        s[k] = sum(Θ .* τ[:, k])
    end

    # Calculate L and J
    L = r ./ sum(r .* x)
    J = q ./ sum(q .* x)

    # Calculate the natural log of the cumulative and residual activity coefficients
    γ_C = 1 .- J .+ log.(J) .- 5 .* q .* (1 .- (J ./ L) .+ log.(J ./ L))
    γ_R = zeros(1, size(ν, 2))
    for i in 1:size(ν, 2)
        γ_R[i] = q[i] * (1 - sum(Θ .* β[:, i] ./ s .- e[:, i] .* log.(β[:, i] ./ s)))
    end

	# Return the activity coefficients
    γ = exp.(γ_C .+ γ_R)

end # function

end # module
