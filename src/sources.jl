# CMB and SZ have nontrivial frequency dependence for σ(ν)
const h = 6.62607e-34 # Plancks constant SI
const k_B = 1.38065e-23 # Boltzman's constant, SI
const T_CMB = 2.725 # CMB Temperature, Kelvins
const c = 3e8 # speed of light
const b_0 = ((2 * k_B^3  * T_CMB^2) / (h * c)^2) # defined in paper

function color_correction(frequencies::Vector{T}, passband::Vector{T},
                          integration_bounds::NTuple{2, Real}, β::Real,
                          ν_r::Real) where {T<:Real}
    ν, passband = get_restricted_regions(frequencies, passband,
                                         integration_bounds)
    numerator = trapz(ν, (ν ./ ν_r).^β .* ν.^(-2) .* passband)
    denom = trapz(ν, (b_ν.(ν) / b_ν(ν_r)) .* ν.^(-2) .* passband)
    return numerator / denom
end

function center_spectral_source_eq30(
    frequencies::Vector{T}, passband::Vector{T},
    integration_bounds::NTuple{2, Real}, β::Real; α::Real=0.,
    ν_B=100.) where {T<:Real} #ν and ν_B must match units!
    # passband weighted mean calculation + spectral index usage
    # ν is our frequency range
    # In the paper the -2 exponent goes with α but in this case we put
    # it with β for efficiency purposes since beam-filling sources set α = 0.
    ν, passband = get_restricted_regions(frequencies, passband,
                                         integration_bounds)
    numerator = trapz(ν, ν.^(β-2) .* ν.^α .* passband)
    denom = (ν_B ^ α) * trapz(ν, passband)
    # the equation in the paper is messed up I believe (wrong exponent)
    return (numerator / denom) ^ (1 / (β - 2))
end

function gamma_eq34(
    frequencies::Vector{T}, passband::Vector{T},
    integration_bounds::NTuple{2, Real}, β::Real; α::Real=0.,
    ν_B=100., ν_R=ν_B, Ω_B=200.) where {T<:Real} #Ω_B in nsr
    # ν, ν_R, and ν_B must match units!
    # set ν_R to just the central frequency of the band
    ν_R = center_spectral_source_eq30(
        frequencies, passband, integration_bounds, β, α=α, ν_B=ν_B)
    # unit-converted prefactor computed from Mathematica earlier
    prefactor = 32.5483
    ν, passband = get_restricted_regions(frequencies, passband,
                                         integration_bounds)
    numerator = prefactor * trapz(ν, ν.^(β) .* ν.^(α-2) .* passband)
    denom = (ν_R ^ β) * Ω_B * 1e-9 * (ν_B ^ α) * trapz(ν, passband)
    return (numerator / denom)
end

function center_CMB_eq18(
    frequencies::Vector{T}, passband::Vector{T},
    integration_bounds::NTuple{2, Real}) where {T<:Real}
    # passband weighted mean calculation + spectral index usage
    # ν is our frequency range
    # Here we have to solve for 'x' and thus 'ν_eff' from the equation,
    # which is nontrivial but can be done numerically
    ν, passband = get_restricted_regions(frequencies, passband,
                                         integration_bounds)
    numerator = trapz(ν, t_ν.(ν) .* passband)
    denom = trapz(ν, passband)
    RHS = (numerator / denom)
    # find the zero of LHS - RHS which should always be between the
    # integration bounds!
    find_zero(ν_eff->t_ν(ν_eff) - RHS, integration_bounds)
end

function center_source_eq16(
    frequencies::Vector{T}, passband::Vector{T},
    integration_bounds::NTuple{2, Real}, Iν::Function,
    Iν_args...) where {T<:Real}

    ν, passband = get_restricted_regions(frequencies, passband,
                                         integration_bounds)
    numerator = trapz(ν, Iν.(ν, Iν_args...) .* (ν).^(-2) .* passband)
    denom = trapz(ν, passband)
    RHS = (numerator / denom)
    # find the zero of LHS - RHS which should always be between the
    # integration bounds!
    # simplify t / b = 1 / ν^2 (constants cancel with RHS)
    # may have to put ν_eff in Ghz for this to work??
    eq_func = ν_eff->((Iν(ν_eff, Iν_args...) / ((ν_eff)^2)) - RHS)
    find_zero(eq_func, integration_bounds)
end

# convert from hz to ghz!
get_x(ν_ghz::Real, T::Real=T_CMB) = h * (ν_ghz*1e9) / (k_B * T)

function CMB_σ(ν::Real)
    x = get_x(ν)
    return  ((x^4. * exp(x)) / (exp(x) - 1.)^2.)
end
b_ν(ν::Real) = b_0 * CMB_σ(ν)

function B_ν(ν::Real, T::Real=T_CMB)
    x = get_x(ν, T)
    return b_0  * (x^3) / (exp(x) - 1)
end

function t_ν(ν::Real)
    x = get_x(ν)
    return ((x^2. * exp(x)) / (exp(x) - 1.)^2.)
end

function dust_Iν(ν_ghz::Real, β::Real, T::Real)
    ν = ν_ghz
    ν_rd = 364.2 #GHz.
    term1 = (2 * ν^2 * k_B / (c^2)) * (ν / ν_rd) ^ (β - 2)
    term2 = B_ν(ν, T) / B_ν(ν_rd, T)
    return term1 * term2
end

function tSZ_Iν(ν::Real, y::Real=1e-4)
    x = get_x(ν)
    term1 = b_0 * T_CMB * ((x^4 * exp(x)) / (exp(x) - 1)^2)
    term2 = x * ((exp(x) + 1) / (exp(x) - 1)) - 4
    return term1 * term2 * y
end

function spectral_Iν(ν::Real, β::Real, ν_r::Real=100., I_r::Real=10.)
    return I_r * (ν / ν_r) ^ β
end


function SZ_σ(ν::Real)
    x = get_x(ν)
    return (x * coth(x / 2.)) - 4.
end

function spectral_σ(i::Real)
    return (x->x^i)
end

function spectral_σ_point_source(i::Real)
    return (x->x^(i - 2.))
end

# Check: CMB should stay the same as before!
function center_CMB(frequencies::Vector{T}, passband::Vector{T},
                    integration_bounds::NTuple{2, Real}) where {T<:Real}
    freq_range, passband = get_restricted_regions(frequencies, passband,
                                                  integration_bounds)
    σ = CMB_σ.(freq_range)
    numerator = trapz(freq_range, σ .* passband .* (freq_range.^(-1)))
    denom = trapz(freq_range, σ .* passband .* freq_range.^(-2.))
    return numerator / denom
end
