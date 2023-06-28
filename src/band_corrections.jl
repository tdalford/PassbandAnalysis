"""interpolation using flat extrapolation scheme"""
function interpolate(x::Vector{S}, y::Vector{T},
                     x_new::Vector{U}) where {S, T, U}
    itp = LinearInterpolation(x, y, extrapolation_bc=Flat())
    return map(itp, x_new)
end

"""Shift the band due to a multiplicative frequency correction (cal_factor)"""
function get_shifted_band(frequencies::Vector{T}, passband::Vector{T},
                          cal_factor::Real) where T
    shifted_frequencies = frequencies .* cal_factor
    return interpolate(shifted_frequencies, passband, frequencies)
end


"""Interpolate an amplitude correction to our band frequencies and apply it."""
function apply_amplitude_correction(
    frequencies::Vector, passband::Vector, correction::Vector,
    correction_frequencies::Vector)
    full_amp_correction = interpolate(correction_frequencies, correction,
                                      frequencies)
    return correct_band(passband, full_amp_correction)
end

function apply_correction_to_bands(
    frequencies::AbstractVector, passbands::AbstractMatrix,
    correction::AbstractVector, correction_freqs::AbstractVector)
    return hcat([normalize(
        frequencies, apply_amplitude_correction(
            frequencies, passbands[:, i], correction, correction_freqs))
                 for i in 1:size(passbands, 2)]...)
end

function apply_corrections(
    passband_data::PassbandArray{<:Real}, bands::Matrix{<:Real},
    corrections::AbstractVector, data_freqs::AbstractVector)
    for (correction, correction_freqs) in zip(corrections, data_freqs)
        bands = apply_correction_to_bands(
            passband_data.frequencies, bands, correction, correction_freqs)
    end
    arr = copy(passband_data); arr.passbands = bands;
    return arr
end

"""Apply a division correction and renormalize band"""
correct_band(passband::Vector{T}, correction::Vector{T}) where T =
    normalize(passband ./ correction)
