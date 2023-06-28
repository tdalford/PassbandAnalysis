const ORIGINAL_CAL_FACTOR = 1.0239
const SET_CAL_FACTOR = 1.02
const Z_FOCUS_MM = 25.0 #16.7
const Z_UNCERT_MM = 50 + 1
const IN_TO_MM = 25.4

mutable struct SimData
    # includes a matrix of correction data (over the focal plane), a vector of
    # frequencies, and the z value that this came from
    frequency_shifts::Matrix{R} where R<:Real
    amplitude_shifts::Matrix{S} where S<:Real
    frequencies::Vector{T} where T<:Real
    z_value::Real
    centers::Matrix{U} where U<:Real

    function SimData(
        frequency_shifts::Matrix{R}, amplitude_shifts::Matrix{S},
        frequencies::Vector{T}, z_value::Real,
        centers::Matrix{U}; n_freq_shifts::Integer=25) where {
            R<:Real, T<:Real, S<:Real, U<:Real}
        # make sure that the # data points in passbands is the # of frequencies
        # we want to store the passbands as columns here!
        axes(frequency_shifts, 1) == axes(frequencies, 1) ||
            error("corrections and frequencies have different numbers")

        # make sure that the number of passbands is the number of attrs
        size(frequency_shifts, 2) == n_freq_shifts || error("incorrect size of
        axis 2 of frequency_shifts. Expected 25, got $(size(frequency_shifts, 2))")

        axes(frequency_shifts) == axes(amplitude_shifts) ||
            error("amplitude shifts and frequency shifts have different axes!")

        new(frequency_shifts, amplitude_shifts, frequencies, z_value, centers)
    end

end

"""Functor to map z coordinate -> SimData"""
function (z_sim_data::Vector{SimData})(z)
    # find the z that corresponds best to this z value
    zvals = map(x->x.z_value, z_sim_data)
    # assume we're in increasing order!
    z_sim_data[findfirst(x -> x >= z, zvals)]
end

# Factor of 1.0183 / 1.0153 is from an FFT sampling error which yields a peak
# frequency that is consistently too low.
frequency_shift_conversion(freq_shifts::AbstractArray) = (1.0183 / 1.015) ./ (
    1 .+ freq_shifts)

function load_python_sim_data(filename::String)
    attr_data_vec = SimData[]
    h5open(filename, "r") do file
    centers = read(file, "centers")
        for i = 0:read(file, "n_runs")-1
            correction_data = read(file, "shift_data_$i")
            # hard coded from reading the data in.
            freq_shift_data = correction_data[1, :, :]
            # apply conversion
            freq_shift_data = frequency_shift_conversion(freq_shift_data)
            amp_data = correction_data[2, :, :]
            frequencies = read(file, "frequencies_$i")
            # convert to mm!
            z_value = read(file, "z_val_$i") * IN_TO_MM
            sim_data = SimData(freq_shift_data, amp_data, frequencies, z_value,
                               centers)
            push!(attr_data_vec, sim_data)
        end
    end
    return attr_data_vec
end

function restrict_z_data(z_sim_data::Vector{SimData})
    restricted_z_sim_data = SimData[]
    for sim_data in z_sim_data
        if (Z_FOCUS_MM - Z_UNCERT_MM  <= sim_data.z_value <=
            Z_FOCUS_MM + Z_UNCERT_MM)
            push!(restricted_z_sim_data, sim_data)
        end
    end
    return restricted_z_sim_data
end

function get_freq_shift_interpolation(
    centers::Matrix{T}, freq_shifts::Vector{U};
    extrapolation_bc=Line()) where {T<:Real, U<:Real}
    # assume the centers object is fully created!
    x = centers[1, :] |> unique
    y = centers[2, :] |> unique
    # in order to correctly interpolate the right dimensions, take the adjoint
    z = reshape(freq_shifts, length(x), length(y))'
    return LinearInterpolation((x, y), z, extrapolation_bc=extrapolation_bc)
end

get_freq_shift_interpolation(sim_data::SimData, frequency::Real) =
    get_freq_shift_interpolation(sim_data.centers, get_freq_item(
        frequency, sim_data))


"""This doesn't work..."""
function get_amplitude_interpolation(
    sim_data::SimData, extrapolation_bc=Line())
    centers = sim_data.centers
    frequencies = sim_data.frequencies
    amplitude_shifts = similar(sim_data.amplitude_shifts)
    for i in 1:size(amplitude_shifts, 2)
        shift = sim_data.amplitude_shifts[:, i]
        # normalize all of these!
        cleaned_shift = normalize(shift)
        amplitude_shifts[:, i] = cleaned_shift
    end

    # assume the centers object is fully created!
    x = centers[1, :] |> unique
    y = centers[2, :] |> unique
    # in order to correctly interpolate the right dimensions, take the adjoint
    z = reshape(amplitude_shifts, length(x), length(y), length(frequencies))
    return LinearInterpolation((x, y, frequencies), z,
                               extrapolation_bc=extrapolation_bc)
end

function get_closest_amplitude_func(sim_data::SimData, loc::NTuple{2, Real})
    # find the pair whose distance is the smallest
    distances = [hypot((loc .- center_loc)...)
                 for center_loc in eachcol(sim_data.centers)]
    return normalize(sim_data.amplitude_shifts[:, argmin(distances)])
end

function get_freq_item(frequency::Real, sim_data::SimData)
    freq_index = findfirst(f -> f >= frequency, sim_data.frequencies)
    return sim_data.frequency_shifts[freq_index, :]
end

function get_shifted_array(passband_data::PassbandArray{T},
                           sim_data::SimData, frequency::Real) where T
    corrected_bands = similar(passband_data.passbands)
    # compute the calibration factor and transfer function iteratively here!
    itp = get_freq_shift_interpolation(sim_data.centers, get_freq_item(
        frequency, sim_data))
    for i in 1:size(passband_data.passbands, 2)
        band = passband_data.passbands[:, i]
        location = (passband_data.attrs[i, "x_dists"], passband_data.attrs[
            i, "y_dists"])
        # find the frequency calibration factor from our interpolator
        cal_factor = itp(location...) / SET_CAL_FACTOR
        shifted_band = get_shifted_band(passband_data.frequencies, band,
                                        cal_factor)
        # now get the amplitude calibration at the same frequencies!
        corrected_shifted_band = get_amplitude_corrected_band(
            passband_data.frequencies, shifted_band, sim_data, location)

        corrected_bands[:, i] = corrected_shifted_band
    end
    return corrected_bands
end

function get_only_freq_shifted_array(passband_data::PassbandArray{T},
                                     sim_data::SimData, frequency::Real) where T
    corrected_bands = similar(passband_data.passbands)
    # compute the calibration factor and transfer function iteratively here!
    itp = get_freq_shift_interpolation(sim_data.centers, get_freq_item(
        frequency, sim_data))
    for i in 1:size(passband_data.passbands, 2)
        band = passband_data.passbands[:, i]
        location = (passband_data.attrs[i, "x_dists"], passband_data.attrs[
            i, "y_dists"])
        # find the frequency calibration factor from our interpolator
        cal_factor = itp(location...) / SET_CAL_FACTOR
        shifted_band = get_shifted_band(passband_data.frequencies, band,
                                        cal_factor)
        corrected_bands[:, i] = shifted_band
    end
    return corrected_bands
end

function get_rmses(passband_data::PassbandArray{<:Real},
                   z_sim_data::Vector{SimData}, center_freq::Real,
                   int_bounds::NTuple{2, Real}; spectral_index::Real=0)
    devs = Float64[]
    for sim_data in z_sim_data
        corrected_bands = get_only_freq_shifted_array(
            passband_data, sim_data, center_freq)
        corrected_centers = get_centers(
            passband_data.frequencies, corrected_bands, int_bounds,
            spectral_index=spectral_index)
        _, rms_center = weighted_rms(
            corrected_centers, get_weights(passband_data))
        push!(devs, rms_center)
    end
    devs
end

function get_corrections(sim_data::SimData)
    radaii = map(x -> hypot(x...), eachcol(sim_data.centers))
    # get the corrections at all radaii
    radius_corrections = Dict()
    for r in sort(unique(radaii))
        # get the corretions corresponding to this radius
        radius_inds = findall(x->x==r, radaii)
    radius_corrections[floor(Integer, r)] = []
        for correc in eachcol(sim_data.amplitude_shifts[:, radius_inds])
            push!(radius_corrections[floor(Integer, r)],
                  correc ./ maximum(correc))
        end
    end
    radius_corrections
end
