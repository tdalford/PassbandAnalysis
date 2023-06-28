using HDF5, DataFrames, Interpolations, Trapz, Printf, Roots, Statistics, StatsBase

const CENTROID_ERROR_MARGIN = 8 #mm
const MAX_CENTROID_ERROR_MARGIN = 5 #mm

# Declare the type PassbandArray, which holds a bunch of passbands along with
# frequencies for the bands and a DataFrame of attributes.
mutable struct PassbandArray{T<:Real}
    passbands::Matrix{T}
    frequencies::Vector{T}
    # includes things like channel, center, width, SNR that we can access!
    attrs::DataFrame
    label::String
    # create an constructor for this that checks that the length of passbands
    # is the length of attrs, and that the width of passbands is the width of
    # attrs
    function PassbandArray{T}(
        passbands::Matrix{T}, frequencies::Vector{T}, attrs::DataFrame,
        label::String; set_distances::Bool=true) where T<:AbstractFloat
        # make sure that the # data points in passbands is the # of frequencies
        # we want to store the passbands as columns here!
        axes(passbands, 1) == axes(frequencies, 1) ||
            error("passbands and frequencies have different numbers")

        # make sure that the number of passbands is the number of attrs
        axes(passbands, 2) == axes(attrs, 1) ||
            error("passbands and attrs have different numbers")

        passband_array = new(passbands, frequencies, attrs, label)
        # Setter if the centroid distances not here!
        if set_distances
            set_dist_from_centroid!(passband_array)
        end
        passband_array
    end

end

function PassbandArray(passbands::Matrix{T}, frequencies::Vector{T},
                       attrs::DataFrame, label::String;
                       set_distances::Bool=true) where T
    return PassbandArray{T}(passbands, frequencies, attrs, label,
                            set_distances=set_distances)
end

"""Copy over all the attrs"""
Base.copy(pa::PassbandArray) = PassbandArray(
    copy(pa.passbands), copy(pa.frequencies), copy(pa.attrs), pa.label,
    set_distances=false)

Base.length(pa::PassbandArray) = size(pa.attrs, 1)

get_centroid(position_data::AbstractArray) = mean(position_data, dims=1)

function safe_stderr(vals::AbstractArray)
    if size(vals, 1) <= 1
        return MAX_CENTROID_ERROR_MARGIN
    end
    #return 8
    return max(std(vals), CENTROID_ERROR_MARGIN)
end

get_centroid_error(position_data::AbstractMatrix) = map(
    safe_stderr, eachcol(position_data))

function set_dist_from_centroid!(passband_data::PassbandArray{<:Real})
    position_data = Matrix(passband_data.attrs[!, ["x_loc", "y_loc"]])
    centroid = get_centroid(position_data)
    centroid_σ = get_centroid_error(position_data)
    dists = position_data .- centroid
    passband_data.attrs[!, "x_dists"] = dists[:, 1]
    passband_data.attrs[!, "y_dists"] = dists[:, 2]

    # add the errors in how well we know the centroid distances
    passband_data.attrs[!, "x_dists_σ"] .= centroid_σ[1, :]
    passband_data.attrs[!, "y_dists_σ"] .= centroid_σ[2, :]
    return
end

"""Concatentate a vector of PassbandArrays into one giant PassbandArray

This leaves the centroid  / distance data as is, doesn't create any new ones
since each array hopefully has its own centroid in it!
"""
function concat(passband_arrays::Vector{PassbandArray{T}}) where T
    # We just need to make sure that the frequencies are the same in each
    frequencies = passband_arrays[1].frequencies
    total_passbands = passband_arrays[1].passbands
    total_attrs = passband_arrays[1].attrs
    label = passband_arrays[1].label

    for passband_array in passband_arrays[2:end]
        all(passband_array.frequencies .≈ frequencies) || error(
        "frequencies need to match between passband arrays!")
        total_passbands = [total_passbands;; passband_array.passbands]
        total_attrs = vcat(total_attrs, passband_array.attrs)
    end
    # Don't create a new centroid!
    PassbandArray(total_passbands, frequencies, total_attrs, label,
                  set_distances=false)

end

###############################################################################
# How do I import this data from the python code?
#
# 1. When saving the data, the passbands can easily form a matrix. Frequencies
# as well, etc
#
# 2. Then convert the attrs into a DataFrame and save as a CSV
# 3. Hard part is just setting the x locations and y locations to the data.
# However, since these are all just dependent on detector, this should be easy
# to work with.
###############################################################################

function load_python_passband_data(
    filename::String; label::String="", set_cal_factor::Bool=true,
    load_attrs::Vector=[], swap_xy::Bool=false)

    passband_data_vec = PassbandArray{Float64}[]
    h5open(filename, "r") do file
        # convert from hz to Ghz
        frequencies = read(file, "frequencies") / 1e9
        # Apply the frequency cal factor conversion!
        if (set_cal_factor)
            frequencies = frequencies * (SET_CAL_FACTOR / ORIGINAL_CAL_FACTOR)
        end

        x_loc = read(file, "x_loc")
        y_loc = read(file, "y_loc")
        key_order = ["det", "center", "width", "SNR", "lower_edge",
                     "upper_edge"]
        for i = 0:read(file, "n_sets")-1
            passbands = read(file, "passband_data_$i")
            if size(passbands)[1] == 0
                continue
            end
            # the data is the transposed of what we want!
            attrs = DataFrame(read(file, "attrs_data_$i")', key_order)
            # convert from python zero-indexing to one-indexing
            loc_accesser = convert(Vector{Int64}, attrs[!, "det"]) .+ 1
            # we want to swap the x and y axes for the position data to match
            # up the axes we choose for the FTS simulations.
            if swap_xy
                attrs[!, "y_loc"] = x_loc[loc_accesser]
                attrs[!, "x_loc"] = y_loc[loc_accesser]
            else
                attrs[!, "x_loc"] = x_loc[loc_accesser]
                attrs[!, "y_loc"] = y_loc[loc_accesser]
            end
            # load all the extra attrs in as well
            for at in load_attrs
                attrs[!, at] = read(file, at)[loc_accesser]
            end
            # set a data element to the run number
            attrs[!, "run_num"] .= i
            passband_array = PassbandArray(passbands, frequencies, attrs,
                                           label)
            push!(passband_data_vec, passband_array)
        end
    end
    return passband_data_vec
end


"""assume arr is sorted"""
function find_subarray_limits(
    arr::Vector{T}, start_val::Real, end_val::Real) where T<:Real
    start_val_index = findfirst(x -> x >= start_val, arr)
    end_val_index = findlast(x -> x <= end_val, arr)
    #start_val_index, end_val_index
    start_val_index, end_val_index
end

function get_restricted_regions(frequencies::Vector{T}, passband::Vector{T},
                                integration_bounds::NTuple{2, Real}) where T
    start_freq_index, end_freq_index = find_subarray_limits(
        frequencies, integration_bounds...)
    freq_range = frequencies[start_freq_index : end_freq_index]
    passband = passband[start_freq_index : end_freq_index]
    return freq_range, passband
end

function center(frequencies::Vector{T}, passband::Vector{T},
                integration_bounds::NTuple{2, Real};
                spectral_index::Real=0.) where {T<:Real}
    # passband weighted mean calculation + spectral index usage
    freq_range, passband = get_restricted_regions(frequencies, passband,
                                                  integration_bounds)
    numerator = trapz(freq_range, freq_range.^(1. + spectral_index) .* passband)
    denom = trapz(freq_range, freq_range.^(0. + spectral_index) .* passband)
    return numerator / denom
end

function center(frequencies::Vector{T}, passband::Vector{T},
                integration_bounds::NTuple{2, Real},
                σ::Function) where {T<:Real}
    # passband weighted mean calculation + spectral index usage
    # ν is our frequency range
    ν, passband = get_restricted_regions(frequencies, passband,
                                         integration_bounds)
    numerator = trapz(ν, σ.(ν) .* passband .* (ν.^(-1.)))
    denom = trapz(ν, σ.(ν) .* passband .* ν.^(-2.))
    return numerator / denom
end

function width(frequencies::Vector{T}, passband::Vector{T},
               integration_bounds::NTuple{2, Real}) where {T<:Real}
    # Dicke bandwidth
    freq_range, passband = get_restricted_regions(frequencies, passband,
                                                  integration_bounds)
    numerator = trapz(freq_range, passband)^2
    denom = trapz(freq_range, passband.^2)
    return numerator / denom
end

function signal_width(frequencies::Vector{T}, passband::Vector{T},
                      integration_bounds::NTuple{2, Real}) where {T<:Real}
    freq_range, passband = get_restricted_regions(frequencies, passband,
                                                  integration_bounds)
    return trapz(freq_range, passband)
end

function get_band_edges_manually(
        passband::Vector{<:Real}; limit::Real=.05)
    # go until we hit .05
    current_index = argmax(passband)
    # get the upper edge
    while (abs(passband[current_index]) > limit)
        current_index += 1
    end
    end_index = current_index

    current_index = argmax(passband)
    # get the lower edge
    while (passband[current_index] > limit)
        current_index -= 1
    end
    start_index = current_index
    return (start_index, end_index)
end

function fit_band_edge(
        frequencies::Vector{<:Real}, passband::Vector{<:Real},
        edge_guess_ind::Real; limit::Real=.05, ind_limit::Int=4)

    @. cubic(x, p) = p[1] * x^3 + p[2] * x^2 + p[3]
    minimumby(iter, f) = reduce(iter) do x, y
        f(x) < f(y) ? x : y
    end

    fit_x = frequencies[edge_guess_ind - ind_limit : edge_guess_ind + ind_limit - 1]
    fit_y = passband[edge_guess_ind - ind_limit : edge_guess_ind + ind_limit - 1]
    res = curve_fit(cubic, fit_x, fit_y, [1., 1., 1.])
    a, b, c = coef(res)
    roots = real.(Polynomials.roots(Polynomials.Polynomial([c - limit, 0, b, a])))
    closest_root = minimumby(roots, (x-> abs(x - frequencies[edge_guess_ind])))
    if closest_root < 10 || closest_root > maximum(frequencies) || abs(
        frequencies[edge_guess_ind] - closest_root) > 4
            return frequencies[edge_guess_ind]
    end
    return closest_root
end

function get_band_edges(
        frequencies::Vector{<:Real}, passband::Vector{<:Real};
        limit::Real=.05)
    try
        start_ind, end_ind = get_band_edges_manually(passband, limit=limit)
        low_edge = fit_band_edge(frequencies, passband, start_ind)
        upper_edge = fit_band_edge(frequencies, passband, end_ind)
        return [low_edge, upper_edge]
    catch e
        return [nothing, nothing]
    end
end

function normalize(frequencies::Vector{T},
                   passband::Vector{T}; start_freq::Real=0,
                   end_freq::Real=maximum(frequencies)) where T<:Real
    start_freq_index, end_freq_index = find_subarray_limits(
        frequencies, start_freq, end_freq)
    return normalize(passband, start_index=start_freq_index,
                     end_index=end_freq_index)
end

function normalize(passband::Vector{T}; start_index::Integer=1,
                   end_index::Integer=length(passband)) where T<:Real
    passband / maximum(passband[start_index:end_index])
end

function normalize_bands(passband_data::PassbandArray; start_freq::Real=0,
                         end_freq::Real=maximum(passband_data.frequencies))
    return hcat([normalize(passband_data.frequencies, passband_data.passbands[
        :, i], start_freq=start_freq, end_freq=end_freq) for i in 1:size(
            passband_data.passbands, 2)]...)
end


function average_bands(passbands::Matrix{<:Real})
    return normalize(vec(mean(passbands, dims=2)))
end

function average_bands(passbands::Matrix{<:Real}, attrs::DataFrame)
    return normalize(vec(mean(passbands, weights(attrs.:SNR.^2), dims=2)))
end

function average_bands(frequencies::Vector{<:Real}, passbands::Matrix{<:Real},
                       attrs::DataFrame; start_freq=1.,
                       end_freq=maximum(frequencies))
    return normalize(frequencies, vec(mean(passbands, weights(attrs.:SNR.^2),
                                           dims=2)),
                     start_freq=start_freq, end_freq=end_freq)
end

function average_bands(array::PassbandArray{<:Real}, start_freq=1.,
                       end_freq=maximum(array.frequencies))
    return normalize(array.frequencies, vec(mean(array.passbands, weights(
        array.attrs.:SNR.^2), dims=2)), start_freq=start_freq, end_freq=end_freq)
end

function average_bands(array::PassbandArray{<:Real},
                       snr_weight_func::Function, start_freq=1.,
                       end_freq=maximum(array.frequencies))
    return normalize(array.frequencies, vec(mean(array.passbands, weights(
        snr_weight_func(array.attrs.:SNR)), dims=2)), start_freq=start_freq,
                     end_freq=end_freq)
end

function get_centers(frequencies::AbstractVector, passbands::AbstractMatrix,
                     int_bounds::NTuple{2, Real};
                     spectral_index::Real=0.)
    return vec([center(frequencies, collect(band),
                   int_bounds, spectral_index=spectral_index)
        for band in eachcol(passbands)])
end

function get_centers(passband_data::PassbandArray{<:Real},
                     int_bounds::NTuple{2, Real}; spectral_index::Real=0.)
    return get_centers(passband_data.frequencies, passband_data.passbands,
                       int_bounds, spectral_index=spectral_index)
end

function get_widths(frequencies::AbstractVector, passbands::AbstractMatrix,
                    int_bounds::NTuple{2, Real})
    return vec([width(frequencies, collect(band), int_bounds)
        for band in eachcol(passbands)])
end

function get_widths(passband_data::PassbandArray{<:Real},
                   int_bounds::NTuple{2, Real})
    return get_widths(passband_data.frequencies, passband_data.passbands,
                      int_bounds)
end

function get_edges(frequencies::AbstractVector, passbands::AbstractMatrix,
                   int_bounds::NTuple{2, Real})
    return hcat([get_band_edges(frequencies, collect(band), int_bounds)
        for band in eachcol(passbands)]...)
end

function get_edges(frequencies::AbstractVector, passbands::AbstractMatrix;
                   limit::Real=0.05)
    return hcat([get_band_edges(frequencies, collect(band), limit=limit)
        for band in eachcol(passbands)]...)
end

function get_edges(passband_data::PassbandArray{<:Real},
                   int_bounds::NTuple{2, Real})
    return get_edges(passband_data.frequencies, passband_data.passbands,
                     int_bounds)
end

function get_sigma_band_limits(
    bands::Matrix{T}, centers::Vector{T}) where T<:AbstractFloat
    # get the 68% confidence intervals of these bands!
    # get the sigmas of the mean
    lower_lim = mean(centers) - std(centers)
    upper_lim = mean(centers) + std(centers)
    # return indices which match these limits
    lower_lim_ind = argmin(map(abs, (centers .- lower_lim)))
    upper_lim_ind = argmin(map(abs, (centers .- upper_lim)))
    return bands[:, lower_lim_ind], bands[:, upper_lim_ind],
        lower_lim, upper_lim
end

function get_upper_lower_bands(
    bands::Matrix{<:AbstractFloat}, interval::AbstractFloat=.68)
    # sort the bands
    sorted_bands = sort(bands, dims=2)
    # get the 1 sigma confidence intervals
    upper_index = round(Integer, size(bands, 2) * interval)
    lower_index = max(round(Integer, size(bands, 2) * (1. - interval)), 1)
    return sorted_bands[:, upper_index], sorted_bands[:, lower_index]
end

function weighted_rms(
    data::Vector{T}, wts::Vector{S}) where {S<:Number, T<:Number}
    mean1 = mean(data, weights(wts))
    std1 = std(data, AnalyticWeights(wts), corrected=true)
    (mean1, std1)
end

weighted_rms(attrs::DataFrame, center_code="center") = weighted_rms(
    attrs[!, center_code], attrs.:SNR.^2)

weighted_rms(arr::PassbandArray, attr::String) = weighted_rms(
    arr.attrs[!, attr], arr.attrs.:SNR.^2)


function cut_cat!(passband_data::PassbandArray{<:Real},
                  cat::String; num_devs::AbstractFloat=4.)
    cat in names(passband_data.attrs) || error("cat must be in names of attrs!")
    cat_attrs = passband_data.attrs[!, cat]
    mean_cat, rms_cat = weighted_rms(cat_attrs, passband_data.attrs.:SNR.^2)
    vals_inside = findall(x->x <= mean_cat + num_devs * rms_cat
                          && x >= mean_cat - num_devs * rms_cat, cat_attrs)
    passband_data.passbands = passband_data.passbands[:, vals_inside]
    passband_data.attrs = passband_data.attrs[vals_inside, :]
    return

end

function cut_cat!(passband_data::PassbandArray{<:Real}, cat::String,
                  filter_func::Function)
    cat_attrs = passband_data.attrs[!, cat]
    good_inds = findall(filter_func, cat_attrs)
    passband_data.attrs = passband_data.attrs[good_inds, :]
    passband_data.passbands = passband_data.passbands[:, good_inds]
    return
end

function cut_cat(passband_data::PassbandArray{<:Real}, cat::String,
                 filter_func::Function)
    cat_attrs = passband_data.attrs[!, cat]
    good_inds = findall(filter_func, cat_attrs)
    # don't set the distances when cutting a category since they should
    # already be set!
    return PassbandArray(
        passband_data.passbands[:, good_inds], passband_data.frequencies,
        passband_data.attrs[good_inds, :], passband_data.label,
        set_distances=false)
end

function filter_data!(passband_data::PassbandArray{<:Real}, filter_func::Function)
    good_inds = [i for (i, row) in enumerate(eachrow(passband_data.attrs))
                     if filter_func(row)]
    #good_inds = findall(filter_func, passband_data.attrs)
    passband_data.attrs = passband_data.attrs[good_inds, :]
    passband_data.passbands = passband_data.passbands[:, good_inds]
end

function cut_bands!(passband_data::PassbandArray{<:Real},
                    num_devs::AbstractFloat=4.0; cut_half_snrs::Bool=false)
    # just cut for SNR > mean!
    mean_SNR = mean(passband_data.attrs.:SNR)
    valid_SNR = findall(x->x > mean_SNR, passband_data.attrs.:SNR)
    #new_data = copy(passband_data)
    if cut_half_snrs
        passband_data.passbands = passband_data.passbands[:, valid_SNR]
        passband_data.attrs = passband_data.attrs[valid_SNR, :]
    end

    for category in ["upper_edge", "lower_edge"]
        cut_cat!(passband_data, category, num_devs=5.)
    end

    # do the edges twice, and cut other categories with threshold
    for category in ["upper_edge", "lower_edge", "width", "center"]
        cut_cat!(passband_data, category, num_devs=num_devs)
    end
    return
end

function cut_repeat_measurements!(passband_data::PassbandArray{<:Real})
    # Cut all detectors with repeat measurements
    # For this we only keep the highest SNR detector measurement and discard
    # all ther rest.
    # Alternatively, we could average together these measurements with weight
    # SNR squared..
    dets = countmap(passband_data.attrs.:det)
    for (det, counts) in dets
        if counts >= 2
            truncated_array = cut_cat(passband_data, "det", x->x==det)
            highest_snr = maximum(truncated_array.attrs.:SNR)
            # get the run num of this SNR -- if multiple just take the 1st
            snr_ind = argmax(truncated_array.attrs.:SNR)
            truncated_array.attrs[snr_ind, "SNR"] ≈ highest_snr || error(
                "SNR should match!")
            snr_run_num = truncated_array.attrs.:run_num[snr_ind]
            filter_data!(passband_data, x->(
                x.:run_num == snr_run_num || x.:det != det))
            truncated_array = cut_cat(passband_data, "det", x->x==det)
            size(truncated_array.passbands, 2) == 1 || error(
                "only want 1 measurement!")
        end
    end

end

get_weights(passband_data::PassbandArray, key::String="SNR") =
    passband_data.attrs[!, key].^2

"""Placeholder for the different passband calculation functions we have.
(e.g. center over various source and spectral index, width)"""
struct BandCalculation
    func::Function
    args::Tuple
    kwargs::Dict
    label::String
end

function BandCalculation(func::Function, args::Tuple, label::String)
    return BandCalculation(func, args, Dict(), label)
end

function calculate_attr(frequencies::Vector{<:Real}, passband::Vector{<:Real},
                        calculation_func::BandCalculation)
    return calculation_func.func(
        frequencies, passband, calculation_func.args...;
        calculation_func.kwargs...)
end

function calculate_attrs(frequencies::Vector{<:Real}, passbands::Matrix{<:Real},
                         calculation_func::BandCalculation)
    return vec([calculate_attr(frequencies, collect(band), calculation_func)
        for band in eachcol(passbands)])
end

function weight_func(snr::Real, max_snr::Real=90)
    weight = snr
    if weight > max_snr
        weight = snr
    end
    return weight^2
end

weight_snrs(snrs::Vector{<:Real}) = weight_func.(snrs)
