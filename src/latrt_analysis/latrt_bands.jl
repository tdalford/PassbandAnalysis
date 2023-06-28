function add_det_map_to_passband_array!(arr::PassbandArray, det_map::DataFrame;
                                        add_attrs::Vector=["msk_220", "msk_280"])
    """
    x_loc and x_pos should match values, use these to match detectors!
    """
    main_attrs = ["pixel", "pixel_centroid_x", "pixel_centroid_y"]
    main_attrs = [main_attrs; add_attrs]
    for k in main_attrs
        arr.attrs[!, k] .= -1.0
    end

    num_not_found = 0
    for pixel in eachrow(det_map)
        matching_pixel = findall(x->x == pixel["x_pos"], arr.attrs.:x_loc)
        if size(matching_pixel, 1) == 0
            num_not_found += 1
            continue
        end
        size(matching_pixel, 1) == 1 || error("size != 1")
        for k in main_attrs
            arr.attrs[matching_pixel, k] = [pixel[k]]
        end
    end
    num_not_found = sum(arr.attrs.:pixel .== -1)
    println("total not found = $(num_not_found)")
end

function set_poly!(data::PassbandArray, poly)
    data.attrs[!, "dist"] = hypot.(data.attrs.:x_loc, data.attrs.:y_loc)
    data.attrs[!, "inside"] = map(v->inpolygon(
        [v.:x_loc, v.:y_loc], poly), eachrow(data.attrs))
    return
end

function ndf_correction(frequency_ghz::Vector{<:Real}, a::Real, b::Real,
                        Rv::Real=.18)
    correction = (1 - Rv) * exp.(-1 * a .* frequency_ghz.^b)
    return correction
end

function ndf_correct_band(band::Vector{<:Real}, ndf_correction::Vector{<:Real},
                          end_ind::Int)
    return normalize((band ./ ndf_correction), end_index=end_ind)
end

function ndf_correct_band(frequencies::Vector{<:Real}, band::Vector{<:Real},
                          ndf_correction::Vector{<:Real};
                          end_freq::Real=maximum(frequencies))
    return normalize(frequencies, (band ./ ndf_correction), end_freq=end_freq)
end

function segment_positions(passband_data::PassbandArray{<:Real},
                           frac::Real=.33)
    center_inds = sortperm(passband_data.attrs.:dist)[1:trunc(Integer, size(
        passband_data.attrs, 1) * frac)]
    center_ave = average_bands(passband_data.passbands[:, center_inds])
    outside_inds = sortperm(
        passband_data.attrs.:dist, order=Base.Order.Reverse)[1:trunc(
            Integer, size(passband_data.attrs, 1) * frac)]
    outside_ave = average_bands(passband_data.passbands[:, outside_inds])
    return center_ave, mean(passband_data.attrs.:center[center_inds]),
        outside_ave, mean(passband_data.attrs.:center[outside_inds])
end

function ndf_correct_bands!(arr::PassbandArray, ndf_correction::Vector{<:Real},
                            int_bounds::NTuple{2, Real};
                            end_freq::Real=maximum(arr.frequencies))
    ndf_corrected_bands = [ndf_correct_band(arr.frequencies, collect(
        band), ndf_correction, end_freq=end_freq) for band in eachcol(arr.passbands)]
    arr.passbands = hcat(ndf_corrected_bands...)
    centers = get_centers(arr, int_bounds)
    widths = get_widths(arr, int_bounds)
    edges = get_edges(arr, int_bounds)
    arr.attrs[!, "NDF_corrected_center"] .= centers
    arr.attrs[!, "NDF_corrected_width"] .= widths
    arr.attrs[!, "NDF_corrected_lower_edge"] .= edges[1, :]
    arr.attrs[!, "NDF_corrected_upper_edge"] .= edges[2, :]
    return
end

function ndf_correct_bands(arr::PassbandArray, ndf_correction::Vector{<:Real},
                           int_bounds::NTuple{2, Real};
                           end_freq::Real=maximum(arr.frequencies),
                           offset::Real=0)
    arr = copy(arr)
    arr.passbands .+= offset
    ndf_correct_bands!(arr, ndf_correction, int_bounds, end_freq=end_freq)
    return arr
end

# Functions for computing out of band leakage!
function get_upper_lower_values(values::Vector{<:Real},
                                interval::AbstractFloat=.68)
    # sort the values
    sorted_values = sort(values)
    # get the 1 sigma confidence intervals
    upper_index = round(Integer, size(values, 1) * interval)
    lower_index = max(round(Integer, size(values, 1) * (1. - interval)), 1)
    return sorted_values[upper_index], sorted_values[lower_index]
end

function get_mean_std(values::Vector{<:Real})
    upper, lower = get_upper_lower_values(values, .68)
    mean = (upper + lower) / 2
    std = (upper - lower) / 2
    return mean, std
end

function get_power_attrs(frequencies::Vector{<:Real}, passband::Vector{<:Real},
                         region_bounds::NTuple{2, Real})
    return signal_width(frequencies, passband, region_bounds)
end

function get_power_bounds(frequencies::Vector{<:Real}, band_samples::Matrix{<:Real},
                          region::NTuple{2, Real})
    rounded_region = round.(region, digits=0)
    powers = []
    for passband in eachcol(band_samples)
        offset_passband = collect(passband)# + -1 * minimum(passband[30:])
        power = get_power_attrs(frequencies, offset_passband, rounded_region)
        push!(powers, power)
    end
    return identity.(vec(powers))
end

function get_region_percent(frequencies::Vector{<:Real}, band_samples::Matrix{
        <:Real}, region::NTuple{2, Real}, powers_band::Vector{<:Real})
    powers = get_power_bounds(frequencies, band_samples, region)
    ratios = 100 * (powers ./ powers_band)
    return get_mean_std(ratios)
end

function get_confidence_intervals(data::Vector{<:Real}, interval::Real=.99)
    # now just change so that our percent limits are the +/- interval
    # ranges then convert to percent!
    mean_val = mean(data)
    sorted = sort(data)
    upper_index = round(Integer, size(sorted, 1) * interval)
    lower_index = max(round(Integer, size(sorted, 1) * (1 - interval)), 1)
    # recompute mean to be perfectly between
    mean_val = (sorted[upper_index] + sorted[lower_index]) / 2
    return mean_val, sorted[lower_index], sorted[upper_index]
end

function get_confidence_intervals(data::Vector{<:Real}, weights::Vector{<:Real},
                                  interval::Real=.99)
    mean_val, std_val = weighted_rms(data, weights)
    alpha = 1 - interval # significance level
    z = norm.ppf(1-alpha/2) # inverse CDF of standard normal distribution
    lower = mean - z*std_val
    upper = mean + z*std_val
    return mean_val, lower, upper
end

function to_pct_diff(mean::Real, lower::Real, upper::Real)
    lower_pct = 100(lower - mean) / mean
    upper_pct = 100(upper - mean) / mean
    return mean, lower_pct, upper_pct
end
