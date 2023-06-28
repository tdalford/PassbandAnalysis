import CSV

function final_band_release(
    passband_data::PassbandArray{<:Real}, z_sim_data::Vector{SimData},
    center_freq::Real, upper_corrections::AbstractVector,
    lower_corrections::AbstractVector, data_freqs::AbstractVector,
    int_bounds::NTuple{2, Real}; return_mean_arr::Bool=true)

    orig_bands = passband_data.passbands
    error_data = Dict()
    # get the mean corrections as (upper + lower) / 2
    mean_corrections = [(u .+ l) ./ 2 for (u, l) in zip(upper_corrections,
                                                        lower_corrections)]
    # have to invert and then average to make this work.
    # mean_corrections = [1 ./ ((1 ./ u) + (1 ./ l)) ./ 2 for (u, l) in zip(
    #                         upper_corrections, lower_corrections)]
    # just apply the systematic corrections to the average... oh well too lazy
    for (corrections, key) in zip([upper_corrections, lower_corrections],
                                  ["upper_systematic", "lower_systematic"])
        final_array = apply_corrections(
            passband_data, orig_bands, corrections, data_freqs)
        # get the average
        average = average_bands(final_array.passbands, final_array.attrs)
        error_data[key] = average
    end

    final_array = apply_corrections(
        passband_data, orig_bands, mean_corrections, data_freqs)
    # now do this for the mean corrections
    # apply for all our possible centers: beam filling and point sources,
    # CMB, width, etc!
    final_data = get_band_average_and_statistical_sample(final_array)
    rmses = get_rmses(final_array, z_sim_data, center_freq, int_bounds)
    error_data["mean_band"] = final_data["average"]
    error_data["statistical_sample"] = final_data["resampled_averages"]
    error_data["rms"] = maximum(rmses)
    if return_mean_arr
        error_data["mean_array"] = final_array
    end
    return error_data
end

function calculate_final_center_uncert(
    h5_data::Dict, calculation_func::BandCalculation, n_iters::Int=1000)
    values = []
    for _=1:n_iters
        freqs, band = monte_carlo_error_prop(
            h5_data["frequencies"], h5_data["mean_band"],
            h5_data["upper_systematic"], h5_data["lower_systematic"],
            h5_data["statistical_sample"],
            h5_data["frequency_shift_uncertainty"])
        value = calculate_attr(freqs, band, calculation_func)
        push!(values, value)
    end
    mean_center = calculate_attr(h5_data["frequencies"],
                                 normalize(h5_data["mean_band"]),
                                 calculation_func)
    return mean_center, std(values)
end

function monte_carlo_band_errs(h5_data::Dict, n_iters::Int=1000)
    bands = []
    for _=1:n_iters
        freqs, band = monte_carlo_error_prop(
            h5_data["frequencies"], h5_data["mean_band"],
            h5_data["upper_systematic"], h5_data["lower_systematic"],
            h5_data["statistical_sample"],
            h5_data["frequency_shift_uncertainty"])
        # interpolate the band
        interp_band = interpolate(
            freqs, band, h5_data["frequencies"])
        push!(bands, interp_band)
    end
    band_mat = hcat(bands...)
    mean_band = mean(band_mat, dims=2)
    std_band = std(band_mat, dims=2)
    # mean_band = mean_band ./ maximum(mean_band)
    # std_band = std_band ./ maximum(mean_band)
    upper_band = mean_band .+ std_band
    lower_band = mean_band .- std_band
    return h5_data["frequencies"], mean_band, upper_band, lower_band
end

"""Function for checking individual error sources"""
function indiv_error_prop(
    frequencies::Vector{T}, mean_band::Vector{T}, upper_sys::Vector{T},
    lower_sys::Vector{T}, statistical_sample::Matrix{T},
    frequency_shift_uncertainty::Real, prop_sys_error::Bool=true,
    prop_stat_error::Bool=true, prop_freq_shift::Bool=true) where {T<:Real}
    #total_lower = min.(lower_sys, upper_sys)
    #total_upper = max.(lower_sys, upper_sys)
    # sample the systematic distribution, get the difference from the mean band
    if prop_sys_error
        sys_band_sample = sample_systematic_band(upper_sys, lower_sys)
        sys_error = (sys_band_sample .- mean_band)
    else
        sys_error = zeros(length(mean_band))
    end
    # get a statistical sample, find the difference from the mean band
    if prop_stat_error
        stat_band_sample = sample_statistical_band(statistical_sample)
        stat_error = (stat_band_sample .- mean_band)
    else
        stat_error = zeros(length(mean_band))
    end
    # add these in quadrature?
    total_error = sys_error .+ stat_error
    error_band = total_error .+ mean_band
    # shift the frequencies
    if prop_freq_shift
        frequency_shift_factor = randn() * frequency_shift_uncertainty
        return (frequencies .- frequency_shift_factor), error_band
    else
        return frequencies, error_band
    end
end

function get_mean_array(final_band_release::Dict)
    upper = final_band_release["upper_systematic"]
    lower = final_band_release["lower_systematic"]
    mean_band = final_band_release["mean_band"]
    # get the difference between this and (upper + lower) / 2
    actual_mean = (upper .+ lower) ./ 2
    mean_correc = (actual_mean ./ mean_band)
    mean_band = normalize(mean_band .* mean_correc)

    # correct each of the passbands to the mean!
    mean_array = final_band_release["mean_array"]
    for i in 1:size(mean_array.passbands, 2)
        mean_array.passbands[:, i] = normalize(
            mean_array.passbands[:, i] .* mean_correc)
    end
    return mean_array
end

function single_temporal_variation(truncated_array::PassbandArray{<:Real})
    # assume the centers have already been calculated
    # get the std in center, width, SNR
    (center_mean, center_rms) = weighted_rms(truncated_array.attrs, "center")
    (width_mean, width_rms) = weighted_rms(truncated_array.attrs, "width")
    (snr_mean, snr_rms) = weighted_rms(truncated_array.attrs, "SNR")
    #snr_mean = mean(truncated_array.attrs.:SNR)
    #snr_std = std(truncated_array.attrs.:SNR)
    count = size(truncated_array.passbands, 2)
    return Dict(:center_mean=>center_mean, :center_rms=>center_rms,
                :width_mean=>width_mean, :width_rms=>width_rms,
                :snr_mean=>snr_mean, :snr_rms=>snr_rms, :count=>count)
end

function temporal_variation(mean_array::PassbandArray{<:Real},
                            count_threshold::Int=3)
    # check for differences in same index bands
    total_spreads = DataFrame(
        center_mean=Float64[], center_rms=Float64[],
        width_mean=Float64[], width_rms=Float64[],
        snr_mean=Float64[], snr_rms=Float64[], count=Int[])
    dets = countmap(mean_array.attrs.:det)
    for (det, counts) in dets
        if counts >= count_threshold
            truncated_array = cut_cat(mean_array, "det", x->x==det)
            variation_dict = single_temporal_variation(truncated_array)
            # we should weight by sqrt(counts) but also include SNR factor...
            push!(total_spreads, variation_dict)
        end
    end
    return total_spreads
end

function modify_band_errs!(band_errs::Dict)
    # update the mean band
    upper = band_errs["upper_systematic"]
    lower = band_errs["lower_systematic"]
    mean_band = band_errs["mean_band"]
    # get the different between this and (upper + lower) / 2
    actual_mean = (upper .+ lower) ./ 2
    mean_correc = (actual_mean ./ mean_band)
    mean_band = normalize(mean_band .* mean_correc)
    # modify the statistical sample
    sample = copy(band_errs["statistical_sample"])
    for i in 1:size(sample, 2)
        sample[:, i] = normalize(sample[:, i] .* mean_correc)
    end
    # Now re-set the mean band, statistical sample
    band_errs["mean_band"] = mean_band
    band_errs["statistical_sample"] = sample
end

function load_act_data()
    data_dict = Dict()
    labels = ["PA4_220", "PA4_150", "PA5_150", "PA5_90", "PA6_150", "PA6_90", "PA7_27",
            "PA7_37"]
    fbase = "/Users/talford/FTS-Simulation/raw_passband_data/"
    for label in labels
        if occursin("PA4", label)
        fname = fbase * "new_PA4/" * label * "Ghz.h5"
        else
            fname = fbase * label * "Ghz.h5"
        end
        # Make sure that we swap the X and Y axes here! Since the ACT
        # convention is backwards from the expected FTS direction.
        data = load_python_passband_data(fname, label=label, swap_xy=true)
        data_dict[label] = concat(data)
        # cut out same dets to only include 1 measurement
        cut_repeat_measurements!(data_dict[label])
        #cut_bands!(data_dict[label], 3.)
    end

    #fname = "/Users/tommyalford/FTS_simulation_results/shift_attr_data_diffraction.h5"
    fname = "/Users/talford/FTS_simulation_results/shift_attr_data.h5"
    z_sim_data = load_python_sim_data(fname);
    z_sim_data_restricted = restrict_z_data(z_sim_data);
    thin_membrane_raw_data = CSV.read(
        "/Users/talford/Downloads/corrections_unnormalized.dat",
        delim=",", DataFrame)
    thin_membrane_correction = thin_membrane_raw_data[:, 3]
    color_correction_freqs = thin_membrane_raw_data[:, 1]
    powers = .76 .^ (color_correction_freqs ./ 300.)
    lower_limit = powers .^ (2e-4 / 3e-4)
    upper_limit = powers .^ (4e-4 / 3e-4);

    correcs = Dict()
    for z in [-50, 0, 50]
        correcs[z] = get_corrections(z_sim_data_restricted(z + Z_FOCUS_MM))
    end

    freqs_defocus = z_sim_data_restricted(0).frequencies
    max_defocus = mean(correcs[0][7])
    min_defocus = mean(correcs[50][7]);
    return Dict("data_dict"=>data_dict,
                "thin_membrane_correction"=>thin_membrane_correction,
                "color_correction_freqs"=>color_correction_freqs,
                "lower_limit"=>lower_limit,
                "upper_limit"=>upper_limit,
                "freqs_defocus"=>freqs_defocus,
                "min_defocus"=>min_defocus,
                "max_defocus"=>max_defocus,
                "z_sim_data_restricted"=>z_sim_data_restricted)
end

function load_act_params()
    int_bounds = Dict("PA7_27"=>(20, 35), "PA7_37"=>(25, 53), "PA6_90"=>(70, 125),
                    "PA6_150"=>(115, 185), "PA5_90"=>(70, 125),
                    "PA5_150"=>(115, 185), "PA4_150"=>(115, 185), "PA4_220"=>(170, 300))
    center_freqs = Dict("PA7_27"=>27, "PA7_37"=>37, "PA6_90"=>90,
                    "PA6_150"=>150, "PA5_90"=>90, "PA5_150"=>150,
                    "PA4_150"=>150, "PA4_220"=>225)
    data_labels = Dict("PA7_27"=>"PA7 27 Ghz", "PA7_37"=>"PA7 37 Ghz", "PA6_90"=>"PA6 90 Ghz",
                    "PA6_150"=>"PA6 150 GHz", "PA5_90"=>"PA5 90 Ghz", "PA5_150"=>"PA5 150 Ghz",
                    "PA4_150"=>"PA4 150 Ghz", "PA4_220"=>"PA4 220 Ghz")
    freq_labels = Dict("PA7_27"=>"27 Ghz", "PA7_37"=>"37 Ghz", "PA6_90"=>"90 Ghz",
                    "PA6_150"=>"150 GHz", "PA5_90"=>"90 Ghz", "PA5_150"=>"150 Ghz",
                    "PA4_150"=>"150 Ghz", "PA4_220"=>"220 Ghz")
    colors = Dict("PA7_27"=>"red", "PA7_37"=>"orange", "PA6_90"=>"green",
                "PA6_150"=>"blue", "PA5_90"=>"green", "PA5_150"=>"blue",
                "PA4_150"=>"blue", "PA4_220"=>"purple")
    norm_start_freqs = Dict("PA7_27"=>10, "PA7_37"=>15, "PA6_90"=>50,
                            "PA6_150"=>60, "PA5_90"=>50, "PA5_150"=>60,
                            "PA4_150"=>60, "PA4_220"=>100)
    plot_start_freqs = Dict("PA7_27"=>10, "PA7_37"=>15, "PA6_90"=>60,
                            "PA6_150"=>80, "PA5_90"=>60, "PA5_150"=>80,
                            "PA4_150"=>80, "PA4_220"=>120)
    arr_names = Dict("PA7_27"=>"PA7_f030", "PA7_37"=>"PA7_f040",
                    "PA6_90"=>"PA6_f090", "PA6_150"=>"PA6_f150",
                    "PA5_90"=>"PA5_f090", "PA5_150"=>"PA5_f150",
                    "PA4_150"=>"PA4_f150", "PA4_220"=>"PA4_f220")
    linestyles = Dict("PA7_27"=>:dashdot, "PA7_37"=>:dashdot,
                    "PA6_90"=>:solid, "PA6_150"=>:solid,
                    "PA5_90"=>:dash, "PA5_150"=>:dash,
                    "PA4_150"=>:dot, "PA4_220"=>:dot)
    return Dict("int_bounds"=>int_bounds,
                "center_freqs"=>center_freqs,
                "data_labels"=>data_labels,
                "freq_labels"=>freq_labels,
                "colors"=>colors,
                "norm_start_freqs"=>norm_start_freqs,
                "plot_start_freqs"=>plot_start_freqs,
                "arr_names"=>arr_names,
                "linestyles"=>linestyles)
end

"""For segmenting corrections that are too are so large they cause
normalization issues. (Primilary the ACT thin membrane correction)."""
function segment_correction(freqs, correction, start_freq, end_freq)
    left_edge, right_edge = find_subarray_limits(freqs, convert(
        eltype(freqs), start_freq), convert(eltype(freqs), end_freq))
    left_correc, right_correc = correction[left_edge], correction[right_edge]
    new_correction = ones(length(correction)) * maximum([left_correc,
                                                         right_correc])
    new_correction[left_edge:right_edge] = correction[left_edge:right_edge]
    return normalize(new_correction)
end

"""Propagate our corrections and data into error properties"""
function get_total_band_errors(act_data::Dict, params::Dict)
    total_band_errors = Dict()
    for label in keys(params["freq_labels"])
        # segment the thin membrane correction since it's too strong
        thin_membrane = segment_correction(act_data[
            "color_correction_freqs"], act_data[
                "thin_membrane_correction"], params["norm_start_freqs"][
                    label], act_data["color_correction_freqs"][end])
        upper_corrections = [thin_membrane, act_data["upper_limit"], act_data[
            "min_defocus"]]
        lower_corrections = [thin_membrane, act_data["lower_limit"], act_data[
            "max_defocus"]]
        data_freqs = [act_data["color_correction_freqs"],
                      act_data["color_correction_freqs"],
                      act_data["freqs_defocus"]]
        band_errs = final_band_release(act_data["data_dict"][label],
                                       act_data["z_sim_data_restricted"],
                                       params["center_freqs"][label],
                                       upper_corrections, lower_corrections,
                                       data_freqs, params["int_bounds"][label])
        modify_band_errs!(band_errs)
        total_band_errors[label] = band_errs
    end
    return total_band_errors
end
