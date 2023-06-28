using Plots

"""For plotting monte carlo error estimates of the ACT passbands"""
function plot_monte_carlo(band_errs::Dict, act_data::Dict, params::Dict;
                          freq_shift::Bool=false)
    #pgfplotsx()
    p = plot();
    for key in sort(collect(keys(band_errs)))
        # sample a bunch of bands
        samples = []
        for i=1:1000
            freqs, sample = monte_carlo_error_prop(
                act_data["data_dict"][key].frequencies, band_errs[key]["mean_band"],
                band_errs[key]["upper_systematic"],
                band_errs[key]["lower_systematic"],
                band_errs[key]["statistical_sample"], band_errs[key]["rms"])
            # interpolate the band
            if freq_shift
                interp_band = interpolate(
                    freqs, sample, act_data["data_dict"][key].frequencies)
            else
                interp_band = sample
            end
            push!(samples, interp_band)
        end

        # then get the mean and errorbars for these
        sample_mat = hcat(samples...)
        centers = get_centers(act_data["data_dict"][key].frequencies, sample_mat,
                              params["int_bounds"][key])
        #@show (key, mean(centers), std(centers))
        mean_band = mean(sample_mat, dims=2)
        std_band = std(sample_mat, dims=2)
        mean_band = mean_band ./ maximum(mean_band)
        std_band = std_band ./ maximum(mean_band)
        upper_band = mean_band .+ std_band
        lower_band = mean_band .- std_band

        start_index, _ = find_subarray_limits(
                act_data["data_dict"][key].frequencies, params["plot_start_freqs"][key],
            300)
        plot!(act_data["data_dict"][key].frequencies[start_index:end],
            mean_band[start_index:end], color=params["colors"][key],
                linestyle=params["linestyles"][key], label=params["arr_names"][key])
        plot!(act_data["data_dict"][key].frequencies[start_index:end],
            lower_band[start_index:end],
            fillrange=upper_band[start_index:end], color=params["colors"][key],
            alpha=.3, xlabel="Frequency (Ghz)")
    end
    return p
end

"""For plotting all the passband measurements we have together, along with the
mean."""
function plot_all_bands(band_errs::Dict, params::Dict, label::String;
                        renormalize_average::Bool=true, plot_kwargs...)
    p = plot()
    plot_all_bands!(band_errs, params, label;
                    renormalize_average=renormalize_average, plot_kwargs...)
    return p
end

function plot_all_bands!(band_errs::Dict, params::Dict, label::String;
                         renormalize_average::Bool=true, plot_kwargs...)
    mean_array = get_mean_array(band_errs)
    total_bands = normalize_bands(mean_array, start_freq=params[
        "norm_start_freqs"][label])
    if renormalize_average
        other_ave = vec(mean(total_bands, dims=2))
        mult_factor = maximum(other_ave)
    else
        mult_factor = 1
    end
    plot!(mean_array.frequencies, total_bands, color=params["colors"][label],
          xlim=params["int_bounds"][label], ylim=(-.1, 1.05);
          plot_kwargs...)
    # plot!(mean_array.frequencies, other_ave)
    plot!(mean_array.frequencies, band_errs["mean_band"] .* mult_factor,
          color=params["colors"][label], label=replace(label, "_"=>" "))
    return
end
