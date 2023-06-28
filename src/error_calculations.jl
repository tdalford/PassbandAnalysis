function bootstrap_random_error(
    passband_data::PassbandArray{<:Real}, n_iters::Int=30)
    # size is (length(frequencies)) *
    averages = Vector{Float64}[]
    orig_attrs = passband_data.attrs
    orig_bands = passband_data.passbands
    for _ in 1:n_iters
        # now sample a subset of these with replacement if specified!
        indices = rand(1:size(orig_attrs, 1), size(orig_attrs, 1))
        jitter_attrs = orig_attrs[indices, :]
        jitter_bands = orig_bands[:, indices]
        # copy the passband data into a new array- don't set the distances!!!
        average_band = average_bands(jitter_bands, jitter_attrs)
        push!(averages, average_band)
    end
    return hcat(averages...)
end

function get_band_average_and_statistical_sample(
    passband_data::PassbandArray{<:Real}, n_iters::Int=100)
    average_band = average_bands(passband_data.passbands, passband_data.attrs)
    random_bands = bootstrap_random_error(passband_data, n_iters)
    data = Dict("average"=>average_band,
                "resampled_averages"=>random_bands)
    return data
end

function get_band_mean_and_errs(
    passband_data::PassbandArray{<:Real}, calculation_func::BandCalculation)
    average_band = average_bands(passband_data.passbands, passband_data.attrs)
    mean = calculation_func.func(
        passband_data.frequencies, average_band, calculation_func.args...;
        calculation_func.kwargs...)
    random_bands = bootstrap_random_error(passband_data)
    random_centers = calculate_attrs(
        passband_data.frequencies, random_bands, calculation_func)
    rand_error = std(random_centers)
    upper_random, lower_random = get_upper_lower_bands(random_bands)
    data = Dict("average_band"=>average_band, "upper_random"=>upper_random,
                "lower_random"=>lower_random,
                "rand_error"=>rand_error, "mean"=>mean, "arr"=>passband_data)
    return data
end


function sample_statistical_band(statistical_sample::Matrix{<:Real})
    # second dimension of this should be the band!
    stat_num = rand(1:size(statistical_sample, 2))
    stat_band = statistical_sample[:, stat_num] |> vec
    return stat_band
end

function sample_systematic_band(upper_sys::Vector{T},
                                lower_sys::Vector{T}) where {T<:Real}
    sys_num = rand()
    sys_sample = (lower_sys .+ sys_num .* (upper_sys .- lower_sys))
    return sys_sample
end

function monte_carlo_error_prop(
    frequencies::Vector{T}, mean_band::Vector{T}, upper_sys::Vector{T},
    lower_sys::Vector{T}, statistical_sample::Matrix{T},
    frequency_shift_uncertainty::Real) where {T<:Real}
    #total_lower = min.(lower_sys, upper_sys)
    #total_upper = max.(lower_sys, upper_sys)
    # sample the systematic distribution, get the difference from the mean band
    sys_band_sample = sample_systematic_band(upper_sys, lower_sys)
    sys_error = (sys_band_sample .- mean_band)
    # get a statistical sample, find the difference from the mean band
    stat_band_sample = sample_statistical_band(statistical_sample)
    stat_error = (stat_band_sample .- mean_band)
    # add these in quadrature?
    total_error = sys_error .+ stat_error
    error_band = total_error .+ mean_band
    # shift the frequencies
    frequency_shift_factor = randn() * frequency_shift_uncertainty
    return (frequencies .- frequency_shift_factor), error_band
end
