function load_UV31(;int_bounds::Dict=Dict("220" => (175, 280),
                                          "280" => (240, 330)),
                   fbase::String= "/Users/talford/Downloads/",
                   fname_220::String="band_data_220_20230418_UV31.h5",
                   fname_280::String="band_data_280_20230418_UV31.h5",
                   det_map_attr_fname::String=(
                       "/Users/talford/Downloads/uv31_det_map_attrs.csv"))
    data_220 = concat(load_python_passband_data(
        fbase * fname_220, set_cal_factor=false, load_attrs=["pol_group"]))
    data_280 = concat(load_python_passband_data(
        fbase * fname_280, set_cal_factor=false, load_attrs=["pol_group"]))

    poly = hcat([[-1.2, -1.2, 6, 10, 15, 15, -1.2], [
        -6, 10, 12, 12, 12, -6, -6]]...)
    poly = [poly[i, :] for i in 1:size(poly, 1)];
    det_map_attrs = CSV.read(det_map_attr_fname, DataFrame);

    for d in [data_220, data_280]
        convert_pol_group_to_ab!(d.attrs)
        cut_repeat_measurements!(d)
        add_det_map_to_passband_array!(d, det_map_attrs)
        set_poly!(d, poly)
        cut_cat!(d, "inside", x->x==1)
        cut_cat!(d, "pixel_centroid_y", x->x!=-1)
        # cut bands?
        cut_bands!(d, 6.)
    end

    data_220.attrs[!, "center_RJ"] = get_centers(
        data_220.frequencies, data_220.passbands, int_bounds["220"])
    data_280.attrs[!, "center_RJ"] = get_centers(
        data_280.frequencies, data_280.passbands, int_bounds["280"])

    edges_220 = get_edges(data_220.frequencies, data_220.passbands,
                          int_bounds["220"])
    edges_280 = get_edges(data_280.frequencies, data_280.passbands,
                          int_bounds["280"])
    data_220.attrs[!, "recompute_lower_edge"] = edges_220[1, :]
    data_220.attrs[!, "recompute_upper_edge"] = edges_220[2, :]
    data_280.attrs[!, "recompute_lower_edge"] = edges_280[1, :]
    data_280.attrs[!, "recompute_upper_edge"] = edges_280[2, :]

    data_220.attrs[!, "SNR_nocut"] =  1 ./ std.(
        eachcol(data_220.passbands[172:end, :]))
    data_280.attrs[!, "SNR_nocut"] =  1 ./ std.(
        eachcol(data_280.passbands[172:end, :]));

    # cut_cat!(data_220, "SNR_nocut", x->x>100)
    # cut_cat!(data_280, "SNR_nocut", x->x>100)
    return data_220, data_280
end

function load_UV8(;int_bounds::Dict=Dict("220" => (170, 290),
                                         "280" => (240, 330)),
                   fbase::String= "/Users/talford/Downloads/",
                   fname_220::String="band_data_220_20230308_UV8_with_beam_map.h5",
                   fname_280::String="band_data_280_20230308_UV8_with_beam_map.h5",
                   det_map_attr_fname::String=(
                       "/Users/talford/Downloads/uv8_det_map_attrs.csv"))
    data_220 = concat(load_python_passband_data(
        fbase * fname_220, set_cal_factor=false))
    data_280 = concat(load_python_passband_data(
        fbase * fname_280, set_cal_factor=false))
    poly = hcat([[-18, -18, -10, -2, -2, -2, -18], [
        -6, 10, 12, 12, 12, -6, -6]]...)
    poly = [poly[i, :] for i in 1:size(poly, 1)];
    det_map_attrs = CSV.read(det_map_attr_fname, DataFrame);

    for d in [data_220, data_280]
        # convert_pol_group_to_ab!(d.attrs)
        cut_repeat_measurements!(d)
        add_det_map_to_passband_array!(d, det_map_attrs)
        set_poly!(d, poly)
        cut_cat!(d, "inside", x->x==1)
        cut_cat!(d, "pixel_centroid_y", x->x!=-1)
    end

    cut_cat!(data_220, "SNR", x->x>20)
    cut_cat!(data_280, "SNR", x->x>20)

    data_220.attrs[!, "center_RJ"] = get_centers(
        data_220.frequencies, data_220.passbands, int_bounds["220"])
    data_280.attrs[!, "center_RJ"] = get_centers(
        data_280.frequencies, data_280.passbands, int_bounds["280"])

    edges_220 = get_edges(data_220.frequencies, data_220.passbands,
                          int_bounds["220"])
    edges_280 = get_edges(data_280.frequencies, data_280.passbands,
                          int_bounds["280"])
    data_220.attrs[!, "recompute_lower_edge"] = edges_220[1, :]
    data_220.attrs[!, "recompute_upper_edge"] = edges_220[2, :]
    data_280.attrs[!, "recompute_lower_edge"] = edges_280[1, :]
    data_280.attrs[!, "recompute_upper_edge"] = edges_280[2, :]

    data_220.attrs[!, "SNR_nocut"] =  1 ./ std.(
        eachcol(data_220.passbands[172:end, :]))
    data_280.attrs[!, "SNR_nocut"] =  1 ./ std.(
        eachcol(data_280.passbands[172:end, :]));
    return data_220, data_280
end
