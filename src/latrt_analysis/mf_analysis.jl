function load_mf_data(;add_attrs::Vector=[])
    fbase = "/Users/talford/Downloads/"
    data_150 = concat(load_python_passband_data(
        fbase * "bands_150_20211129_uncorrected.h5", set_cal_factor=false));
    data_90 = concat(load_python_passband_data(
        fbase * "bands_90_20211129_uncorrected.h5", set_cal_factor=false));

    # load in the det map attrs
    # det_map_attrs = CSV.read(
    #     "/Users/talford/Downloads/mf_data_pixel_attrs.csv", DataFrame);

    # set_hard_coded_poly90!(data_90)
    # set_hard_coded_poly150!(data_150)


    # attrs_90 = filter(:msk_220 => ==(1), det_map_attrs)
    # attrs_150 = filter(:msk_280 => ==(1), det_map_attrs)
    # data_90.attrs = attrs_90
    # data_150.attrs = attrs_150;

    # new way of loading in the det map attrs with better saved map:
    det_map_attrs = CSV.read(
        "/Users/talford/Downloads/cv4_det_map_attrs.csv", DataFrame)
    add_det_map_to_passband_array!(data_90, det_map_attrs, add_attrs=add_attrs)
    add_det_map_to_passband_array!(data_150, det_map_attrs, add_attrs=add_attrs)

    cut_cat!(data_90, "pixel_centroid_y", x->x!=-1)
    cut_cat!(data_150, "pixel_centroid_y", x->x!=-1)
    # cut out any masks with don't match (from bias groups)
    cut_cat!(data_90, "msk_90", x->x==1)
    cut_cat!(data_150, "msk_150", x->x==1)

    # checked: the center is already corrected calculated as the RJ center!
    return data_90, data_150
end

function set_hard_coded_poly90!(data_90::PassbandArray)
    poly_90 = hcat([[-1.2, -1.2, 6, 10, 10, 6, -1.2],
                    [0, 10, 12, 10, 0, 2, 0]]...)
    poly_90 = [poly_90[i, :] for i in 1:size(poly_90, 1)];
    set_poly!(data_90, poly_90)
    cut_cat!(data_90, "inside", x->x==1)
end

function set_hard_coded_poly150!(data_150::PassbandArray)
    poly_150 = hcat([[-1, -1, 6, 6, -1], [0, 10, 12, 2, 0]]...)
    poly_150 = [poly_150[i, :] for i in 1:size(poly_150, 1)]
    set_poly!(data_150, poly_150)
    cut_cat!(data_150, "inside", x->x==1)
end
