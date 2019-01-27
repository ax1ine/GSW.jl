using GSW

    println(libgswteos)

    sa=35.
    sp=35.
    ct=20.
    t=20.
    pt=20.
    p=10.
    h=1.15103e5
    saturation_fraction=0.5

    println(gsw_rho(sa,ct,p));
    println(gsw_adiabatic_lapse_rate_from_ct(sa, ct, p));
    println(gsw_adiabatic_lapse_rate_ice(t, p));
    println(gsw_alpha(sa, ct, p));
    println(gsw_alpha_on_beta(sa, ct, p));
    println(gsw_alpha_wrt_t_exact(sa, t, p));
    println(gsw_alpha_wrt_t_ice(t, p));
    println(gsw_beta_const_t_exact(sa, t, p));
    println(gsw_beta(sa, ct, p));
    println(gsw_cabbeling(sa, ct, p));
    println(gsw_c_from_sp(sp, t, p));
    println(gsw_chem_potential_water_ice(t, p));
    println(gsw_chem_potential_water_t_exact(sa, t, p));
    println(gsw_cp_ice(t, p));
    println(gsw_cp_t_exact(sa, t, p));
