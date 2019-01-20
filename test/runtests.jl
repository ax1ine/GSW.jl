using GSW
using NCDatasets

@show("Testing...")

sa=35.
sp=35.
ct=20.
t=20.
pt=20.
p=10.
h=1.15103e5
saturation_fraction=0.5

println(GSW.gsw_adiabatic_lapse_rate_from_ct(sa,ct,p));
println(GSW.gsw_rho(sa,ct,p));
println(GSW.gsw_adiabatic_lapse_rate_from_ct(sa, ct, p));
println(GSW.gsw_adiabatic_lapse_rate_ice(t, p));
println(GSW.gsw_alpha(sa, ct, p));
println(GSW.gsw_alpha_on_beta(sa, ct, p));
println(GSW.gsw_alpha_wrt_t_exact(sa, t, p));
println(GSW.gsw_alpha_wrt_t_ice(t, p));
println(GSW.gsw_beta_const_t_exact(sa, t, p));
println(GSW.gsw_beta(sa, ct, p));
println(GSW.gsw_cabbeling(sa, ct, p));
println(GSW.gsw_c_from_sp(sp, t, p));
println(GSW.gsw_chem_potential_water_ice(t, p));
println(GSW.gsw_chem_potential_water_t_exact(sa, t, p));
println(GSW.gsw_cp_ice(t, p));
println(GSW.gsw_cp_t_exact(sa, t, p));

"""
ct_sa = Ref{Cdouble}(0.0)
ct_pt = Ref{Cdouble}(0.0)
ct_sa[], ct_pt[] = GSW.gsw_ct_first_derivatives(sa, pt)
@show ct_sa[], ct_pt[]

ct_sa_wrt_t = Ref{Cdouble}(0.0)
ct_t_wrt_t  = Ref{Cdouble}(0.0)
ct_p_wrt_t  = Ref{Cdouble}(0.0)
ct_sa_wrt_t[], ct_t_wrt_t[], ct_p_wrt_t[] = GSW.gsw_ct_first_derivatives_wrt_t_exact(sa, t, p)
@show ct_sa_wrt_t[], ct_t_wrt_t[], ct_p_wrt_t[]
"""

println(GSW.gsw_ct_freezing(sa, t, saturation_fraction));
println(GSW.gsw_ct_freezing_poly(sa, p, saturation_fraction));
println(GSW.gsw_ct_from_enthalpy(sa, h, p));
println(GSW.gsw_ct_from_enthalpy_exact(sa, h, p));

ct = [20, 20]
t  = [20, 20]
println(GSW.gsw_sa_freezing_estimate(p, saturation_fraction, ct, t))
