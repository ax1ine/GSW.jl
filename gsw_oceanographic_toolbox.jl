"""
Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10
Documentation: http://www.teos-10.org/pubs/gsw/html/gsw_contents.html
These declarations facilitate the use of TEOS-10 functions with Julia 1.0
"""

#path to precompiled teos-10 library (x64)
if Sys.islinux()
  const libgswteos = joinpath(@__DIR__, "gsw/libgswteos-10.so")
end
if Sys.iswindows()
  const libgswteos = joinpath(@__DIR__, "gsw/libgswteos-10.dll")
end


function gsw_add_barrier(input_data, lon, lat, long_grid, lat_grid, dlong_grid, dlat_grid)
  ccall(("gsw_add_barrier",libgswteos),Cvoid,(Ptr{Cdouble},Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),input_data,lon,lat,long_grid,lat_grid,dlong_grid,dlat_grid,output_data)
  return output_data[]
end

# Replaces NaN's with non-nan mean of the 4 adjacent neighbours
#extern void   gsw_add_mean(double *data_in, double *data_out);
function gsw_add_mean(data_in)
  ccall(("gsw_add_mean",libgswteos),Cvoid,(Ptr{Cdouble},Ptr{Cdouble}),data_in,data_out)
  return data_out[]
end

# Adiabatic lapse rate from Conservative Temperature
#extern double gsw_adiabatic_lapse_rate_from_ct(double sa, double ct, double p);
function gsw_adiabatic_lapse_rate_from_ct(sa, ct, p)::Cdouble
  return ccall(("gsw_adiabatic_lapse_rate_from_ct",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Adiabatic lapse rate of ice
#extern double gsw_adiabatic_lapse_rate_ice(double t, double p);
function gsw_adiabatic_lapse_rate_ice(t, p)::Cdouble
  return ccall(("gsw_adiabatic_lapse_rate_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Thermal expansion coefficient with respect to CT
#extern double gsw_alpha(double sa, double ct, double p);
function gsw_alpha(sa, ct, p)::Cdouble
  return ccall(("gsw_alpha",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Alpha divided by beta
#extern double gsw_alpha_on_beta(double sa, double ct, double p);
function gsw_alpha_on_beta(sa, ct, p)::Cdouble
  return ccall(("gsw_alpha_on_beta",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Thermal expansion coefficient with respect to in-situ temperature
#extern double gsw_alpha_wrt_t_exact(double sa, double t, double p);
function gsw_alpha_wrt_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_alpha_wrt_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Thermal expansion coefficient of ice with respect to in-situ temperature
#extern double gsw_alpha_wrt_t_ice(double t, double p);
function gsw_alpha_wrt_t_ice(t, p)::Cdouble
  return ccall(("gsw_alpha_wrt_t_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Saline contraction coefficient at constant in-situ temperature
#extern double gsw_beta_const_t_exact(double sa, double t, double p);
function gsw_beta_const_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_beta_const_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Saline contraction coefficient at constant CT
#extern double gsw_beta(double sa, double ct, double p);
function gsw_beta(sa, ct, p)::Cdouble
  return ccall(("gsw_beta",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Cabbeling coefficient (75-term equation)
#extern double gsw_cabbeling(double sa, double ct, double p);
function gsw_cabbeling(sa, ct, p)::Cdouble
  return ccall(("gsw_cabbeling",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Conductivity from Practical Salinity (incl. for SP < 2)
#extern double gsw_c_from_sp(double sp, double t, double p);
function gsw_c_from_sp(sp, t, p)::Cdouble
  return ccall(("gsw_c_from_sp",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sp,t,p)
end

# Chemical potential of water in ice
#extern double gsw_chem_potential_water_ice(double t, double p);
function gsw_chem_potential_water_ice(t, p)::Cdouble
  return ccall(("gsw_chem_potential_water_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Chemical potential of water in seawater
#extern double gsw_chem_potential_water_t_exact(double sa, double t, double p);
function gsw_chem_potential_water_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_chem_potential_water_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Isobaric heat capacity of ice
#extern double gsw_cp_ice(double t, double p);
function gsw_cp_ice(t, p)::Cdouble
  return ccall(("gsw_cp_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Isobaric heat capacity
#extern double gsw_cp_t_exact(double sa, double t, double p);
function gsw_cp_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_cp_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# First derivatives of Conservative Temperature
#extern void   gsw_ct_first_derivatives (double sa, double pt, double *ct_sa, double *ct_pt);
function gsw_ct_first_derivatives(sa, pt)
  ccall(("gsw_ct_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,pt,ct_sa,ct_pt)
  return ct_sa[], ct_pt[]
end

# First derivatives of Conservative Temperature with respect to t
#extern void   gsw_ct_first_derivatives_wrt_t_exact(double sa, double t, double p, double *ct_sa_wrt_t, double *ct_t_wrt_t, double *ct_p_wrt_t);
function gsw_ct_first_derivatives_wrt_t_exact(sa, t, p)
  ccall(("gsw_ct_first_derivatives_wrt_t_exact",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,t,p,ct_sa_wrt_t,ct_t_wrt_t,ct_p_wrt_t)
  return ct_sa_wrt_t[], ct_t_wrt_t[], ct_p_wrt_t[]
end

# Conservative Temperature freezing point
#extern double gsw_ct_freezing(double sa, double p, double saturation_fraction);
function gsw_ct_freezing(sa, t, saturation_fraction)::Cdouble
  return ccall(("gsw_ct_freezing",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,saturation_fraction)
end

# First derivatives of Conservative Temperature freezing temperature of seawater
#extern void   gsw_ct_freezing_first_derivatives(double sa, double p, double saturation_fraction, double *ctfreezing_sa, double *ctfreezing_p);
function gsw_ct_freezing_first_derivatives(sa, p, saturation_fraction)
  ccall(("gsw_ct_freezing_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,p,saturation_fraction,ctfreezing_sa,ctfreezing_p)
  return ctfreezing_sa[],ctfreezing_p[]
end

# First derivatives of Conservative Temperature freezing temperature of seawater (polynomial)
#extern void   gsw_ct_freezing_first_derivatives_poly(double sa, double p, double saturation_fraction, double *ctfreezing_sa, double *ctfreezing_p);
function gsw_ct_freezing_first_derivatives_poly(sa, p, saturation_fraction)
  ccall(("gsw_ct_freezing_first_derivatives_poly",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,p,saturation_fraction,ctfreezing_sa,ctfreezing_p)
  return ctfreezing_sa[],ctfreezing_p[]
end

# Conservative Temperature freezing point (poly)
#extern double gsw_ct_freezing_poly(double sa, double p, double saturation_fraction);
function gsw_ct_freezing_poly(sa, p, saturation_fraction)::Cdouble
  return ccall(("gsw_ct_freezing_poly",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,p,saturation_fraction)
end

# Conservative Temperature from specific enthalpy of seawater (75-term equation)
#extern double gsw_ct_from_enthalpy(double sa, double h, double p);
function gsw_ct_from_enthalpy(sa, h, p)::Cdouble
  return ccall(("gsw_ct_from_enthalpy",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,h,p)
end

# Conservative Temperature as a function of enthalpy
#extern double gsw_ct_from_enthalpy_exact(double sa, double h, double p);
function gsw_ct_from_enthalpy_exact(sa, h, p)::Cdouble
  return ccall(("gsw_ct_from_enthalpy_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,h,p)
end

#extern double gsw_ct_from_entropy(double sa, double entropy);
function gsw_ct_from_entropy(sa, entropy)::Cdouble
  return ccall(("gsw_ct_from_entropy",libgswteos),Cdouble,(Cdouble,Cdouble),sa,entropy)
end

# Conservative Temperature from potential temperature
#extern double gsw_ct_from_pt(double sa, double pt);
function gsw_ct_from_pt(sa, pt)::Cdouble
  return ccall(("gsw_ct_from_pt",libgswteos),Cdouble,(Cdouble,Cdouble),sa,pt)
end

# Conservative Temperature from density
#extern void   gsw_ct_from_rho(double rho, double sa, double p, double *ct, double *ct_multiple);
function gsw_ct_from_rho(rho, sa, p)
  ccall(("gsw_ct_from_rho",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),rho,sa,p,ct,ct_multiple)
  return ct[],ct_multiple[]
end

# Conservative Temperature from in-situ temperature
#extern double gsw_ct_from_t(double sa, double t, double p);
function gsw_ct_from_t(sa, t, p)::Cdouble
  return ccall(("gsw_ct_from_t",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Conservative Temperature of maximum density of seawater
#extern double gsw_ct_maxdensity(double sa, double p);
function gsw_ct_maxdensity(sa, p)::Cdouble
  return ccall(("gsw_ct_maxdensity",libgswteos),Cdouble,(Cdouble,Cdouble),sa,p)
end

#extern void   gsw_ct_second_derivatives(double sa, double pt, double *ct_sa_sa, double *ct_sa_pt, double *ct_pt_pt);
function gsw_ct_second_derivatives(sa, pt)
  ccall(("gsw_ct_second_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,pt,ct_sa_sa,ct_sa_pt,ct_pt_pt)
  return ct_sa_sa[],ct_sa_pt[],ct_pt_pt[]
end

# Absolute Salinity Anomaly atlas value (excluding the Baltic Sea)
#extern double gsw_deltasa_atlas(double p, double lon, double lat);
function gsw_deltasa_atlas(p, lon, lat)::Cdouble
  return ccall(("gsw_deltasa_atlas",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),p,lon,lat)
end

# Absolute Salinity Anomaly from Practical Salinity
#extern double gsw_deltasa_from_sp(double sp, double p, double lon, double lat);
function gsw_deltasa_from_sp(sp, p, lon, lat)::Cdouble
  return ccall(("gsw_deltasa_from_sp",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sp,p,lon,lat)
end

# Dilution coefficient of seawater
#extern double gsw_dilution_coefficient_t_exact(double sa, double t, double p);
function gsw_dilution_coefficient_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_dilution_coefficient_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Dynamic enthalpy
#extern double gsw_dynamic_enthalpy(double sa, double ct, double p);
function gsw_dynamic_enthalpy(sa, ct, p)::Cdouble
  return ccall(("gsw_dynamic_enthalpy",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Enthalpy
#extern double gsw_enthalpy_ct_exact(double sa, double ct, double p);
function gsw_enthalpy_ct_exact(sa, ct, p)::Cdouble
  return ccall(("gsw_enthalpy_ct_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Difference of enthalpy between two pressures
#extern double gsw_enthalpy_diff(double sa, double ct, double p_shallow, double p_deep);
function gsw_enthalpy_diff(sa, ct, p, p_shallow, p_deep)::Cdouble
  return ccall(("gsw_enthalpy_diff",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble,Cdouble),sa,ct,p,p_shallow,p_deep)
end

# Enthalpy
#extern double gsw_enthalpy(double sa, double ct, double p);
function gsw_enthalpy(sa, ct, p)::Cdouble
  return ccall(("gsw_enthalpy",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# First derivatives of enthalpy
#extern void   gsw_enthalpy_first_derivatives_ct_exact(double sa, double ct, double p, double *h_sa, double *h_ct);
function gsw_enthalpy_first_derivatives_ct_exact(sa, ct, p)
  ccall(("gsw_enthalpy_first_derivatives_ct_exact",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,h_sa,h_ct)
  return h_sa[],h_ct[]
end

# First derivatives of enthalpy
#extern void   gsw_enthalpy_first_derivatives(double sa, double ct, double p, double *h_sa, double *h_ct);
function gsw_enthalpy_first_derivatives(sa, ct, p)
  ccall(("gsw_enthalpy_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,h_sa,h_ct)
  return h_sa[],h_ct[]
end

# Enthalpy of ice
#extern double gsw_enthalpy_ice(double t, double p);
function gsw_enthalpy_ice(t, p)::Cdouble
  return ccall(("gsw_enthalpy_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Second derivatives of enthalpy
#extern void   gsw_enthalpy_second_derivatives_ct_exact(double sa, double ct, double p, double *h_sa_sa, double *h_sa_ct, double *h_ct_ct);
function gsw_enthalpy_second_derivatives_ct_exact(sa, ct, p)
  ccall(("gsw_enthalpy_second_derivatives_ct_exact",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,h_sa_sa,h_sa_ct,h_ct_ct)
  return h_sa_sa[],h_sa_ct[],h_ct_ct[]
end

#extern void   gsw_enthalpy_second_derivatives(double sa, double ct, double p, double *h_sa_sa, double *h_sa_ct, double *h_ct_ct);
function gsw_enthalpy_second_derivatives(sa, ct, p)
  ccall(("gsw_enthalpy_second_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,h_sa_sa,h_sa_ct,h_ct_ct)
  return h_sa_sa[],h_sa_ct[],h_ct_ct[]
end

# Enthalpy(35.16504,0,p)
#extern double gsw_enthalpy_sso_0(double p);
function gsw_enthalpy_sso_0(p)::Cdouble
  return ccall(("gsw_enthalpy_sso_0",libgswteos),Cdouble,(Cdouble,),p)
end

# Enthalpy
#extern double gsw_enthalpy_t_exact(double sa, double t, double p);
function gsw_enthalpy_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_enthalpy_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

#
#extern void   gsw_entropy_first_derivatives(double sa, double ct, double *eta_sa, double *eta_ct);
function gsw_entropy_first_derivatives(sa, ct)
  ccall(("gsw_entropy_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,eta_sa,eta_ct)
  return eta_sa[],eta_ct[]
end

# Entropy from Conservative Temperature
#extern double gsw_entropy_from_ct(double sa, double ct);
function gsw_entropy_from_ct(sa, ct)::Cdouble
  return ccall(("gsw_entropy_from_ct",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# Entropy from potential temperature
#extern double gsw_entropy_from_pt(double sa, double pt);
function gsw_entropy_from_pt(sa, pt)::Cdouble
  return ccall(("gsw_entropy_from_pt",libgswteos),Cdouble,(Cdouble,Cdouble),sa,pt)
end

# Entropy from in-situ temperature
#extern double gsw_entropy_from_t(double sa, double t, double p);
function gsw_entropy_from_t(sa, t, p)::Cdouble
  return ccall(("gsw_entropy_from_t",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Entropy of ice
#extern double gsw_entropy_ice(double t, double p);
function gsw_entropy_ice(t, p)::Cdouble
  return ccall(("gsw_entropy_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Entropy minus the terms that are a function of only SA
#extern double gsw_entropy_part(double sa, double t, double p);
function gsw_entropy_part(sa, t, p)::Cdouble
  return ccall(("gsw_entropy_part",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Entropy_part evaluated at 0 dbar
#extern double gsw_entropy_part_zerop(double sa, double pt0);
function gsw_entropy_part_zerop(sa, pt0)::Cdouble
  return ccall(("gsw_entropy_part_zerop",libgswteos),Cdouble,(Cdouble,Cdouble),sa,pt0)
end

# Second derivatives of entropy
#extern void   gsw_entropy_second_derivatives(double sa, double ct, double *eta_sa_sa, double *eta_sa_ct, double *eta_ct_ct);
function gsw_entropy_second_derivatives(sa, ct)
  ccall(("gsw_entropy_second_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,eta_sa_sa,eta_sa_ct,eta_ct_ct)
  return eta_sa_sa[],eta_sa_ct[],eta_ct_ct[]
end

# Ratio of Absolute to Preformed Salinity, minus 1
#extern double gsw_fdelta(double p, double lon, double lat);
function gsw_fdelta(p, lon, lat)::Cdouble
  return ccall(("gsw_fdelta",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),p,lon,lat)
end

# SA, CT and ice mass fraction from bulk SA and bulk enthalpy
#extern void   gsw_frazil_properties(double sa_bulk, double h_bulk, double p, double *sa_final, double *ct_final, double *w_ih_final);
function gsw_frazil_properties(sa_bulk, n_bulk, p)
  ccall(("gsw_frazil_properties",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa_bulk,n_bulk,p,sa_final,ct_final,w_ih_final)
  return sa_final[],ct_final[],w_ih_final[]
end

# SA, CT and ice mass fraction from bulk SA and bulk potential enthalpy
#extern void   gsw_frazil_properties_potential(double sa_bulk, double h_pot_bulk, double p, double *sa_final, double *ct_final, double *w_ih_final);
function gsw_frazil_properties_potential(sa_bulk, h_pot_bulk, p)
  ccall(("gsw_frazil_properties_potential",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa_bulk,h_pot_bulk,p,sa_final,ct_final,w_ih_final)
  return sa_final[],ct_final[],w_ih_final[]
end

# SA, CT and ice mass fraction from bulk SA and bulk potential enthalpy (polynomial)
#extern void   gsw_frazil_properties_potential_poly(double sa_bulk, double h_pot_bulk, double p, double *sa_final, double *ct_final, double *w_ih_final);
function gsw_frazil_properties_potential_poly(sa_bulk, h_pot_bulk, p)
  ccall(("gsw_frazil_properties_potential_poly",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa_bulk,h_pot_bulk,p,sa_final,ct_final,w_ih_final)
  return sa_final[],ct_final[],w_ih_final[]
end

# Ratios of SA, CT and P changes during frazil ice formation
#extern void   gsw_frazil_ratios_adiabatic(double sa, double p, double w_ih, double *dsa_dct_frazil, double *dsa_dp_frazil, double *dct_dp_frazil);
function gsw_frazil_ratios_adiabatic(sa, p, w_ih)
  ccall(("gsw_frazil_ratios_adiabatic",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,p,w_ih,dsa_dct_frazil,dsa_dp_frazil,dct_dp_frazil)
  return dsa_dct_frazil[],dsa_dp_frazil[],dct_dp_frazil[]
end

# Ratios of SA, CT and P changes during frazil ice formation (polynomial)
#extern void   gsw_frazil_ratios_adiabatic_poly(double sa, double p, double w_ih, double *dsa_dct_frazil, double *dsa_dp_frazil, double *dct_dp_frazil);
function gsw_frazil_ratios_adiabatic_poly(sa, p, w_ih)
  ccall(("gsw_frazil_ratios_adiabatic_poly",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,p,w_ih,dsa_dct_frazil,dsa_dp_frazil,dct_dp_frazil)
  return dsa_dct_frazil[],dsa_dp_frazil[],dct_dp_frazil[]
end

#
#extern double *gsw_geo_strf_dyn_height(double *sa, double *ct, double *p, double p_ref, int n_levels, double *dyn_height);
function gsw_geo_strf_dyn_height(sa, ct, p, p_ref, n_levels, dyn_height)::Cdouble
  return ccall(("gsw_geo_strf_dyn_height",libgswteos),Cdouble,(Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Cdouble,Int,Ptr{Cdouble}),sa, ct, p, p_ref, n_levels, dyn_height)
end

#extern int gsw_geo_strf_dyn_height_1(double *sa, double *ct, double *p, double p_ref, int n_levels, double *dyn_height, double max_dp_i, int interp_method);
function gsw_geo_strf_dyn_height_1(sa, ct, p, p_ref, n_levels, dyn_height, max_dp_i, interp_method)::Cdouble
  return ccall(("gsw_geo_strf_dyn_height_1",libgswteos),Int,(Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Cdouble,Int,Ptr{Cdouble},Cdouble, Int),sa, ct, p, p_ref, n_levels, dyn_height, max_dp_i, interp_method)
end

#
#extern double *gsw_geo_strf_dyn_height_pc(double *sa, double *ct,double *delta_p, int n_levels, double *geo_strf_dyn_height_pc,double *p_mid);
function gsw_geo_strf_dyn_height_pc(sa, ct, delta_p, n_levels, geo_strf_dyn_height_pc, p_mid)::Cdouble
  return ccall(("gsw_geo_strf_dyn_height_pc",libgswteos),Cdouble,(Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Int,Ptr{Cdouble},Ptr{Cdouble}),sa, ct, delta_p, n_levels, geo_strf_dyn_height_pc, p_mid)
end

# The TEOS-10 Gibbs function of ice and its derivatives
#extern double gsw_gibbs_ice (int nt, int np, double t, double p);
function gsw_gibbs_ice(nt, np, t, p)::Cdouble
  return ccall(("gsw_gibbs_ice",libgswteos),Cdouble,(Int,Int,Cdouble,Cdouble),nt,np,t,p)
end

# Part of gibbs_ice(1,0,t,p)
#extern double gsw_gibbs_ice_part_t(double t, double p);
function gsw_gibbs_ice_part_t(t, p)::Cdouble
  return ccall(("gsw_gibbs_ice_part_t",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Part of gibbs_ice(1,0,pt0,0)
#extern double gsw_gibbs_ice_pt0(double pt0);
function gsw_gibbs_ice_pt0(pt0)::Cdouble
  return ccall(("gsw_gibbs_ice_pt0",libgswteos),Cdouble,(Cdouble,),pt0)
end

#
#extern double gsw_gibbs_ice_pt0_pt0(double pt0);
function gsw_gibbs_ice_pt0_pt0(pt0)::Cdouble
#  return ccall(("gsw_gibbs_ice_pt0_pt0",libgswteos),Cdouble,(Cdouble),pt0)
end

# The TEOS-10 Gibbs function of seawater and its derivatives
#extern double gsw_gibbs(int ns, int nt, int np, double sa, double t, double p);
function gsw_gibbs(ns, nt, np, sa, t, p)::Cdouble
  return ccall(("gsw_gibbs",libgswteos),Cdouble,(Int,Int,Int,Cdouble,Cdouble,Cdouble),ns,nt,np,sa,t,p)
end

# Gibbs(0,2,0,SA,t,0)
#extern double gsw_gibbs_pt0_pt0(double sa, double pt0);
function gsw_gibbs_pt0_pt0(sa, pt0)::Cdouble
  return ccall(("gsw_gibbs_pt0_pt0",libgswteos),Cdouble,(Cdouble,Cdouble),sa,pt0)
end

# Gravitational acceleration
#extern double gsw_grav(double lat, double p);
function gsw_grav(lat, p)::Cdouble
  return ccall(("gsw_grav",libgswteos),Cdouble,(Cdouble,Cdouble),lat,p)
end

# Helmholtz energy of ice
#extern double gsw_helmholtz_energy_ice(double t, double p);
function gsw_helmholtz_energy_ice(t, p)::Cdouble
  return ccall(("gsw_helmholtz_energy_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Hill ratio at a Practical Salinity of 2
#extern double gsw_hill_ratio_at_sp2(double t);
function gsw_hill_ratio_at_sp2(t)::Cdouble
  return ccall(("gsw_hill_ratio_at_sp2",libgswteos),Cdouble,(Cdouble,),t)
end

# Ice mass fraction to freeze seawater
#extern void   gsw_ice_fraction_to_freeze_seawater(double sa, double ct,double p, double t_ih, double *sa_freeze, double *ct_freeze,double *w_ih);
function gsw_ice_fraction_to_freeze_seawater(sa, ct, p, t_ih)
  ccall(("gsw_ice_fraction_to_freeze_seawater",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,t_ih,sa_freeze,ct_freeze,w_ih)
  return sa_freeze[],ct_freeze[],w_ih[]
end

# Internal energy
#extern double gsw_internal_energy(double sa, double ct, double p);
function gsw_internal_energy(sa, ct, p)::Cdouble
  return ccall(("gsw_internal_energy",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Internal energy of ice
#extern double gsw_internal_energy_ice(double t, double p);
function gsw_internal_energy_ice(t, p)::Cdouble
  return ccall(("gsw_internal_energy_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

#extern void   gsw_ipv_vs_fnsquared_ratio(double *sa, double *ct, double *p,double p_ref, int nz, double *ipv_vs_fnsquared_ratio,double *p_mid);
function gsw_ipv_vs_fnsquared_ratio(sa, ct, p, p_ref, nz)
  ccall(("gsw_ipv_vs_fnsquared_ratio",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Cdouble, Int, Ptr{Cdouble}, Ptr{Cdouble}),sa, ct, p, p_ref, nz, ipv_vs_fnsquared_ratio, p_mid)
  return ipv_vs_fnsquared_ratio[], p_mid[]
end

# Isothermal compressibility of ice
#extern double gsw_kappa_const_t_ice(double t, double p);
function gsw_kappa_const_t_ice(t, p)::Cdouble
  return ccall(("gsw_kappa_const_t_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Isentropic compressibility
#extern double gsw_kappa(double sa, double ct, double p);
function gsw_kappa(sa, ct, p)::Cdouble
  return ccall(("gsw_kappa",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Isentropic compressibility of ice
#extern double gsw_kappa_ice(double t, double p);
function gsw_kappa_ice(t, p)::Cdouble
  return ccall(("gsw_kappa_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Isentropic compressibility
#extern double gsw_kappa_t_exact(double sa, double t, double p);
function gsw_kappa_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_kappa_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Latent heat of evaporation of water from seawater (isobaric evaporation enthalpy) with CT as input temperature
#extern double gsw_latentheat_evap_ct(double sa, double ct);
function gsw_latentheat_evap_ct(sa, ct)::Cdouble
  return ccall(("gsw_latentheat_evap_ct",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# Latent heat of evaporation of water from seawater (isobaric evaporation enthalpy) with in-situ temperature, t, as input
#extern double gsw_latentheat_evap_t(double sa, double t);
function gsw_latentheat_evap_t(sa, t)::Cdouble
  return ccall(("gsw_latentheat_evap_t",libgswteos),Cdouble,(Cdouble,Cdouble),sa,t)
end

# Latent heat of melting of ice into seawater (isobaric melting enthalpy)
#extern double gsw_latentheat_melting(double sa, double p);
function gsw_latentheat_melting(sa, p)::Cdouble
  return ccall(("gsw_latentheat_melting",libgswteos),Cdouble,(Cdouble,Cdouble),sa,p)
end

#extern void   gsw_linear_interp_sa_ct(double *sa, double *ct, double *p, int np,double *p_i, int npi, double *sa_i, double *ct_i);
function gsw_linear_interp_sa_ct(sa, ct, p, np, p_i, npi)
  ccall(("gsw_linear_interp_sa_ct",libgswteos),Cvoid,(Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Int, Ptr{Cdouble}, Int, Ptr{Cdouble}, Ptr{Cdouble}),ssa, ct, p, np, p_i, npi, sa_i, ct_i)
  return sa_i[], ct[i]
end

# SA to CT ratio when ice melts into seawater near equilibrium
#extern double gsw_melting_ice_equilibrium_sa_ct_ratio(double sa, double p);
function gsw_melting_ice_equilibrium_sa_ct_ratio(sa, p)::Cdouble
  return ccall(("gsw_melting_ice_equilibrium_sa_ct_ratio",libgswteos),Cdouble,(Cdouble,Cdouble),sa,p)
end

# SA to CT ratio when ice melts into seawater near equilibrium (polynomial)
#extern double gsw_melting_ice_equilibrium_sa_ct_ratio_poly(double sa, double p);
function gsw_melting_ice_equilibrium_sa_ct_ratio_poly(sa, p)::Cdouble
  return ccall(("gsw_melting_ice_equilibrium_sa_ct_ratio_poly",libgswteos),Cdouble,(Cdouble,Cdouble),sa,p)
end

# SA and CT when ice melts in seawater
#extern void   gsw_melting_ice_into_seawater(double sa, double ct, double p,double w_ih, double t_ih, double *sa_final, double *ct_final,double *w_ih_final);
function gsw_melting_ice_into_seawater(sa, ct, p, w_ih, t_ih)
  ccall(("gsw_melting_ice_into_seawater",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,w_ih,t_ih,sa_final,ct_final,w_ih_final)
  return sa_final[],ct_final[],w_ih_final[]
end

# SA to CT ratio when ice melts in seawater
#extern double gsw_melting_ice_sa_ct_ratio(double sa, double ct, double p,double t_ih);
function gsw_melting_ice_sa_ct_ratio(sa, ct, p, t_ih)::Cdouble
  return ccall(("gsw_melting_ice_sa_ct_ratio",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sa,ct,p,t_ih)
end

# SA to CT ratio when ice melts in seawater (polynomial)
#extern double gsw_melting_ice_sa_ct_ratio_poly(double sa, double ct, double p,double t_ih);
function gsw_melting_ice_sa_ct_ratio_poly(sa, ct, p, t_ih)::Cdouble
  return ccall(("gsw_melting_ice_sa_ct_ratio_poly",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sa,ct,p,t_ih)
end

# SA to CT ratio when sea ice melts into seawater equilibrium
#extern double gsw_melting_seaice_equilibrium_sa_ct_ratio(double sa, double p);
function gsw_melting_seaice_equilibrium_sa_ct_ratio(sa, p)::Cdouble
  return ccall(("gsw_melting_seaice_equilibrium_sa_ct_ratio",libgswteos),Cdouble,(Cdouble,Cdouble),sa,p)
end

#	SA to CT ratio when sea ice melts into seawater equilibrium (polynomial)
#extern double gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(double sa,double p);
function gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(sa, p)::Cdouble
  return ccall(("gsw_melting_seaice_equilibrium_sa_ct_ratio_poly",libgswteos),Cdouble,(Cdouble,Cdouble),sa,p)
end

# SA and CT when sea ice melts in seawater
#extern void   gsw_melting_seaice_into_seawater(double sa, double ct, double p,double w_seaice, double sa_seaice, double t_seaice,double *sa_final, double *ct_final);
function gsw_melting_seaice_into_seawater(sa, ct, p, w_seaice, sa_seaice, t_seaice)
  ccall(("gsw_melting_seaice_into_seawater",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,w_seaice,sa_seaice,t_seaice,sa_final,ct_final)
  return sa_final[],ct_final[]
end

# SA to CT ratio when sea ice melts in seawater
#extern double gsw_melting_seaice_sa_ct_ratio(double sa, double ct, double p,double sa_seaice, double t_seaice);
function gsw_melting_seaice_sa_ct_ratio(sa, ct, p, sa_seaice, t_seaice)::Cdouble
  return ccall(("gsw_melting_seaice_sa_ct_ratio",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble,Cdouble),sa,ct,p,sa_seaice,t_seaice)
end

# SA to CT ratio when sea ice melts in seawater (polynomial)
#extern double gsw_melting_seaice_sa_ct_ratio_poly(double sa, double ct,double p, double sa_seaice, double t_seaice);
function gsw_melting_seaice_sa_ct_ratio_poly(sa, ct, p, sa_seaice, t_seaice)::Cdouble
  return ccall(("gsw_melting_seaice_sa_ct_ratio_poly",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble,Cdouble),sa,ct,p,sa_seaice,t_seaice)
end

#extern void   gsw_nsquared(double *sa, double *ct, double *p, double *lat,int nz, double *n2, double *p_mid);
function gsw_nsquared(sa,ct,p,lat,nz)
  ccall(("gsw_nsquared",libgswteos),Cvoid,(Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Int, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,lat,nz,n2,p_mld)
  return n2[],p_mld[]
end

#extern double gsw_o2sol(double sa, double ct, double p, double lon, double lat);
function gsw_o2sol(sa,ct,p,lon,lat)::Cdouble
  return ccall(("gsw_o2sol",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble,Cdouble),sa,ct,p,lon,lat)
end

#extern double gsw_o2sol_sp_pt(double sp, double pt);
function gsw_o2sol_sp_pt(sp,pt)::Cdouble
  return ccall(("gsw_o2sol_sp_pt",libgswteos),Cdouble,(Cdouble,Cdouble),sp,pt)
end

# Potential enthalpy from potential temperature of ice
#extern double gsw_pot_enthalpy_from_pt_ice(double pt0_ice);
function gsw_pot_enthalpy_from_pt_ice(pt0_ice)::Cdouble
  return ccall(("gsw_pot_enthalpy_from_pt_ice",libgswteos),Cdouble,(Cdouble,),pt0_ice)
end

# Potential enthalpy from potential temperature of ice (polynomial)
#extern double gsw_pot_enthalpy_from_pt_ice_poly(double pt0_ice);
function gsw_pot_enthalpy_from_pt_ice_poly(pt0_ice)::Cdouble
  return ccall(("gsw_pot_enthalpy_from_pt_ice_poly",libgswteos),Cdouble,(Cdouble,),pt0_ice)
end

# Potential enthalpy of ice at which seawater freezes
#extern double gsw_pot_enthalpy_ice_freezing(double sa, double p);
function gsw_pot_enthalpy_ice_freezing(sa, p)::Cdouble
  return ccall(("gsw_pot_enthalpy_ice_freezing",libgswteos),Cdouble,(Cdouble,Cdouble),sa,p)
end

# First derivatives of potential enthalpy of ice at which seawater freezes
#extern void   gsw_pot_enthalpy_ice_freezing_first_derivatives(double sa,double p, double *pot_enthalpy_ice_freezing_sa,double *pot_enthalpy_ice_freezing_p);
function gsw_pot_enthalpy_ice_freezing_first_derivatives(sa, p)
  ccall(("gsw_pot_enthalpy_ice_freezing_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,p,pot_enthalpy_ice_freezing_sa,pot_enthalpy_ice_freezing_p)
  return pot_enthalpy_ice_freezing_sa[],pot_enthalpy_ice_freezing_p[]
end

# First derivatives of potential enthalpy of ice at which seawater freezes (polynomial)
#extern void   gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(double sa,double p, double *pot_enthalpy_ice_freezing_sa,double *pot_enthalpy_ice_freezing_p);
function gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(sa, p)
  ccall(("gsw_pot_enthalpy_ice_freezing_first_derivatives_poly",libgswteos),Cvoid,(Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,p,pot_enthalpy_ice_freezing_sa,pot_enthalpy_ice_freezing_p)
  return pot_enthalpy_ice_freezing_sa[],pot_enthalpy_ice_freezing_p[]
end

# Potential enthalpy of ice at which seawater freezes (polynomial)
#extern double gsw_pot_enthalpy_ice_freezing_poly(double sa, double p);
function gsw_pot_enthalpy_ice_freezing_poly(sa, p)::Cdouble
  return ccall(("gsw_pot_enthalpy_ice_freezing_poly",libgswteos),Cdouble,(Cdouble,Cdouble),sa,p)
end

# Potential density
#extern double gsw_pot_rho_t_exact(double sa, double t, double p, double p_ref);
function gsw_pot_rho_t_exact(sa, t, p, p_ref)::Cdouble
  return ccall(("gsw_pot_rho_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sa,t,p,p_ref)
end

# Pressure coefficient of ice
#extern double gsw_pressure_coefficient_ice(double t, double p);
function gsw_pressure_coefficient_ice(t, p)::Cdouble
  return ccall(("gsw_pressure_coefficient_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Pressure of seawater at the freezing temperature (for given CT)
#extern double gsw_pressure_freezing_ct(double sa, double ct,double saturation_fraction);
function gsw_pressure_freezing_ct(sa, ct, saturation_fraction)::Cdouble
  return ccall(("gsw_pressure_freezing_ct",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,saturation_fraction)
end

# Calculates an initial estimate of pt0_ice when it is less than about -100 deg C.
#extern double gsw_pt0_cold_ice_poly(double pot_enthalpy_ice);
function gsw_pt0_cold_ice_poly(pot_enthalpy_ice)::Cdouble
  return ccall(("gsw_pt0_cold_ice_poly",libgswteos),Cdouble,(Cdouble,),pot_enthalpy_ice)
end

# Potential temperature with a reference pressure of zero dbar
#extern double gsw_pt0_from_t(double sa, double t, double p);
function gsw_pt0_from_t(sa, t, p)::Cdouble
  return ccall(("gsw_pt0_from_t",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Potential temperature of ice with reference pressure of 0 dbar
#extern double gsw_pt0_from_t_ice(double t, double p);
function gsw_pt0_from_t_ice(t, p)::Cdouble
  return ccall(("gsw_pt0_from_t_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# First derivatives of potential temperature
#extern void   gsw_pt_first_derivatives (double sa, double ct, double *pt_sa,double *pt_ct);
function gsw_pt_first_derivatives(sa, ct)
  ccall(("gsw_pt_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,pt_sa,pt_ct)
  return pt_sa[],pt_ct[]
end

# Potential temperature from Conservative Temperature
#extern double gsw_pt_from_ct(double sa, double ct);
function gsw_pt_from_ct(sa, ct)::Cdouble
  return ccall(("gsw_pt_from_ct",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# Potential temperature from entropy
#extern double gsw_pt_from_entropy(double sa, double entropy);
function gsw_pt_from_entropy(sa, entropy)::Cdouble
  return ccall(("gsw_pt_from_entropy",libgswteos),Cdouble,(Cdouble,Cdouble),sa,entropy)
end

# Potential temperature from potential enthalpy of ice
#extern double gsw_pt_from_pot_enthalpy_ice(double pot_enthalpy_ice);
function gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_ice)::Cdouble
  return ccall(("gsw_pt_from_pot_enthalpy_ice",libgswteos),Cdouble,(Cdouble,),pot_enthalpy_ice)
end

# Derivative of potential temperature of ice with respect to potential enthalpy
#extern double gsw_pt_from_pot_enthalpy_ice_poly_dh(double pot_enthalpy_ice);
function gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice)::Cdouble
  return ccall(("gsw_pt_from_pot_enthalpy_ice_poly_dh",libgswteos),Cdouble,(Cdouble,),pot_enthalpy_ice)
end

# Potential temperature from potential enthalpy of ice (polynomial)
#extern double gsw_pt_from_pot_enthalpy_ice_poly(double pot_enthalpy_ice);
function gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice)::Cdouble
  return ccall(("gsw_pt_from_pot_enthalpy_ice_poly",libgswteos),Cdouble,(Cdouble,),pot_enthalpy_ice)
end

# Potential temperature
#extern double gsw_pt_from_t(double sa, double t, double p, double p_ref);
function gsw_pt_from_t(sa, t, p, p_ref)::Cdouble
  return ccall(("gsw_pt_from_t",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sa,t,p,p_ref)
end

# Potential temperature of ice
#extern double gsw_pt_from_t_ice(double t, double p, double p_ref);
function gsw_pt_from_t_ice(t, p, p_ref)::Cdouble
  return ccall(("gsw_pt_from_t_ice",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),t,p,p_ref)
end

# Second derivatives of potential temperature
#extern void   gsw_pt_second_derivatives (double sa, double ct, double *pt_sa_sa,double *pt_sa_ct, double *pt_ct_ct);
function gsw_pt_second_derivatives(sa, ct)
  ccall(("gsw_pt_second_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,pt_sa_sa,pt_sa_ct,pt_ct_ct)
  return pt_sa_sa[],pt_sa_ct[],pt_ct_ct[]
end

# in-situ density, thermal expansion & saline contraction coefficients
#extern void   gsw_rho_alpha_beta (double sa, double ct, double p, double *rho,double *alpha, double *beta);
function gsw_rho_alpha_beta(sa, ct, p)
  ccall(("gsw_rho_alpha_beta",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa,ct,p,rho,alpha,beta)
  return rho[],alpha[],beta[]
end

# in-situ density and potential density
#extern double gsw_rho(double sa, double ct, double p);
function gsw_rho(sa, ct, p)::Cdouble
  return ccall(("gsw_rho",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

#extern void   gsw_rho_first_derivatives(double sa, double ct, double p,double *drho_dsa, double *drho_dct, double *drho_dp);
function gsw_rho_first_derivatives(sa, ct, p)
  ccall(("gsw_rho_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa, ct, p, drho_dsa, drho_dct, drho_dp)
  return drho_dsa[], drho_dct[], drho_dp[]
end

#extern void   gsw_rho_first_derivatives_wrt_enthalpy (double sa, double ct,double p, double *rho_sa, double *rho_h);
function gsw_rho_first_derivatives_wrt_enthalpy(sa, ct, p)
  ccall(("gsw_rho_first_derivatives_wrt_enthalpy",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa, ct, p, rho_sa, rho_h)
  return rho_sa[], rho_h[]
end

# in-situ density of ice
#extern double gsw_rho_ice(double t, double p);
function gsw_rho_ice(t, p)::Cdouble
  return ccall(("gsw_rho_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

#extern void   gsw_rho_second_derivatives(double sa, double ct, double p,double *rho_sa_sa, double *rho_sa_ct, double *rho_ct_ct,double *rho_sa_p, double *rho_ct_p);
function gsw_rho_second_derivatives(sa, ct, p)
  ccall(("gsw_rho_second_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa, ct, p, rho_sa_sa, rho_sa_ct, rho_ct_ct, rho_sa_p, rho_ct_p)
  return rho_sa_sa[], rho_sa_ct[], rho_ct_ct[], rho_sa_p[], rho_ct_p[]
end

#extern void   gsw_rho_second_derivatives_wrt_enthalpy(double sa, double ct,double p, double *rho_sa_sa, double *rho_sa_h, double *rho_h_h);
function gsw_rho_second_derivatives_wrt_enthalpy(sa, ct, p)
  ccall(("gsw_rho_second_derivatives_wrt_enthalpy",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa, ct, p, rho_sa_sa, rho_sa_h, rho_h_h)
  return rho_sa_sa[], rho_sa_h[], rho_h_h[]
end

# in-situ density (Laboratory function, for use with densimeter measurements)
#extern double gsw_rho_t_exact(double sa, double t, double p);
function gsw_rho_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_rho_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

#extern void   gsw_rr68_interp_sa_ct(double *sa, double *ct, double *p, int mp,double *p_i, int mp_i, double *sa_i, double *ct_i);
function gsw_rr68_interp_sa_ct(sa, ct, p, mp, p_i, mp_i)
  ccall(("gsw_rr68_interp_sa_ct",libgswteos),Cvoid,(Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Int, Ptr{Cdouble}, Int, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, mp, p_i, mp_i, sa_i, ct_i)
  return sa_i[], ct_i[]
end

# Absolute Salinity Anomaly Ratio (excluding the Baltic Sea)
#extern double gsw_saar(double p, double lon, double lat);
function gsw_saar(p, lon, lat)::Cdouble
  return ccall(("gsw_saar",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),p,lon,lat)
end

#extern double gsw_sa_freezing_estimate(double p, double saturation_fraction,double *ct, double *t)
function gsw_sa_freezing_estimate(p, saturation_fraction, ct, t)::Cdouble
  return ccall(("gsw_sa_freezing_estimate",libgswteos),Cdouble,(Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),p,saturation_fraction,ct,t)
end

# Absolute Salinity of seawater at the freezing point (for given CT)
#extern double gsw_sa_freezing_from_ct(double ct, double p,double saturation_fraction);
function gsw_sa_freezing_from_ct(ct, p, saturation_fraction)::Cdouble
  return ccall(("gsw_sa_freezing_from_ct",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),ct,p,saturation_fraction)
end

# Absolute Salinity of seawater at the freezing point (for given CT) (polynomial)
#extern double gsw_sa_freezing_from_ct_poly(double ct, double p,double saturation_fraction);
function gsw_sa_freezing_from_ct_poly(ct, p, saturation_fraction)::Cdouble
  return ccall(("gsw_sa_freezing_from_ct_poly",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),ct,p,saturation_fraction)
end

# Absolute Salinity of seawater at the freezing point (for given t)
#extern double gsw_sa_freezing_from_t(double t, double p,double saturation_fraction);
function gsw_sa_freezing_from_t(t, p, saturation_fraction)::Cdouble
  return ccall(("gsw_sa_freezing_from_t",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),t,p,saturation_fraction)
end

# Absolute Salinity of seawater at the freezing point (for given t) (polynomial)
#extern double gsw_sa_freezing_from_t_poly(double t, double p,double saturation_fraction);
function gsw_sa_freezing_from_t_poly(t, p, saturation_fraction)::Cdouble
  return ccall(("gsw_sa_freezing_from_t_poly",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),t,p,saturation_fraction)
end

# Absolute Salinity from density
#extern double gsw_sa_from_rho(double rho, double ct, double p);
function gsw_sa_from_rho(rho, ct, p)::Cdouble
  return ccall(("gsw_sa_from_rho",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),rho,ct,p)
end

# Calculates Absolute Salinity in the Baltic Sea
#extern double gsw_sa_from_sp_baltic(double sp, double lon, double lat);
function gsw_sa_from_sp_baltic(sp, lon, lat)::Cdouble
  return ccall(("gsw_sa_from_sp_baltic",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sp,lon,lat)
end

# Absolute Salinity from Practical Salinity
#extern double gsw_sa_from_sp(double sp, double p, double lon, double lat);
function gsw_sa_from_sp(sp, p, lon, lat)::Cdouble
  return ccall(("gsw_sa_from_sp",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sp,p,lon,lat)
end

# Absolute Salinity from Preformed Salinity
#extern double gsw_sa_from_sstar(double sstar, double p,double lon,double lat);
function gsw_sa_from_sstar(sstar, p, lon, lat)::Cdouble
  return ccall(("gsw_sa_from_sstar",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sstar,p,lon,lat)
end

# Check for any values that are out of the TEOS-10 range ... [0..1]
#extern int    gsw_sa_p_inrange(double sa, double p);
function gsw_sa_p_inrange(sa, p)::Int
  return ccall(("gsw_sa_p_inrange",libgswteos),Int,(Cdouble,Cdouble),sa,p)
end

#extern void   gsw_seaice_fraction_to_freeze_seawater(double sa, double ct,double p, double sa_seaice, double t_seaice, double *sa_freeze,double *ct_freeze, double *w_seaice);
function gsw_seaice_fraction_to_freeze_seawater(sa, ct, p, sa_seaice, t_seaice)
  ccall(("gsw_seaice_fraction_to_freeze_seawater",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa, ct, p, sa_seaice, t_seaice, sa_freeze, ct_freeze, w_seaice)
  return sa_freeze[], ct_freeze[], w_seaice[]
end

# sigma0 with reference pressure of 0 dbar
#extern double gsw_sigma0(double sa, double ct);
function gsw_sigma0(sa, ct)::Cdouble
  return ccall(("gsw_sigma0",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# sigma1 with reference pressure of 1000 dbar
#extern double gsw_sigma1(double sa, double ct);
function gsw_sigma1(sa, ct)::Cdouble
  return ccall(("gsw_sigma1",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# sigma2 with reference pressure of 2000 dbar
#extern double gsw_sigma2(double sa, double ct);
function gsw_sigma2(sa, ct)::Cdouble
  return ccall(("gsw_sigma2",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# sigma3 with reference pressure of 3000 dbar
#extern double gsw_sigma3(double sa, double ct);
function gsw_sigma3(sa, ct)::Cdouble
  return ccall(("gsw_sigma3",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# sigma4 with reference pressure of 4000 dbar
#extern double gsw_sigma4(double sa, double ct);
function gsw_sigma4(sa, ct)::Cdouble
  return ccall(("gsw_sigma4",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# Sound speed
#extern double gsw_sound_speed(double sa, double ct, double p);
function gsw_sa_from_rho(sa, ct, p)::Cdouble
  return ccall(("gsw_sa_from_rho",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

#extern double gsw_sound_speed(double sa, double ct, double p);
function gsw_sound_speed(sa, ct, p)::Cdouble
  return ccall(("gsw_sound_speed",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Sound speed of ice (compression waves)
#extern double gsw_sound_speed_ice(double t, double p);
function gsw_sound_speed_ice(t, p)::Cdouble
  return ccall(("gsw_sound_speed_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

# Sound speed
#extern double gsw_sound_speed_t_exact(double sa, double t, double p);
function gsw_sound_speed_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_sound_speed_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

#extern void   gsw_specvol_alpha_beta(double sa, double ct, double p,double *specvol, double *alpha, double *beta);
function gsw_specvol_alpha_beta(sa, ct, p)
  ccall(("gsw_specvol_alpha_beta",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, specvol, alpha, beta)
  return specvol[], alpha[], beta[]
end

# Specific volume anomaly relative to SSO & 0°C
#extern double gsw_specvol_anom_standard(double sa, double ct, double p);
function gsw_specvol_anom_standard(sa, ct, p)::Cdouble
  return ccall(("gsw_specvol_anom_standard",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# Specific volume
#extern double gsw_specvol(double sa, double ct, double p);
function gsw_specvol(sa, ct, p)::Cdouble
  return ccall(("gsw_specvol",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

#extern void   gsw_specvol_first_derivatives(double sa, double ct, double p,double *v_sa, double *v_ct, double *v_p);
function gsw_specvol_first_derivatives(sa, ct, p)
  ccall(("gsw_specvol_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa, v_ct, v_p)
  return v_sa[], v_ct[], v_p[]
end

#extern void   gsw_specvol_first_derivatives_wrt_enthalpy(double sa, double ct,double p, double *v_sa, double *v_h);
function gsw_specvol_first_derivatives_wrt_enthalpy(sa, ct, p)
  ccall(("gsw_specvol_first_derivatives_wrt_enthalpy",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), sa, ct, p, v_sa, v_h)
  return v_sa[], v_h[]
end

# Specific volume of ice
#extern double gsw_specvol_ice(double t, double p);
function gsw_specvol_ice(t, p)::Cdouble
  return ccall(("gsw_specvol_ice",libgswteos),Cdouble,(Cdouble,Cdouble),t,p)
end

#extern void   gsw_specvol_second_derivatives (double sa, double ct, double p,double *v_sa_sa, double *v_sa_ct, double *v_ct_ct,double *v_sa_p, double *v_ct_p);
function gsw_specvol_second_derivatives(sa, ct, p)
  ccall(("gsw_specvol_second_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa, ct, p, v_sa_sa, v_sa_ct, v_ct_ct, v_sa_p, v_ct_p)
  return v_sa_sa[], v_sa_ct[], v_ct_ct[], v_sa_p[], v_ct_p[]
end

#extern void   gsw_specvol_second_derivatives_wrt_enthalpy(double sa, double ct,double p, double *v_sa_sa, double *v_sa_h, double *v_h_h);
function gsw_specvol_second_derivatives_wrt_enthalpy(sa, ct, p)
  ccall(("gsw_specvol_second_derivatives_wrt_enthalpy",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),sa, ct, p, v_sa_sa, v_sa_h, v_h_h)
  return v_sa_sa[], v_sa_h[], v_h_h[]
end

# Specific volume specvol(35.16504,0,p)
#extern double gsw_specvol_sso_0(double p);
function gsw_specvol_sso_0(p)::Cdouble
  return ccall(("gsw_specvol_sso_0",libgswteos),Cdouble,(Cdouble,),p)
end

# Specific volume
#extern double gsw_specvol_t_exact(double sa, double t, double p);
function gsw_specvol_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_specvol_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# Practical Salinity from conductivity, C (inc. for SP < 2)
#extern double gsw_sp_from_c(double c, double t, double p);
function gsw_sp_from_c(c, t, p)::Cdouble
  return ccall(("gsw_sp_from_c",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),c,t,p)
end

# Calculates Absolute Salinity in the Baltic Sea
#extern double gsw_sp_from_sa_baltic(double sa, double lon, double lat);
function gsw_sp_from_sa_baltic(sa, lon, lat)::Cdouble
  return ccall(("gsw_sp_from_sa_baltic",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,lon,lat)
end

# Practical Salinity from Absolute Salinity
#extern double gsw_sp_from_sa(double sa, double p, double lon, double lat);
function gsw_sp_from_sa(sa, p, lon, lat)::Cdouble
  return ccall(("gsw_sp_from_sa",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sa,p,lon,lat)
end

# Practical Salinity from Knudsen Salinity
#extern double gsw_sp_from_sk(double sk);
function gsw_sp_from_sk(sk)::Cdouble
  return ccall(("gsw_sp_from_sk",libgswteos),Cdouble,(Cdouble,),sk)
end

# Practical Salinity from Reference Salinity
#extern double gsw_sp_from_sr(double sr);
function gsw_sp_from_sr(sr)::Cdouble
  return ccall(("gsw_sp_from_sr",libgswteos),Cdouble,(Cdouble,),sr)
end

# Practical Salinity from Preformed Salinity
#extern double gsw_sp_from_sstar(double sstar, double p,double lon,double lat);
function gsw_sp_from_sstar(sstar, p, lon, lat)::Cdouble
  return ccall(("gsw_sp_from_sstar",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sstar,p,lon,lat)
end

#extern double gsw_sp_salinometer(double rt, double t);
function gsw_sp_salinometer(rt, t)::Cdouble
  return ccall(("gsw_sp_salinometer",libgswteos),Cdouble,(Cdouble,Cdouble),rt,t)
end

# Spiciness with p_ref of 0 dbar
#extern double gsw_spiciness0(double sa, double ct);
function gsw_spiciness0(sa, ct)::Cdouble
  return ccall(("gsw_spiciness0",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# Spiciness with p_ref of 1000 dbar
#extern double gsw_spiciness1(double sa, double ct);
function gsw_spiciness1(sa, ct)::Cdouble
  return ccall(("gsw_spiciness1",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# Spiciness with p_ref of 2000 dbar
#extern double gsw_spiciness2(double sa, double ct);
function gsw_spiciness2(sa, ct)::Cdouble
  return ccall(("gsw_spiciness2",libgswteos),Cdouble,(Cdouble,Cdouble),sa,ct)
end

# Reference Salinity from Practical Salinity
#extern double gsw_sr_from_sp(double sp);
function gsw_sr_from_sp(sp)::Cdouble
  return ccall(("gsw_sr_from_sp",libgswteos),Cdouble,(Cdouble,),sp)
end

# Preformed Salinity from Absolute Salinity
#extern double gsw_sstar_from_sa(double sa, double p, double lon, double lat);
function gsw_sstar_from_sa(sa, p, lon, lat)::Cdouble
  return ccall(("gsw_sstar_from_sa",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sa,p,lon,lat)
end

# Preformed Salinity from Practical Salinity
#extern double gsw_sstar_from_sp(double sp, double p, double lon, double lat);
function gsw_sstar_from_sp(sp, p, lon, lat)::Cdouble
  return ccall(("gsw_sstar_from_sp",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),sp,p,lon,lat)
end

# Temperature derivative of chemical potential of water
#extern double gsw_t_deriv_chem_potential_water_t_exact(double sa, double t,double p);
function gsw_t_deriv_chem_potential_water_t_exact(sa, t, p)::Cdouble
  return ccall(("gsw_t_deriv_chem_potential_water_t_exact",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,t,p)
end

# in-situ freezing temperature of seawater
#extern double gsw_t_freezing(double sa, double p, double saturation_fraction);
function gsw_t_freezing(sa, p, saturation_fraction)::Cdouble
  return ccall(("gsw_t_freezing",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,p,saturation_fraction)
end

#extern void   gsw_t_freezing_first_derivatives_poly(double sa, double p,double saturation_fraction, double *tfreezing_sa,double *tfreezing_p);
function gsw_t_freezing_first_derivatives_poly(sa, p, saturation_fraction)
  ccall(("gsw_t_freezing_first_derivatives_poly",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,p,saturation_fraction,tfreezing_sa,freezing_p)
  return tfreezing_sa[],freezing_p[]
end

#extern void   gsw_t_freezing_first_derivatives(double sa, double p,double saturation_fraction, double *tfreezing_sa,double *tfreezing_p);
function gsw_t_freezing_first_derivatives(sa, p, saturation_fraction)
  ccall(("gsw_t_freezing_first_derivatives",libgswteos),Cvoid,(Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),sa,p,saturation_fraction,tfreezing_sa,freezing_p)
  return tfreezing_sa[],freezing_p[]
end

# in-situ freezing temperature of seawater (polynomial)
#extern double gsw_t_freezing_poly(double sa, double p,double saturation_fraction);
function gsw_t_freezing_poly(sa, p, saturation_fraction)::Cdouble
  return ccall(("gsw_t_freezing_poly",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,p,saturation_fraction)
end

# in-situ temperature from Conservative Temperature
#extern double gsw_t_from_ct(double sa, double ct, double p);
function gsw_t_from_ct(sa, ct, p)::Cdouble
  return ccall(("gsw_t_from_ct",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

# in-situ temperature from potential temperature of ice with p_ref of 0 dbar
#extern double gsw_t_from_pt0_ice(double pt0_ice, double p);
function gsw_t_from_pt0_ice(pt0_ice, p)::Cdouble
  return ccall(("gsw_t_from_pt0_ice",libgswteos),Cdouble,(Cdouble,Cdouble),pt0_ice,p)
end

# Thermobaric coefficient
#extern double gsw_thermobaric(double sa, double ct, double p);
function gsw_thermobaric(sa, ct, p)::Cdouble
  return ccall(("gsw_thermobaric",libgswteos),Cdouble,(Cdouble,Cdouble,Cdouble),sa,ct,p)
end

#extern void   gsw_turner_rsubrho(double *sa, double *ct, double *p, int nz,double *tu, double *rsubrho, double *p_mid);
function gsw_turner_rsubrho(rarray, nx, iarray)
  ccall(("gsw_turner_rsubrho",libgswteos),Cvoid,(Ptr{Cdouble}, Cdouble, Ptr{Cdouble}), rarray, nx, iarray)
  return iarray[]
end

#
#extern int    gsw_util_indx(double *x, int n, double z);
function gsw_util_indx(x, n, z)::Int
  return ccall(("gsw_util_indx",libgswteos),Int,(Ptr{Cdouble},Int,Cdouble),x, n, z)
end

#extern double *gsw_util_interp1q_int(int nx, double *x, int *iy, int nxi,double *x_i, double *y_i);
function gsw_util_interp1q_int(nx, x, iy, nxi)
  return ccall(("gsw_util_interp1q_int",libgswteos),Cdouble,(Int, Ptr{Cdouble}, Ptr{Int}, Int, Ptr{Cdouble}, Ptr{Cdouble}), nx, x, iy, nxi, x_i, y_i)
end

#extern double *gsw_util_linear_interp(int nx, double *x, int ny, double *y,int nxi, double *x_i, double *y_i);
function gsw_util_linear_interp(nx, x, ny, y, nxi)
  return ccall(("gsw_util_linear_interp",libgswteos),Cdouble,(Int, Ptr{Cdouble}, Int, Ptr{Cdouble}, Int, Ptr{Cdouble}, Ptr{Cdouble}), nx, x, ny, y, nxi, x_i, y_i)
end

#extern void   gsw_util_sort_real(double *rarray, int nx, int *iarray);
function gsw_util_sort_real(rarray, nx, iarray)
  ccall(("gsw_util_sort_real",libgswteos),Cvoid,(Ptr{Cdouble}, Cdouble, Ptr{Cdouble}), rarray, nx, iarray)
  return iarray[]
end

#extern double gsw_util_xinterp1(double *x, double *y, int n, double x0);
function gsw_util_xinterp1(x, y, n, x0)::Cdouble
  return ccall(("gsw_util_xinterp1",libgswteos),Cdouble,(Ptr{Cdouble},Ptr{Cdouble}, Int, Cdouble),x, y, n, x0)
end


#extern int gsw_util_pchip_interp(double *x, double *y, int n,double *xi, double *yi, int ni);
function gsw_util_pchip_interp(x, y, n, xi, yi, ni)::Cdouble
  return ccall(("gsw_util_pchip_interp",libgswteos),Int,(Ptr{Cdouble},Ptr{Cdouble}, Int, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble), x, y, n, xi, yi, ni)
end

# Height from pressure
#extern double gsw_z_from_p(double p, double lat);
function gsw_z_from_p(p, lat)::Cdouble
  return ccall(("gsw_z_from_p",libgswteos),Cdouble,(Cdouble,Cdouble),p,lat)
end

# Pressure from height
#extern double gsw_p_from_z(double z, double lat);
function gsw_p_from_z(z, lat)::Cdouble
  return ccall(("gsw_p_from_z",libgswteos),Cdouble,(Cdouble,Cdouble),z,lat)
end

"""
Copyright (C) 2018 Alexander Smirnov (axline@mail.ru)

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License 3 as published by the Free Software
Foundation. See the GNU General Public License for more details
(http://www.gnu.org/licenses/).

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.
"""
