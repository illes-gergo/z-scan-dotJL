@kwdef struct userinputs
  Nx::Int = 512
  Nt::Int = 512
  cry::Int = Int(LT)
  sigma_t::Float64 = 1e-12
  sigma_x::Float64 = 1e-3
  lambda0::Float64 = 1030e-9
  I0::Float64 = 50e13
  STR::String = "test_calculation"
  dz::Float64 = 1e-6
  z_end::Float64 = 1e-3
  mpa_order::Int = 4
  betaN::Float64 = 1e-30
  x::Vector = range(-sigma_x, sigma_x, Nx) * 8
  t::Vector = range(-sigma_t, sigma_t, Nt) * 8
end
