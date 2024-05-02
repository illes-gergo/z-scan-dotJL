@kwdef struct userinputs
  Nx::Int = 400
  Nt::Int = 400
  cry::Int = Int(LT)
  sigma_t::Float64 = 50e-15
  sigma_x::Float64 = 35e-6
  lambda0::Float64 = 800e-9
  I0::Float64 = 200e13
  STR::String = "test_calculation"
  dz::Float64 = 1e-6
  z_end::Float64 = 0.5e-3 + dz
  mpa_order::Int = 4
  betaN::Float64 = 1e-59
  x::Vector = range(-sigma_x, sigma_x, Nx) * 10
  t::Vector = range(-sigma_t, sigma_t, Nt) * 10
  zRayleigh::Float64 = -1
end
