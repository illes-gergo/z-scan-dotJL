import Base.+, Base.*, Base./

struct userinputs
  Nx::Int
  Nt::Int
  cry::Int
  sigma_t::Float64
  sigma_x::Float64
  lambda0::Float64
  I0::Float64
  STR::String
  gamma::Float64
  dz::Float64
  z_end::Float64
  x::Vector
  t::Vector
end

struct gaussVars
  E0::Float64
  t::Array{Float64,2}
  sigma_t::Float64
  x::Array{Float64,2}
  sigma_x::Float64
  omega0::Float64
  gamma::Float64
  lambda0::Float64
  cry::Int
end

mutable struct compositeInput
  Akxo::Array{ComplexF64,2}
  ATHz_kx_o::Array{ComplexF64,2}
  ASH::Array{ComplexF64,2}
end

function +(a::compositeInput, b::compositeInput)
  return compositeInput(a.Akxo .+ b.Akxo, a.ATHz_kx_o .+ b.ATHz_kx_o, a.ASH .+ b.ASH)
end

function *(a::compositeInput, b::compositeInput)
  return compositeInput(a.Akxo .* b.Akxo, a.ATHz_kx_o .* b.ATHz_kx_o, a.ASH .* b.ASH)
end
function *(a::Float64, b::compositeInput)
  return compositeInput(a .* b.Akxo, a .* b.ATHz_kx_o, a .* b.ASH)
end
function *(a::Int, b::compositeInput)
  return compositeInput(a .* b.Akxo, a .* b.ATHz_kx_o, a .* b.ASH)
end

function +(a::compositeInput, b::Float64)
  return compositeInput(a.Akxo .+ b, a.ATHz_kx_o .+ b, a.ASH .+ b)
end

function /(a::compositeInput, b::Int)
  return compositeInput(a.Akxo ./ b, a.ATHz_kx_o ./ b, a.ASH ./ b)
end

struct fourierOperations
  fft_t_o
  fft_x_kx
  ifft_o_t
  ifft_kx_x
  fast_conv_plan
  fast_conv_fft_plan
end

struct naturalConstants
  e0::Float64
  c0::Float64
  function naturalConstants()
    new(8.854187817e-12, 3e8)
  end
end

struct pumpFieldConstants
  kx_omega::Array{Float64,2}
  kz_omega::Array{Float64,2}
end

struct THzFieldConstants
  alpha::Array{Float64,2}
  kz_omegaTHz::Array{Float64,2}
end

struct SHFieldConstants
  kx_omegaSHG::Array{Float64,2}
  kz_omegaSHG::Array{Float64,2}
end

struct runTimeConstants
  kxMax::Float64
  cx::Array{Float64,2}
  d_eff::Float64
  khi_eff::Float64
  dOmega::Float64
  padding::Array{Float64,2}
  SHG_SHIFT::Int
  ckx::Array{Float64,2}
  comega::Array{Float64,2}
  comegaTHz::Array{Float64,2}
  comegaSHG::Array{Float64,2}
  omegaMax::Float64
  lambda0::Float64
  omega0::Float64
  cry::Int
end

struct miscInputs
  FOPS::fourierOperations
  NC::naturalConstants
  RTC::runTimeConstants
  TFC::THzFieldConstants
  PFC::pumpFieldConstants
  SFC::SHFieldConstants
end
