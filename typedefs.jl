import Base.+, Base.*, Base./

include("valtozok.jl")

@enum CRYSTAL LN = 0 LT = 1 ZnTe = 2 GaAs = 4

@kwdef struct gaussVars
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

@kwdef struct compositeInput
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

@kwdef struct fourierOperations
  fft_t_o
  fft_x_kx
  ifft_o_t
  ifft_kx_x
  fast_conv_plan
  fast_conv_fft_plan
end

@kwdef struct naturalConstants
  e0::Float64 = 8.854187817e-12
  c0::Float64 = 3e8
end

@kwdef struct pumpFieldConstants
  kz_omega::Array{Float64,2}
end

@kwdef struct THzFieldConstants
  alpha::Array{Float64,2}
  kz_omegaTHz::Array{Float64,2}
end

@kwdef struct SHFieldConstants
  kz_omegaSHG::Array{Float64,2}
end

@kwdef struct runTimeConstants
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

@kwdef struct runTimeConstantsZSCAN
  kxMax::Float64
  cx::Array{Float64,2}
  d_eff::Float64
  dOmega::Float64
  padding::Array{Float64,2}
  SHG_SHIFT::Int
  ckx::Array{Float64,2}
  comega::Array{Float64,2}
  comegaSHG::Array{Float64,2}
  omegaMax::Float64
  lambda0::Float64
  omega0::Float64
  cry::Int
  m_eff::Float64
  Erinf::Float64
  tsc::Float64
  n2::Float64
end

@kwdef struct miscInputs
  FOPS::fourierOperations
  NC::naturalConstants
  RTC::runTimeConstantsZSCAN
  PFC::pumpFieldConstants
  SFC::SHFieldConstants
  UIN::userinputs
end

@kwdef struct zscanInput
  Akxo::Array{ComplexF64}
  ASH::Array{ComplexF64}
end

function *(a::Number, b::zscanInput)
  b.Akxo .*= a
  b.ASH .*= a
  return b
end

function +(a::zscanInput, b::zscanInput)
  a.Akxo .+= b.Akxo
  a.ASH .+= b.ASH
  return a
end

function /(b::zscanInput, a::Number)
  b.Akxo ./= a
  b.ASH ./= a
  return b
end
