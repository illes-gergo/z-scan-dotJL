using LazyGrids, FFTW, FourierTools, Base.Threads, Dates, HDF5

# FFT -> /omegaMAX ; IFFT -> * omegaMAX

#default(levels=100)
#gr()
include("valtozok.jl")
include("gauss_impulzus.jl")
include("diffegy_megoldo.jl")
include("differencial_egyenletek.jl")
include("fuggvenyek.jl")
include("typedefs.jl")

function runcalc()
  inputs = userinputs()

  c0 = 3e8
  khi_eff = 2 * deffTHz(inputs.cry)
  d_eff = deff(inputs.cry)
  e0 = 8.854187817e-12

  SHG_SHIFT = floor(Int, inputs.Nt / 4)

  tMax = inputs.t[end] - inputs.t[1]
  dt = inputs.t[2] - inputs.t[1]

  xMax = inputs.x[end] - inputs.x[1]
  dx = inputs.x[2] - inputs.x[1]

  dOmega = 2 * pi / tMax
  omega_ = range(0, length=inputs.Nt, step=dOmega)

  omegaMax = omega_[end] - omega_[1]
  omega = omega_ .- omegaMax ./ 2

  dkx = 2 * pi / xMax
  kx = range(0, length=inputs.Nx, step=dkx) .- inputs.Nx / 2 * dkx

  kxMax = kx[end] - kx[1]

  omega_diff = round(Int, 2 * pi * c0 / inputs.lambda0 / dOmega)
  omega0 = omega_diff .* dOmega
  lambda0 = 2 * pi * c0 / omega0


  (ct, cx) = ndgrid(inputs.t, inputs.x)
  (comega_, ckx) = ndgrid(omega, kx)

  comega = comega_ .+ omega0

  comegaSHG = comega_ .+ 2 * omega0

  comegaTHz = comega_ .- omega[1]

  clambda = c0 ./ comega * 2 * pi

  cLambdaSHG = c0 ./ comegaSHG * 2 * pi

  n = neo(clambda, 300, inputs.cry)

  k_omega = n .* comega ./ c0

  k_omegaSHG = neo(cLambdaSHG, 300, inputs.cry) .* comegaSHG ./ c0


  E0 = sqrt(2 * inputs.I0 / neo(lambda0, 300, inputs.cry) / e0 / c0)
  gaussInput = gaussVars(E0, ct, inputs.sigma_t, cx, inputs.sigma_x, omega0, inputs.gamma, lambda0, inputs.cry)
  #  Axt = gauss_impulzus_omega0(E0, inputs.sigma_t, inputs.sigma_x, lambda0, inputs.gamma, ct, cx)
  Axt = gauss_impulzus_omega0(gaussInput)

  #=  fft_t_o = plan_fft(Axt, 1)
    fft_x_kx = plan_fft(Axt, 2)
    ifft_o_t = plan_ifft(Axt, 1)
    ifft_kx_x = plan_ifft(Axt, 2) =#
  alpha = aTHzo(comegaTHz, 300, inputs.cry)
  padding = zeros(inputs.Nt, inputs.Nx)
  #fast_conv_plan, fast_conv_fft_plan = plan_fast_conv(Axt, Axt, padding)
  RTC = runTimeConstants(kxMax, cx, d_eff, khi_eff, dOmega, padding, SHG_SHIFT, ckx, comega, comegaTHz, comegaSHG, omegaMax, lambda0, omega0, inputs.cry)
  PFC = pumpFieldConstants(kz_omega=k_omega)
  SFC = SHFieldConstants(kz_omegaSHG=k_omegaSHG)

  FOPS = fourierOperations(plan_fft(Axt, 1), plan_fft(Axt, 2), plan_ifft(Axt, 1), plan_ifft(Axt, 2), plan_fast_conv(Axt, Axt, RTC)...)


  misc = miscInputs(FOPS, naturalConstants(), RTC, PFC, SFC)

  Axo = fftshift(FOPS.fft_t_o * Axt, 1) ./ omegaMax 
  Akxo = fftshift(FOPS.fft_x_kx * Axo / kxMax, 2)
  z = Array{Float64}(undef, floor(Int, inputs.z_end / inputs.dz))
  ASH = copy(zeros(size(Akxo)))

  A_kompozit = compositeInput(Akxo, ASH)

  z[1] = 0

  for ii in 1:(length(z)-1)
    A_kompozit, z[ii+1] = RK4M(z_scan_SHG, z[ii], A_kompozit, inputs.dz, misc)
    display(ii)
  end
end
