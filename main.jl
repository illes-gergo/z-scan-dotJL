using LazyGrids, FFTW, FourierTools, Base.Threads, Dates, HDF5

# FFT -> /omegaMAX ; IFFT -> * omegaMAX

#default(levels=100)
#gr()
include("gauss_impulzus.jl")
include("diffegy_megoldo.jl")
include("differencial_egyenletek.jl")
include("fuggvenyek.jl")
include("typedefs.jl")
include("valtozok.jl")

function runcalc()
  inputs::userinputs = setInput()

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
  kx0 = sin(inputs.gamma) ./ inputs.lambda0 .* 2 .* pi * neo(inputs.lambda0, 300, inputs.cry)

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
  kx_omega = real.(k_omega .* sin(inputs.gamma))
  kz_omega = real.(k_omega .* cos(inputs.gamma))

  k_omegaSHG = neo(cLambdaSHG, 300, inputs.cry) .* comegaSHG ./ c0
  kx_omegaSHG = real.(k_omegaSHG .* sin(inputs.gamma))
  kz_omegaSHG = real.(k_omegaSHG .* cos(inputs.gamma))

  nTHz = nTHzo(comegaTHz, 300, inputs.cry)

  k_omegaTHz = nTHz .* comegaTHz ./ c0
  kz_omegaTHz = real.(sqrt.(Complex.(k_omegaTHz .^ 2 - ckx .^ 2)))

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
  TFC = THzFieldConstants(alpha, kz_omegaTHz)
  PFC = pumpFieldConstants(kx_omega, kz_omega)
  SFC = SHFieldConstants(kx_omegaSHG, kz_omegaSHG)

  FOPS = fourierOperations(plan_fft(Axt, 1), plan_fft(Axt, 2), plan_ifft(Axt, 1), plan_ifft(Axt, 2), plan_fast_conv(Axt, Axt, RTC)...)


  misc = miscInputs(FOPS, naturalConstants(), RTC, TFC, PFC, SFC)

  Axo = fftshift(FOPS.fft_t_o * Axt, 1) ./ omegaMax .* exp.(+1im .* kx_omega .* cx)
  Akxo = fftshift(FOPS.fft_x_kx * Axo / kxMax, 2)
  z = Array{Float64}(undef, floor(Int, inputs.z_end / inputs.dz))
  #error("stop")
  ATHz_kx_o = zeros(size(Akxo))

  ASH = copy(ATHz_kx_o)

  A_kompozit = compositeInput(Akxo, ATHz_kx_o, ASH)

  z[1] = 0

  plotInteraction::Bool = false
  #STR = Dates.format(now(), "YYYY-MM-DD hh-mm-ss")

  #global Axo_prew = zeros(size(Axo))
  FID = h5open(inputs.STR * ".hdf5", "w")
  entryCounter::Int = 1
  #STR = "elojel_minusz"
  #error()

  for ii in 1:(length(z)-1)
    A_kompozit, z[ii+1] = RK4M(thz_feedback_n2_SHG, z[ii], A_kompozit, inputs.dz, misc)
    if mod(ii, 11) == 0
      plotInteraction = true
    else
      plotInteraction = false
    end

    #if (mod(ii, 100) == 0 || ii == 1 ) && false
    if (ii == length(z) - 1 || mod(ii, 20) == 0 || ii == 1)

      Aop_kx_o = A_kompozit.Akxo
      Axo = FOPS.ifft_kx_x * ifftshift(Aop_kx_o, 2) .* kxMax .* exp.(-1im .* kx_omega .* cx - 1im .* kz_omega .* z[ii+1])
      Axt = FOPS.ifft_o_t * ifftshift(Axo .* omegaMax, 1)
      Aop_kx_oSH = A_kompozit.ASH
      AxoSH = FOPS.ifft_kx_x * ifftshift(Aop_kx_oSH, 2) .* kxMax .* exp.(-1im .* kx_omegaSHG .* cx - 1im .* kz_omegaSHG .* z[ii+1])
      AxtSH = FOPS.ifft_o_t * ifftshift(AxoSH .* omegaMax, 1)

      #p1 = heatmap(x, t, abs.(Axt .* exp.(1im .* omega0 .* ct) .* 1e-8), linewidth=0, colormap=:jet)
      #p3 = heatmap(x, omega, abs.(Axo - Axo_prew), colormap=:jet)
      #p4 = heatmap(x, t, abs.(AxtSH .* exp.(2im .* omega0 .* ct) .* 1e-8), linewidth=0, colormap=:jet)
      ATHz_kx_o = A_kompozit.ATHz_kx_o
      ATHz_xo = FOPS.ifft_kx_x * ifftshift(ATHz_kx_o .* exp.(-1im .* k_omegaTHz .* z[ii+1]), 2) .* kxMax
      ATHz_xt = FOPS.ifft_o_t * ATHz_xo * omegaMax
      #p2 = heatmap(x, t, real.(ATHz_xt) * 1e-5, linewidth=0, colormap=:jet)
      _, max_indices = findmax(abs.(Axt))
      #(scatter!([x[max_indices[2]]], [t[max_indices[1]]]))
      # display(plot(p1, p2, p3, p4, layout=(2, 2), size=[1200, 900]))
      #global Axo_prew = copy(Axo)
      #display(heatmap(x, t, abs.(ATHz_kx_o), linewidth=0, colormap=:jet))
              FID["/"*string(entryCounter)*"/Eop"] = Axt
              FID["/"*string(entryCounter)*"/Aop"] = Axo
              FID["/"*string(entryCounter)*"/ASH"] = AxoSH
              FID["/"*string(entryCounter)*"/ATHz_xo"] = ATHz_xo
              FID["/"*string(entryCounter)*"/ATHz_xt"] = ATHz_xt
           
      entryCounter += 1
    end
    display(ii)
  end
  FID["/maxEntry"] = entryCounter - 1
  FID["/gamma"] = rad2deg(inputs.gamma)
  FID["/z"] = z
  FID["/omega"] = omega
  FID["/omega0"] = omega0
  FID["/x"] = collect(inputs.x)
  FID["/t"] = collect(inputs.t)
  close(FID)
end
