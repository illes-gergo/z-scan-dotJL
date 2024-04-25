include("typedefs.jl")

function imp_terjedes(t, Y, PFC::pumpFieldConstants, RTC::runTimeConstants)
  dAdz = -1im .* PFC.kx_omega .* RTC.ckx ./ PFC.kz_omega .* Y .+ 1im .* RTC.ckx .^ 2 ./ 2 ./ PFC.kz_omega .* Y
  return dAdz
end

function imp_terjedesSH(t, Y, SFC::SHFieldConstants, RTC::runTimeConstants)
  dAdz = -1im .* SFC.kx_omegaSHG .* RTC.ckx ./ SFC.kz_omegaSHG .* Y + 1im .* RTC.ckx .^ 2 ./ 2 ./ SFC.kz_omegaSHG .* Y
  return dAdz
end

function thz_generation(t, Y, misc::miscInputs)
  Eop = misc.FOPS.ifft_kx_x * ifftshift(Y, 2) * misc.RTC.kxMax .* exp.(-1im .* misc.PFC.kx_omega .* misc.RTC.cx - 1im .* misc.PFC.kz_omega .* t)
  conv_part = fast_forward_convolution(Eop, conj(Eop), misc.RTC, misc.FOPS) * misc.NC.e0 * misc.RTC.khi_eff * misc.RTC.dOmega
  return fftshift(misc.FOPS.fft_x_kx * (conv_part), 2) ./ misc.RTC.kxMax
end

function thz_egyszeru(t, Y)
  Aop = Y[:, :, 1]
  ATHz = Y[:, :, 2]
  dAopdz = @spawn imp_terjedes(t, Aop)
  dTHz_gen = @spawn begin
    temp_val = -1im .* comegaTHz .^ 2 ./ 2 ./ kz_omegaTHz ./ e0 ./ c0 .^ 2 .* thz_generation(t, Aop) .* exp.(1im .* kz_omegaTHz .* t) - alpha / 2 .* ATHz
    temp_val[isnan.(temp_val)] .= 0
    return temp_val
  end
  wait.([dAopdz, dTHz_gen])
  return cat(dAopdz.result, dTHz_gen.result, zeros(size(Aop)), dims=3)
end

function plan_fast_conv(a, b, RT::runTimeConstants)
  a_ = vcat(RT.padding, a)
  b_ = vcat(b, RT.padding)
  _, fast_conv_plan = FourierTools.plan_conv(a_, b_, 1)
  fast_conv_fft_plan = plan_fft(b_, 1)
  return fast_conv_plan, fast_conv_fft_plan
end

function fast_forward_convolution(a, b, RT::runTimeConstants, FO::fourierOperations)
  a_ = vcat(RT.padding, a)
  b_ = vcat(b, RT.padding)
  return FO.fast_conv_plan(circshift(a_, (1, 0)), FO.fast_conv_fft_plan * (reverse((b_), dims=(1))))[floor(Int, end / 2)+1:end, :]
end

function fast_backward_convolution(a, b, RT::runTimeConstants, FO::fourierOperations)
  a_ = vcat(RT.padding, a)
  b_ = vcat(b, RT.padding)
  return reverse(FO.fast_conv_plan(circshift(reverse(a_, dims=1), (1, 0)), FO.fast_conv_fft_plan * (reverse((b_), dims=(1)))), dims=1)[floor(Int, end / 2)+1:end, :]
end

function fast_forward_convolution_SHG(a, b, RT::runTimeConstants, FO::fourierOperations)
  a_ = circshift(vcat(RT.padding, a), (RT.SHG_SHIFT, 0))
  b_ = circshift(vcat(b, RT.padding), (-RT.SHG_SHIFT, 0))
  return FO.fast_conv_plan(circshift(a_, (1, 0)), FO.fast_conv_fft_plan * (reverse((b_), dims=(1))))[floor(Int, end / 2)+1:end, :]
end

function fast_backward_convolution_SHG(a, b, RT::runTimeConstants, FO::fourierOperations)
  a_ = circshift(vcat(RT.padding, a), (-RT.SHG_SHIFT, 0))
  b_ = circshift(vcat(b, RT.padding), (-RT.SHG_SHIFT, 0))
  return reverse(FO.fast_conv_plan(circshift(reverse(a_, dims=1), (1, 0)), FO.fast_conv_fft_plan * (reverse((b_), dims=(1)))), dims=1)[floor(Int, end / 2)+1:end, :]
end

function thz_cascade(t, Aop, ATHz, misc::miscInputs)
  Eop = @spawn misc.FOPS.ifft_kx_x * ifftshift(Aop, 2) * misc.RTC.kxMax .* exp.(-1im .* misc.PFC.kx_omega .* misc.RTC.cx - 1im .* misc.PFC.kz_omega .* t)
  ETHz = @spawn misc.FOPS.ifft_kx_x * ifftshift(ATHz .* misc.RTC.kxMax .* exp.(-1im .* misc.TFC.kz_omegaTHz .* t), 2)
  wait.([Eop, ETHz])
  temp_val1 = @spawn misc.NC.e0 .* misc.RTC.khi_eff .* fast_forward_convolution(Eop.result, conj(ETHz.result), misc.RTC, misc.FOPS)
  temp_val2 = @spawn misc.NC.e0 .* misc.RTC.khi_eff .* fast_backward_convolution(Eop.result, ETHz.result, misc.RTC, misc.FOPS)
  wait.([temp_val1, temp_val2])

  return fftshift(misc.FOPS.fft_x_kx * ((temp_val1.result .+ temp_val2.result) .* exp.(+1im .* misc.PFC.kx_omega .* misc.RTC.cx)) ./ misc.RTC.kxMax .* misc.RTC.dOmega, 2)
end

function thz_feedback(t, Y)
  Aop = Y[:, :, 1]
  ATHz = Y[:, :, 2]
  dAop_lin = @spawn imp_terjedes(t, Aop)
  dTHz_gen = @spawn begin
    temp_val = -1im .* comegaTHz .^ 2 ./ 2 ./ kz_omegaTHz ./ e0 ./ c0 .^ 2 .* thz_generation(t, Aop) .* exp.(1im .* kz_omegaTHz .* t) - alpha / 2 .* ATHz
    temp_val[isnan.(temp_val)] .= 0
    return temp_val
  end
  dAopCsc = @spawn thz_cascade(t, Aop, ATHz) .* exp.(1im .* kz_omega .* t)

  wait.([dAop_lin, dTHz_gen, dAopCsc])

  return cat(dAop_lin.result .- dAopCsc.result .* 1im .* comega .^ 2 ./ 2 ./ kz_omega ./ e0 ./ c0 .^ 2, dTHz_gen.result, zeros(size(Aop)), dims=3)
end

function n2calc(t, Aop, misc::miscInputs)
  mult1 = @spawn -2 .* misc.PFC.kz_omega .* misc.NC.e0 .* misc.NC.c0 .^ 2 ./ misc.RTC.comega .^ 2
  Eop = @spawn misc.FOPS.ifft_o_t * ifftshift(misc.FOPS.ifft_kx_x * ifftshift(Aop, 2) * misc.RTC.kxMax .* exp.(-1im .* misc.PFC.kx_omega .* misc.RTC.cx - 1im .* misc.PFC.kz_omega .* t), 1) * misc.RTC.omegaMax
  mult2 = misc.NC.e0 * misc.RTC.omega0 * neo(misc.RTC.lambda0, 300, misc.RTC.cry) * n2value(misc.RTC.cry) / 2
  wait(Eop)
  aEop2 = abs.(Eop.result) .^ 2
  wait(mult1)
  return fftshift(misc.FOPS.fft_x_kx * (mult1.result .* fftshift(misc.FOPS.fft_t_o * (mult2 .* aEop2 .* Eop.result), 1) ./ misc.RTC.omegaMax .* exp.(+1im .* misc.PFC.kx_omega .* misc.RTC.cx) ./ misc.RTC.kxMax), 2)
end

function thz_feedback_n2(t, Y)
  Aop = Y[:, :, 1]
  ATHz = Y[:, :, 2]
  dAop_lin = @spawn imp_terjedes(t, Aop)
  dTHz_gen = @spawn begin
    temp_val = -1im .* comegaTHz .^ 2 ./ 2 ./ kz_omegaTHz ./ e0 ./ c0 .^ 2 .* thz_generation(t, Aop) .* exp.(1im .* kz_omegaTHz .* t) - alpha / 2 .* ATHz
    temp_val[isnan.(temp_val)] .= 0
    return temp_val
  end
  dAopCsc = @spawn thz_cascade(t, Aop, ATHz) #.* exp.(1im .* kz_omega .* t)
  dAopn2 = @spawn n2calc(t, Aop)

  sum_dAop = @spawn begin
    wait.([dAopCsc, dAopn2])
    return (dAopCsc.result .- dAopn2.result) .* exp.(1im .* kz_omega .* t)
  end

  wait.([sum_dAop, dTHz_gen, dAop_lin])
  return cat(dAop_lin.result .- sum_dAop.result .* 1im .* comega .^ 2 ./ 2 ./ kz_omega ./ e0 ./ c0 .^ 2, dTHz_gen.result, zeros(size(Aop)), dims=3)
end

function SHG_GEN(t, Aop, misc::miscInputs)
  Eop = misc.FOPS.ifft_kx_x * ifftshift(Aop, 2) * misc.RTC.kxMax .* exp.(-1im .* misc.PFC.kx_omega .* misc.RTC.cx - 1im .* misc.PFC.kz_omega .* t)
  temp_val = misc.NC.e0 .* misc.RTC.d_eff .* fast_backward_convolution_SHG(Eop, Eop, misc.RTC, misc.FOPS)
  return fftshift(misc.FOPS.fft_x_kx * (temp_val .* exp.(+1im .* misc.SFC.kx_omegaSHG .* misc.RTC.cx)) ./ misc.RTC.kxMax .* misc.RTC.dOmega, 2)
end

function SH_OP_INTERACTION(t, Aop, ASH, misc::miscInputs)
  Eop = @spawn misc.FOPS.ifft_kx_x * ifftshift(Aop, 2) * misc.RTC.kxMax .* exp.(-1im .* misc.PFC.kx_omega .* misc.RTC.cx - 1im .* misc.PFC.kz_omega .* t)
  ESH = @spawn misc.FOPS.ifft_kx_x * ifftshift(ASH, 2) * misc.RTC.kxMax .* exp.(-1im .* misc.SFC.kx_omegaSHG .* misc.RTC.cx - 1im .* misc.SFC.kz_omegaSHG .* t)
  wait.([Eop, ESH])
  conv_part = fast_forward_convolution_SHG(ESH.result, conj(Eop.result),misc.RTC,misc.FOPS) * misc.NC.e0 * misc.RTC.d_eff * misc.RTC.dOmega

  return fftshift(misc.FOPS.fft_x_kx * (conv_part .* exp.(+1im .* misc.PFC.kx_omega .* misc.RTC.cx)) ./ misc.RTC.kxMax, 2)
end

function thz_feedback_n2_SHG(t, Y::compositeInput, misc::miscInputs)
  Aop = Y.Akxo
  ATHz = Y.ATHz_kx_o
  ASH = Y.ASH
  dAop_lin = @spawn imp_terjedes(t, Aop, misc.PFC, misc.RTC)
  dTHz_gen = @spawn begin
    temp_val = -1im .* misc.RTC.comegaTHz .^ 2 ./ 2 ./ misc.TFC.kz_omegaTHz ./ misc.NC.e0 ./ misc.NC.c0 .^ 2 .* thz_generation(t, Aop, misc) .* exp.(1im .* misc.TFC.kz_omegaTHz .* t) - misc.TFC.alpha / 2 .* ATHz
    temp_val[isnan.(temp_val)] .= 0
    return temp_val
  end
  dAopCsc = @spawn thz_cascade(t, Aop, ATHz, misc) #=zeros(size(Aop))=#
  dAopn2 = @spawn n2calc(t, Aop, misc) #=zeros(size(Aop))=#
  dAopSH = @spawn SH_OP_INTERACTION(t, Aop, ASH, misc)
  dASHNL = @spawn SHG_GEN(t, Aop, misc)
  dASHlin = @spawn imp_terjedesSH(t, ASH, misc.SFC, misc.RTC)

  wait(dASHNL)
  dASHNL.result = dASHNL.result .* exp.(1im .* misc.SFC.kz_omegaSHG .* t)
  sum_dAop = @spawn begin
    wait.([dAopCsc, dAopn2, dAopSH])
    return (dAopCsc.result .- dAopn2.result .+ dAopSH.result .* 2) .* exp.(1im .* misc.PFC.kz_omega .* t)
  end

  wait.([sum_dAop, dTHz_gen, dAop_lin, dASHlin])
  return compositeInput(dAop_lin.result .- sum_dAop.result .* 1im .* misc.RTC.comega .^ 2 ./ 2 ./ misc.PFC.kz_omega ./ misc.NC.e0 ./ misc.NC.c0 .^ 2,
    dTHz_gen.result, dASHlin.result .- dASHNL.result .* 1im .* misc.RTC.comegaSHG .^ 2 ./ 2 ./ misc.SFC.kz_omegaSHG ./ misc.NC.e0 ./ misc.NC.c0 .^ 2)
end
