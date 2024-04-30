include("main.jl")

function sweep()
  zR = -5:0.25:5
  transmission = zeros(size(zR))

  @threads for i in eachindex(zR)
    transmission[i] = runcalc(zR=zR[i])
  end

  plot(zR, transmission)
end

sweep()
