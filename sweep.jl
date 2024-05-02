include("main.jl")

function sweep()
  zR = -10:0.5:10
  transmission = zeros(size(zR))

  @threads for i in eachindex(zR)
    transmission[i] = runcalc(zR=zR[i])
  end

  p = plot(zR, transmission)
  savefig(p, "result.png", width=640, height=480, scale=5)
  display(p)
  return nothing
end

