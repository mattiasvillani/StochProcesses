using Plots, GaussianProcesses, KernelFunctions, LaTeXStrings, Distributions
using Utils

figFolder = "/home/mv/Dropbox/Teaching/StochProcesses/Slides/figs/"
gifFolder = "/home/mv/Dropbox/Teaching/StochProcesses/GIFs/"

Plots.reset_defaults()
gr(legend = nothing, grid = false, color = colors[2], lw = 2, legendfontsize=10,
    xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, 
    markersize = 4, markerstrokecolor = :auto)


harmonic(t, f, A, ϕ) = A*cos(2π*f*t + ϕ)
tGrid = -5:0.01:5
f =.5; A = 2; ϕ = 0
plot(tGrid, harmonic.(tGrid, f, A, ϕ), xlab = L"t", 
    label = L"A = %$A, \phi = %$ϕ, f_0 = %$f", legend = :topright, 
    c = colors[10], title = L"\cos(2\pi f_0 t + \phi)")
f =.5; A = 2; ϕ = 1
plot!(tGrid, harmonic.(tGrid, f, A, ϕ), xlab = L"t", 
    label = L"A = %$A, \phi = %$ϕ, f_0 = %$f", c = colors[9])
f =.5; A = 1; ϕ = 1
plot!(tGrid, harmonic.(tGrid, f, A, ϕ), xlab = L"t", 
    label = L"A = %$A, \phi = %$ϕ, f_0 = %$f", c = colors[2])
f = 1.0; A = 1; ϕ = 1
plot!(tGrid, harmonic.(tGrid, f, A, ϕ), xlab = L"t", 
    label = L"A = %$A, \phi = %$ϕ, f_0 = %$f", c = colors[1])
savefig(figFolder*"singleharmonics.pdf")

# Animate stochastic of a single harmonic
f = 1
anim = @animate for i in 1:50
    A = abs(round(randn(), digits = 4)); ϕ = round(2π*rand(), digits = 4);
    plot(tGrid, harmonic.(tGrid, f, A, ϕ), xlab = L"t", 
    title = L"A = %$A, \phi = %$ϕ, f_0 = %$f", legend = nothing, 
    c = colors[2], ylab = L"\cos(2\pi f_0 t + \phi)")
    ylims!((-3,3))
    xlims!((-5,5))
end 
gif(anim, gifFolder*"randomharmonics.gif", fps = 1)

# Plot sum of harmonics
p = plot(xlab = L"t", legend = :topright)
harmonicSum = zeros(length(tGrid))
for i = 1:5
    A = abs(round(randn(), digits = 4))
    ϕ = round(2π*rand(), digits = 4)
    f = round(-0.5 + rand(), digits = 4)
    harmonicSum .+= harmonic.(tGrid, f, A, ϕ)
    if i == 1
        plot!(p, tGrid, harmonic.(tGrid, f, A, ϕ), 
        xlab = L"t", label = "random harmonics", c = colors[1])
    else
        plot!(p, tGrid, harmonic.(tGrid, f, A, ϕ), 
            xlab = L"t", label = nothing, c = colors[1])
    end
    ylims!((-5,5))
    xlims!((-5,5))
    if i == 5
        plot!(p, tGrid, harmonicSum, c = colors[4], label = "sum of harmonics")
    end
end
p
savefig(figFolder*"sumharmonics.pdf")


anim = @animate for i = 1:(5*5)
    if rem(i,5) == 1
        global p = plot(xlab = L"t", legend = :topright)
        global harmonicSum = zeros(length(tGrid))
    end
    A = abs(round(randn(), digits = 4))
    ϕ = round(2π*rand(), digits = 4)
    f = round(-0.5 + rand(), digits = 4)
    harmonicSum .+= harmonic.(tGrid, f, A, ϕ)
    if rem(i,5) == 1
        plot!(p, tGrid, harmonic.(tGrid, f, A, ϕ), 
        xlab = L"t", label = "random harmonics", c = colors[1])
    else
        plot!(p, tGrid, harmonic.(tGrid, f, A, ϕ), 
            xlab = L"t", label = nothing, c = colors[1])
    end
    ylims!((-5,5))
    xlims!((-5,5))
    if rem(i,5) == 0
        plot!(p, tGrid, harmonicSum, c = colors[4], label = "sum of harmonics")
    end
end
gif(anim, gifFolder*"randommultipleharmonics.gif", fps = 1)


function simAR1(T, ϕ, σ, μ)
    x = zeros(T)
    x[1] = μ[1]
    for t = 2:T
        x[t] = μ[t] + ϕ*(x[t-1] - μ[t]) + σ*randn()
    end
    return x
end

# Ensemble mean - non constant mean
gr(legendfontsize=12)
nRep = 10
T = 10 + 50
σ = 0.1
μ(t, T) = 0 .+ 1*(t/T) .- 1*(t/T)^2
μvect = μ.(1:T, T)

# iid
ϕ = 0.0
X = zeros(nRep, T)
for rep = 1:nRep
    X[rep,:] = simAR1(T, ϕ, σ, μvect)
end
X = X[:,11:end]
p = plot(1:50, X[1,:], c = colors[1], lw = 1, legend = :topright, label = "realization 1", 
    title = L"\mathrm{iid}", xlab = L"t")
plot!(p, 1:50, X[2,:], c = colors[3], lw = 1, label = "realization 2")
for t in 1:50
    scatter!(p, repeat([t],nRep), X[:,t], c = :lightgray, markersize = 3, label = nothing)
end
plot!(1:50, μvect[11:end], c = colors[10], label = L"m(t)")
plot!(p, 1:50, mean(X, dims = 1)[:], c = :black, lw = 3, label = L"\hat m_n(t)")
savefig(figFolder*"ensemble_iid.pdf")

ϕ = 0.9
X = zeros(nRep, T)
for rep = 1:nRep
    X[rep,:] = simAR1(T, ϕ, σ, μvect)
end
X = X[:,11:end]
p = plot(1:50, X[1,:], c = colors[1], lw = 1, legend = :topright, label = "realization 1", 
    title = L"\mathrm{AR}(\phi = 0.9)", xlab = L"t")
plot!(p, 1:50, X[2,:], c = colors[3], lw = 1, label = "realization 2")
for t in 1:50
    scatter!(p, repeat([t],nRep), X[:,t], c = :lightgray, markersize = 3, label = nothing)
end
plot!(1:50, μvect[11:end], c = colors[10], label = L"m(t)")
plot!(p, 1:50, mean(X, dims = 1)[:], c = :black, lw = 3, label = L"\hat m_n(t)")
savefig(figFolder*"ensemble_AR09.pdf")
p


# Ensemble mean - constant mean
gr(legendfontsize=12)
nRep = 10
T = 10 + 50
σ = 0.1
μvect = repeat([2], T)

# iid
ϕ = 0.0
X = zeros(nRep, T)
for rep = 1:nRep
    X[rep,:] = simAR1(T, ϕ, σ, μvect)
end
X = X[:,11:end]
p = plot(1:50, X[1,:], c = colors[1], lw = 1, legend = :topright, label = "realization 1", 
    title = L"\mathrm{iid}", xlab = L"t")
plot!(p, 1:50, X[2,:], c = colors[3], lw = 1, label = "realization 2")
for t in 1:50
    scatter!(p, repeat([t],nRep), X[:,t], c = :lightgray, markersize = 3, label = nothing)
end
plot!(1:50, μvect[11:end], c = colors[10], label = L"m(t)")
plot!(p, 1:50, mean(X, dims = 1)[:], c = :black, lw = 3, label = L"\hat m_n(t)")
savefig(figFolder*"ensemble_iid_stationary.pdf")

ϕ = 0.9
X = zeros(nRep, T)
for rep = 1:nRep
    X[rep,:] = simAR1(T, ϕ, σ, μvect)
end
X = X[:,11:end]
p = plot(1:50, X[1,:], c = colors[1], lw = 1, legend = :topright, label = "realization 1", 
    title = L"\mathrm{AR}(\phi = 0.9)", xlab = L"t")
plot!(p, 1:50, X[2,:], c = colors[3], lw = 1, label = "realization 2")
for t in 1:50
    scatter!(p, repeat([t],nRep), X[:,t], c = :lightgray, markersize = 3, label = nothing)
end
plot!(1:50, μvect[11:end], c = colors[10], label = L"m(t)")
plot!(p, 1:50, mean(X, dims = 1)[:], c = :black, lw = 3, label = L"\hat m_n(t)")
savefig(figFolder*"ensemble_AR09_stationary.pdf")
p


# Sinc 
r_sinc(τ) = sin(π*τ)/(π*τ)

plot(-10:0.01:10, r_sinc.(-10:0.01:10), xlab = "time lag, "*L"\tau", ylab = L"r(\tau)", 
    title = L"\mathrm{sinc}(\tau) = \sin(\pi\tau)/\pi\tau")
savefig(figFolder*"ACF_sinc.pdf")

box(f) = abs(f) <=0.5 ? 1 : 0
plot(-1:0.01:1, box.(-1:0.01:1), xlab = "frequency, "*L"f", ylab = L"R(f)", 
    title = "Spectral density for ACF: "*L"r(\tau) = \sin(\pi\tau)/\pi\tau", 
    ylim = [0,1.5], color = colors[4])
savefig(figFolder*"specdens_box.pdf")

r_sqrexp(τ, ℓ) = exp(-τ^2/(2*ℓ^2))
R_sqrexp(f, ℓ) = sqrt(2*π*ℓ^2)*exp(-2*π^2*f^2*ℓ^2)
ℓ = 1
plot(-5:0.01:5, r_sqrexp.(-5:0.01:5, ℓ), xlab = "time lag, "*L"\tau", 
    ylab = L"r(\tau)", 
    title = "squared exponential: "*L"r(\tau) = \sqrt{2 \pi \ell^2}\exp(-2 \pi^2 \ell^2 \tau^2)")
savefig(figFolder*"ACF_sqrexp.pdf")

plot(-1:0.01:1, R_sqrexp.(-1:0.01:1, ℓ), xlab = "frequency, "*L"f", ylab = L"R(f)", 
    title = "Spectral density for squared exponential ACF: ", 
    color = colors[4])
savefig(figFolder*"specdens_sqrexp.pdf")



ℓ = 1
k₁ = SqExponentialKernel()
τ_grid = 0:0.01:10
sqrexp_acf = kernelmatrix(k₁, τ_grid/ℓ, [0])
plot(τ_grid, sqrexp_acf, xlab = "time lag, "*L"\tau", ylab = L"r(\tau)",
    title =  "Squared exponential, "*L"\ell = 1")

# Function to plot draws from GP
function PlotGPDraws(nSamples, kernel, σ, ℓ; cols = colors[1:nSamples], 
        x_grid = 0:0.01:10, jitter = 1e-7)
    K₁ = kernelmatrix(kernel, x_grid/ℓ, x_grid/ℓ) + jitter*I(length(x_grid))
    x = rand(MvNormal(K₁), nSamples)
    p = plot(τ_grid, x, c = cols, xlab = "time, "*L"t", ylab = "Function, "*L"f(t)")
    return p
end

# Simulate from SqExponentialKernel() 
ℓ = 1
p = PlotGPDraws(3, SqExponentialKernel(), 1, ℓ; cols = colors[[2 4 6]])
title!(p, "Squared exponential, "*L"\ell = %$ℓ")
savefig(figFolder*"sim_sqrexp.pdf")

# Simulate from MaternKernel(ν = 1/2)
ℓ = 1
ν = 1/2
p = PlotGPDraws(3, MaternKernel(ν = ν), 1, ℓ; cols = colors[[2 4 6]])
title!(p, "Ornstein-Uhlenbeck - Matern "*L"(\nu=1/2)"*", "*L"\ell = %$ℓ")
savefig(figFolder*"sim_matern12.pdf")

# Simulate from MaternKernel(ν = 3/2)
ℓ = 1
ν = 3/2
p = PlotGPDraws(3, MaternKernel(ν = ν), 1, ℓ; cols = colors[[2 4 6]])
title!(p, "Matern "*L"(\nu=3/2)"*", "*L"\ell = %$ℓ")
savefig(figFolder*"sim_matern32.pdf")

# Simulate from MaternKernel(ν = 5/2)
ℓ = 1
ν = 5/2
p = PlotGPDraws(3, MaternKernel(ν = ν), 1, ℓ; cols = colors[[2 4 6]])
title!(p, "Matern "*L"(\nu=5/2)"*", "*L"\ell = %$ℓ")
savefig(figFolder*"sim_matern52.pdf")
