### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ fe8d2920-3fd5-11eb-11fc-4f99600b2f63
begin
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add("Distributions")
	Pkg.add(Pkg.PackageSpec(url="https://github.com/bci4cpl/Tempotrons.jl"))
	Pkg.add("Latexify")
	Pkg.add("Plots")
	using PlutoUI
	using Distributions
	using Tempotrons
	using Latexify
	using Plots
	
	# Latexify default
	set_default(cdot = false, fmt = FancyNumberFormatter(3))
	
	md"""
	# Synaptic Current Integration - Demo
	"""
end

# ╔═╡ 8899f7d0-3fe3-11eb-0e51-1171553266e6
md"""
```math
I_s\left(t\right) = I_0 e^{\frac{-t}{\tau_s}} Θ(t)
\quad,\qquad
K\left(t\right) = \frac{R}{\tau_m}\int_0^t e^{-\frac{t - t'}{\tau_m}}I_s\left(t'\right)\mathrm{d}t' = \frac{I_0 R}{\alpha - 1}\left(e^{-\frac{t}{\tau_m}} - e^{-\frac{t}{\tau_s}}\right)\Theta\left(t\right)
```
```math
\alpha \equiv \frac{\tau_m}{\tau_s}
\quad,\qquad
I_0 = 1\mathrm{mV}\cdot \frac{\alpha^{\frac{\alpha}{\alpha - 1}}}{R}
```
```math
I\left(t\right) = \sum_{i=1}^N w_i \sum_{t_i^j} I_s\left(t - t_i^j\right)
\quad,\qquad
V\left(t\right) = \sum_{i=1}^N w_i \sum_{t_i^j} K\left(t - t_i^j\right)
```
"""

# ╔═╡ 93e6f4b0-3fd6-11eb-19fb-65186c7c7eb1
md"""
Membrane time constant ``\tau_m`` $(@bind τₘ Slider(10:0.5:20, default = 15, show_value = true))ms
"""

# ╔═╡ bbe0b872-3fd6-11eb-25a2-8bf3eadd76c0
md"""
Synaptic time constant ``\tau_s`` $(@bind τₛ Slider(2:0.1:8, default = τₘ/4, show_value = true))ms
"""

# ╔═╡ c5f856b0-3fd6-11eb-1ddb-c1d6fa2de6d8
md"""
 $(@bind redraw_sample Button("Redraw sample"))
$(@bind redraw_weights Button("Redraw weights"))
"""

# ╔═╡ 33009700-3fd6-11eb-3c35-79bd39086ce6
begin
	# Set parameters
	
	N = 5
	ν = 10 	# Hz
	T = 200 # msec
	C = 1 	# μF/cm²
	θ = 30 	# mV
end;

# ╔═╡ 44e857ae-3fda-11eb-0e5c-6b7d11002ac3
begin
	# Plot functions
	
	function plot_input(smp; kwargs...)
		return plot([1, 1] * vcat(smp...)', vcat([i*ones(length(smp[i])) for i = 1:N]...)' .+ [-0.45, 0.45];
			linewidth = 2,
			xlim = (0, T), ylim = (0.5, N + 0.5),
			yticks = 1:N,
			xlabel = "t [ms]", ylabel = "Neuron. No.",
			label = "",
			kwargs...)
	end
	
end;

# ╔═╡ 6e20b6d0-43c8-11eb-345d-d5c31ee75824
begin
	tmp = Tempotron(N, τₘ = τₘ, τₛ = τₛ, θ = θ, V₀ = 0.0)
end;

# ╔═╡ de778c60-3fd6-11eb-3271-090e8ee8aeff
begin
	redraw_sample
	
	sample = [rand(Uniform(0,T), rand(Poisson(ν*T/1000))) for i = 1:N]
	
	n_spikes_max = sum(length.(sample))
	
	spikes = sort([(tᵢʲ = sample[i][j], i = i, j = j) for i = 1:length(sample) for j = 1:length(sample[i])], by = s -> s.tᵢʲ)
	

end;

# ╔═╡ 972089f0-3fd8-11eb-2be5-bf2c2e26f879
md"""
Number of input spikes: $(@bind n_spikes Slider(0:n_spikes_max, default = n_spikes_max, show_value = true))
"""

# ╔═╡ c9ce2bf0-3fd8-11eb-007d-1bee38f0d5e8
begin
	s_sample = n_spikes == 0 ? [[] for i = 1:N] : filter.(t -> t ≤ spikes[n_spikes].tᵢʲ, sample)
	s_sample = SpikesInput(s_sample)
end;

# ╔═╡ 73fe8270-3fd7-11eb-1e04-2bd3de2777d4
begin
	redraw_weights
	
	W = θ*(1.2rand(Float64, N) .- 0.2)
	tmp.w .= W
end;

# ╔═╡ 29976ef2-3fdb-11eb-3ba8-8d53a40da2d7
begin
	W
	latexify(:(w = $(W')))
end

# ╔═╡ efc3f440-3fd6-11eb-3f88-c7dd655770cc
begin
	# Set dependent parameters
	
	R = τₘ/C
	α = τₘ/τₛ
	I₀ = α^(α/(α - 1))/R
	Iₛ(t) = t < 0 ? 0.0 : I₀*exp(-t/τₛ)
	κ = I₀*R/(α - 1)
	K(t) = t < 0 ? 0.0 : κ*(exp(-t/τₘ) - exp(-t/τₛ))
end;

# ╔═╡ aa131d2e-3fd7-11eb-264c-e1391eca9a29
begin
	# Get current and voltage
	
	Δt = 0.1 		#ms
	time = 0:Δt:T 	#ms
	
	I(t) = isempty(s_sample) ? 0.0 : sum([W[i]*Iₛ(t - tᵢʲ) for i = 1:N for tᵢʲ ∈ s_sample[i]])
	
	V(t) = isempty(s_sample) ? 0.0 : sum([psp.ΔV(t) for psp ∈ Tempotrons.get_psps(tmp, s_sample)])
	
	Vₙ = tmp(s_sample, t = collect(time)).V
end;

# ╔═╡ f9788390-3fd9-11eb-32c6-c96b30195a56
begin
	# Plots
	gr(size = (675, 550))
	
	p_inp = plot_input(s_sample, 
		color = 1, 
		title = "Input Spikes")
	
	p_cur = plot(time, I.(time), 
		xlim = (0, T),
		xlabel = "t [ms]", ylabel = "I [μA/cm²]",
		color = 2,
		title = "Total Current",
		label = "")
	
	p_vol = plot(time, V.(time), 
		xlim = (0, T),
		xlabel = "t [ms]", ylabel = "V [mV]",
		color = 3,
		title = "Voltage (unresetted)",
		label = "")
	
	p_vol_n = plot(time, Vₙ, 
		xlim = (0, T),
		xlabel = "t [ms]", ylabel = "V [mV]",
		color = 4,
		title = "Voltage (resetted)",
		label = "")
	plot!([0, T], θ*[1, 1],
		color = :black, linestyle = :dash,
		label = "")
	
	plot(p_inp, p_cur, p_vol, p_vol_n, layout = (4, 1), link = :x)
end

# ╔═╡ Cell order:
# ╟─fe8d2920-3fd5-11eb-11fc-4f99600b2f63
# ╟─8899f7d0-3fe3-11eb-0e51-1171553266e6
# ╟─93e6f4b0-3fd6-11eb-19fb-65186c7c7eb1
# ╟─bbe0b872-3fd6-11eb-25a2-8bf3eadd76c0
# ╟─c5f856b0-3fd6-11eb-1ddb-c1d6fa2de6d8
# ╟─29976ef2-3fdb-11eb-3ba8-8d53a40da2d7
# ╟─972089f0-3fd8-11eb-2be5-bf2c2e26f879
# ╟─f9788390-3fd9-11eb-32c6-c96b30195a56
# ╟─44e857ae-3fda-11eb-0e5c-6b7d11002ac3
# ╟─33009700-3fd6-11eb-3c35-79bd39086ce6
# ╟─6e20b6d0-43c8-11eb-345d-d5c31ee75824
# ╟─de778c60-3fd6-11eb-3271-090e8ee8aeff
# ╟─c9ce2bf0-3fd8-11eb-007d-1bee38f0d5e8
# ╟─73fe8270-3fd7-11eb-1e04-2bd3de2777d4
# ╟─efc3f440-3fd6-11eb-3f88-c7dd655770cc
# ╟─aa131d2e-3fd7-11eb-264c-e1391eca9a29
