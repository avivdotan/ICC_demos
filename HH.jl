### A Pluto.jl notebook ###
# v0.12.3

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

# ╔═╡ b63cade0-22bf-11eb-1847-bb0e34bcaeaf
begin
	#Imports
	using PlutoUI
	using Latexify
	using DifferentialEquations
	using Plots
end

# ╔═╡ 20e11820-22c0-11eb-36e6-1716f7dd2a93
"""
Calculate the parameters for n,h,m differential equations
for details, see e.g.:
www.math.pitt.edu/~bdoiron/assets/ermentrout-and-terman-ch-1.pdf
"""
function gates_params(V)

	V₁ = V + 60

	αₙ = (10 - V₁)/(100(exp((10 - V₁)/10) - 1))
	βₙ = 0.125exp(-V₁/80)
	n∞ = αₙ/(αₙ + βₙ)
	τₙ = 1/(αₙ + βₙ)

	αₘ = (25 - V₁)/(10(exp((25 - V₁)/10) - 1))
	βₘ = 4exp(-V₁/18)
	m∞ = αₘ/(αₘ + βₘ)
	τₘ = 1/(αₘ + βₘ)

	αₕ = 0.07exp(-V₁/20)
	βₕ = 1/(exp((30 - V₁)/10) + 1)
	h∞ = αₕ/(αₕ + βₕ)
	τₕ = 1/(αₕ + βₕ)
	
	return (n∞, τₙ, m∞, τₘ, h∞, τₕ)
end;

# ╔═╡ 9a409100-22c0-11eb-0fb3-e9792f351784
"""
HH Hodgkin-Huxley Model
"""
function HH!(du, u, p, t)

	V, n, m, h = u
	g_Na, V_Na, g_K, V_K, g_l, V_l, C, I = p
	n∞, τₙ, m∞, τₘ, h∞, τₕ = gates_params(V)

	Iᵢₒₙ = g_Na*h*m^3 * (V - V_Na) + 	# I_Na
		   g_K*n^4 * (V - V_K) + 		# I_K
		   g_l * (V - V_l) 				# I_l
	
	du[1] = (-Iᵢₒₙ + I(t))/C 			# dV/dt = ...
	du[2] = (n∞ - n)/τₙ 				# dn/dt = ...
	du[3] = (m∞ - m)/τₘ 				# dm/dt = ...
	du[4] = (h∞ - h)/τₕ 				# dh/dt = ...

end;

# ╔═╡ ef25b370-22c1-11eb-018c-05ae50bd2c15
md"""
## Hodgkin–Huxley Model -- A Simulation

```math
C\frac{\mathrm{d}V}{\mathrm{d}t} = -\left[
\bar{g}_{\mathrm{Na}}m^3h\left(V - V_{\mathrm{Na}}\right) + 
\bar{g}_{\mathrm{K}}n^4\left(V - V_{\mathrm{K}}\right) + 
\bar{g}_{L}n^4\left(V - V_{L}\right)
\right] + 
I(t)
```

```math
\frac{\mathrm{d}n}{\mathrm{d}t} = \frac{n_{\infty}\left(V\right) - n}{\tau_n\left(V\right)}
\qquad
\frac{\mathrm{d}m}{\mathrm{d}t} = \frac{m_{\infty}\left(V\right) - m}{\tau_m\left(V\right)}
\qquad
\frac{\mathrm{d}h}{\mathrm{d}t} = \frac{h_{\infty}\left(V\right) - h}{\tau_h\left(V\right)}
```

Current shape $(@bind Iₛ Select(["none" => "None", "step" => "Step", "sin" => "Sinusoidal", "rand" => "Random"], default = "step"))
"""

# ╔═╡ a36a7a40-22c8-11eb-18a3-93d9128b3017
begin
	
	I₀_min = -10
	I₀_max = 150
	ΔI₀ = 1
	
	# Set LaTeX display
	set_default(cdot = false, fmt=FancyNumberFormatter(3))
	
	if Iₛ == "none"
		# I(t) = 0
		I⁼ = :(0)
	
		md"""
		``I(t) = `` $(latexify(I⁼))
		"""
	elseif Iₛ == "step"
		# I(t) = 50 < t < 150 ? I₀ : 0
		I⁼ = :(20*ms < t < 150*ms ? I₀ : 0)
		const ms = 1
	
		md"""
		``I(t) = `` $(latexify(I⁼))
		
		``I_0`` $(@bind I₀ Slider(I₀_min:ΔI₀:I₀_max, default = 15, show_value = true)) [μA/cm²]
		"""
	elseif Iₛ == "sin"
		# f = 30 			#[Hz]
		# ω = 2π*f/1000 	#[rad/sec]
		# I(t) = I₀*(1 + cos(ω*t))/2
		
		I⁼ = :(I₀ + I₁*cos(ω*t))
	
		md"""
		``I(t) = `` $(latexify(I⁼))
		
		``I_0`` $(@bind I₀ Slider(0:ΔI₀:(I₀_max/2), default = 30, show_value = true)) [μA/cm²]
		
		``I_1`` $(@bind I₁ Slider(0:ΔI₀:(I₀_max/2), default = 30, show_value = true)) [μA/cm²]
		
		``\omega`` $(@bind ω Slider(0:0.001:2π*0.1, default = 0.18, show_value = true)) [rad/sec]
		"""
	
	elseif Iₛ == "rand"
		# I(t) = I₀*rand()
		I⁼ = :(I₀*rand())
	
		md"""
		``I(t) = I_0\cdot \mathrm{rand()}``
		
		``I_0`` $(@bind I₀ Slider(I₀_min:ΔI₀:I₀_max, default = 15, show_value = true)) [μA/cm²]
		"""
	end
end

# ╔═╡ 2a82ed52-22c4-11eb-34ea-c14874528737
begin
	
	# Define the current
	I⁼
	if Iₛ == "step"
		I₀
	elseif Iₛ == "sin"
		I₀, I₁, ω
	elseif Iₛ == "rand"
		I₀
	end
	eval(:(I(t) = $I⁼))
	
	# Set initial conditions
	V₀ = -60
	n₀, ~, m₀, ~, h₀, ~ = gates_params(V₀)
	u₀ = [V₀, n₀, m₀, h₀]
	
	# Set solution parameters
	tspan = (0.0, 200.0)
	p = [g_Na, V_Na, g_K, V_K, g_l, V_l, C, I]
	
	# Solve the equations system
	prob = ODEProblem(HH!, u₀, tspan, p)
	sol = solve(prob)
	
end;

# ╔═╡ Cell order:
# ╟─b63cade0-22bf-11eb-1847-bb0e34bcaeaf
# ╟─20e11820-22c0-11eb-36e6-1716f7dd2a93
# ╟─9a409100-22c0-11eb-0fb3-e9792f351784
# ╟─ef25b370-22c1-11eb-018c-05ae50bd2c15
# ╟─a36a7a40-22c8-11eb-18a3-93d9128b3017
# ╟─2a82ed52-22c4-11eb-34ea-c14874528737
