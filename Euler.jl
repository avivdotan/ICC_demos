### A Pluto.jl notebook ###
# v0.12.10

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

# ╔═╡ 11fa4afe-131c-11eb-0207-1587f2b52fd4
begin
	using PlutoUI
	
	"""
	Log-scale Slider
	"""
	struct LogSlider
		range::AbstractRange
		default::Number
		show_value::Bool
	end

	LogSlider(range::AbstractRange; default=missing, show_value=false) = LogSlider(range, (default === missing) ? first(range) : default, show_value)

	function Base.show(io::IO, ::MIME"text/html", slider::LogSlider)
		print(io, """<input 
			type="range" 
			min="$(first(slider.range))" 
			step="$(step(slider.range))" 
			max="$(last(slider.range))" 
			value="$(slider.default)"
			$(slider.show_value ? "oninput=\"this.nextElementSibling.value=Math.pow(10,this.value).toPrecision(1)\"" : "")
			>""")

		if slider.show_value
			print(io, """<output>$(10.0^slider.default)</output>""")
		end
	end

	Base.get(slider::LogSlider) = 10.0^slider.default
end;

# ╔═╡ 2c685200-1314-11eb-09e6-9d0936417ada
begin
	τ_min = 1
	τ_max = 30
	Δτ = 1
	
	y∞_min = -30
	y∞_max = 30
	Δy∞ = 1
	
	md"""
	### Euler method demo
	
	
	Time constant ``\tau`` $(@bind τ Slider(τ_min:Δτ:τ_max, default = 3, show_value = true))

	Terminal value ``y_{\infty}`` $(@bind y∞ Slider(y∞_min:Δy∞:y∞_max, default = 18, show_value = true))
	"""
end

# ╔═╡ 0606b130-204b-11eb-09ba-47c21dee5dbc
begin
	
	# Set cell dependencies
	τ
	y∞
	
	# The differential equation is dy/dx = f(x, y) = ...
	f⁼ = :(-(y/$τ) + $(y∞/τ))
	
	# Analytical solution yₜ(x) = ...
	yₜ⁼ = :(y₀*exp(-(x/$τ)) + $y∞ * (1 - exp(-(x/$τ))))
	
	# Set display of equations
	using Latexify
	set_default(cdot = false, fmt=FancyNumberFormatter(3))
	
	y₀_min = -30
	y₀_max = 30
	Δy₀ = 1
	
	logΔx_min = -2
	logΔx_max = 1
	ΔlogΔx = 0.1
	
	md"""
	The differential equation is: ``\frac{\mathrm{d}x}{\mathrm{d}y} = `` $(latexify(f⁼))
	
	The analytical solution is: ``y_t \left(x\right) = `` $(latexify(yₜ⁼))
	
	Initial condition ``y_0`` $(@bind y₀ Slider(y₀_min:Δy₀:y₀_max, default = 0, show_value = true))

	Step size ``\Delta x`` $(@bind logΔx LogSlider(logΔx_min:ΔlogΔx:logΔx_max, default = 0, show_value = true))
	"""
end

# ╔═╡ 7922ce90-1314-11eb-0077-fd83044eba3c
begin
	# Build functions for the equation and the analytical solution
	eval(:(f(x, y) = $f⁼))
	eval(:(yₜ(x) = $yₜ⁼))
	
	# Simulation parameters
	Δx = 10.0^logΔx
	xᶠ = 50
	x = 0:Δx:xᶠ
	y = fill(NaN, size(x))
	y[1] = y₀
end;

# ╔═╡ 9c183d40-1314-11eb-0b66-5f4d06a8b721
begin
	# Plots
	using Plots
	
	xₜ = 0:0.001:xᶠ
	plot(xₜ, yₜ.(xₜ), label = "Analytical solution")
	plot!(x, y, label = "Numerical approximation")
	xlims!(0, xᶠ)
	ylims!(minimum([y₀_min, y∞_min, minimum(y)]), maximum([y₀_max, y∞_max, maximum(y)]))
	xlabel!("x")
	ylabel!("y")
end

# ╔═╡ 8a45b570-1314-11eb-03ab-0111c5a17e22
begin
	# Euler method
	for n = 1:(length(y) - 1)
		y[n + 1] = y[n] + f(x[n], y[n])*Δx
	end
end;

# ╔═╡ Cell order:
# ╠═11fa4afe-131c-11eb-0207-1587f2b52fd4
# ╠═2c685200-1314-11eb-09e6-9d0936417ada
# ╟─0606b130-204b-11eb-09ba-47c21dee5dbc
# ╟─7922ce90-1314-11eb-0077-fd83044eba3c
# ╟─8a45b570-1314-11eb-03ab-0111c5a17e22
# ╟─9c183d40-1314-11eb-0b66-5f4d06a8b721
