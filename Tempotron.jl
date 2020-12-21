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

# ╔═╡ a2fffdb0-3ca8-11eb-0a28-7bb88a2f2292
begin
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add(Pkg.PackageSpec(url="https://github.com/bci4cpl/Tempotrons.jl"))
	Pkg.add("Distributions")
	Pkg.add("LatexPrint")
	Pkg.add("Latexify")
	Pkg.add("Plots")
	using PlutoUI
	using Tempotrons
	using Distributions
	using LatexPrint
	using Latexify
	using Plots
	using Plots.PlotMeasures
	
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
			$(slider.show_value ? "oninput=\"this.nextElementSibling.value=Math.pow(10,this.value).toPrecision(3)\"" : "")
			>""")

		if slider.show_value
			print(io, """<output>$(round(10.0^slider.default, sigdigits = 3))</output>""")
		end
	end

	Base.get(slider::LogSlider) = 10.0^slider.default
	
	md"""
	# Tempotron - Step by Step Demo
	"""
end

# ╔═╡ 16b92fa0-3cbe-11eb-129d-456f9aeebd93
md"""
 $(@bind redraw_samples Button("Redraw samples"))
$(@bind redraw_weights Button("Redraw weights"))
\
``\eta`` $(@bind log_η LogSlider(-1:0.1:1, default = 0.3, show_value = true))

$(@bind redraw_sample Button("Step"))
"""

# ╔═╡ d9ca7af0-3ca8-11eb-1d71-b571a6089f88
begin
	
	# Set time vector
	dt      = 0.1        	#[ms] time step for numerical integration
	t_final = 500           #[ms] duration of numerical simulation
	t       = 0:dt:t_final  #[ms] time vector
	
end;

# ╔═╡ 37e9d3f0-3caa-11eb-33b2-79045f34aa25
begin
	# Generate Poisson samples
	redraw_samples
	
	N = 5 			# number of input neurons
	n_samples = 20 	# number of samples
	ν = 10 			#[Hz] input firing rate
	
	Samples = [(x = InputGen.poisson_spikes_input(N, ν = ν, T = t_final), 
		  		y = Bool(2(s - 1)÷n_samples))
		 for s = 1:n_samples]
	
end;

# ╔═╡ 115700f0-3caf-11eb-287d-15083df59991
begin
	# Define weights
	redraw_weights
	
	tmp = Tempotron(N, θ = 30, V₀ = 0)
	tmp.w .= (tmp.θ - tmp.V₀)*(0.9rand(Float64, size(tmp.w)) .- 0.1)
end;

# ╔═╡ 5e7bb330-3cbe-11eb-0972-35792137cf8d
begin
	# Set plot parameters
	max_V = 2(tmp.θ - tmp.V₀)
end;

# ╔═╡ 2148f98e-3f7b-11eb-031f-9d22f9d77b47
begin
	# Plot functions
	
	function plot_input(X; kwargs...)
		n_neur = length(X)
		d_neur = Int(ceil(n_neur/6, digits = 0))
		return plot([1, 1] * vcat(X...)',
    		vcat([i .* ones(length(X[i])) for i = 1:n_neur]...)' .+ [-0.45, 0.45];
			linewidth = 2,
			xlim = (0, t_final), ylim = (0.5, length(X) + 0.5),
			yticks = d_neur:d_neur:n_neur,
			xlabel = "t [ms]", ylabel = "Neuron. No.",
			label = "",
			kwargs...)
	end
	
	function plot_voltage(V_b, V_a, θ, color; kwargs...)
		p_volt = plot([0, t_final], θ*[1, 1]; 
			color = :black, linestyle = :dash, 
			xlabel = "t [ms]", ylabel = "V [mV]",
			xlim = (0, t_final), ylim = (-Inf, max_V),
			label = "", kwargs...)
		plot!(t, V_b, 
			color = color,
			label = "before update")
		plot!(t, V_a,
			linestyle = :dash,
			color = color,
			label = "after update")
		return p_volt
	end
	
end;

# ╔═╡ bce5d730-3cc7-11eb-338a-15172a215b22
begin
	η = 10.0^log_η
end;

# ╔═╡ 4ff22490-3cad-11eb-3ed5-1b2c43d4221e
begin
	redraw_sample
	
	# Choose a random sample
	sample = rand(Samples)
	
end;

# ╔═╡ 20d32260-3cbf-11eb-2fae-0f0b8823ad80
begin
	
	# Get the voltage function
	PSPs_max, spk = Tempotrons.get_binary_training_potential(tmp, sample.x)[2:3]
	if spk
		PSPs = PSPs_max
	else
		PSPs = sort(Tempotrons.get_psps(tmp, sample.x), by = x -> x.time)
	end
	Vₜᵣ(x) = isempty(PSPs) ? 0.0 : sum(psp -> psp.ΔV(x), PSPs)
	V_train = Vₜᵣ.(t)
	V_true = tmp(sample.x, t = collect(t)).V
	
	k_max = argmax(V_train)
	t_max = t[k_max]
	V_max = V_train[k_max]
	
	w₀ = copy(tmp.w)
	train!(tmp, sample.x, sample.y, optimizer = Optimizers.SGD(η))
	Δw = tmp.w .- w₀
	
	PSPs_max, spk2 = Tempotrons.get_binary_training_potential(tmp, sample.x)[2:3]
	if spk2
		PSPs = PSPs_max
	else
		PSPs = sort(Tempotrons.get_psps(tmp, sample.x), by = x -> x.time)
	end
	Vₜᵣ(x) = isempty(PSPs) ? 0.0 : sum(psp -> psp.ΔV(x), PSPs)
	V_train_new = Vₜᵣ.(t)
	V_true_new = tmp(sample.x, t = collect(t)).V
	
end;

# ╔═╡ efbd08f0-3cc1-11eb-2174-af477223385e
begin
	if spk == sample.y
		error_text = "\$\$\\varepsilon = 0\$\$"
		delta_text = "\\Delta w_i = 0"
		grad_text = "\\Delta\\vec{w} = " * latex_form(Δw)
	else
		error_exp = spk ? "V_{\\mathrm{max}}-\\theta" : "\\theta-V_{\\mathrm{max}}"
		error_val = round(abs(V_max - tmp.θ), sigdigits = 3)
		error_text = "\$\$\\varepsilon = " * error_exp * " = $error_val\\mathrm{mV}\$\$"
		
		delta_sign = sample.y ? "+" : "-"
		delta_exp = delta_sign * "\\eta\\sum_{t_i^j<t_{\\mathrm{max}}} K\\left(t_{\\mathrm{max}} - t_i^j\\right)"
		grad_val = latex_form(Δw)
		delta_text = "\\Delta w_i = " * delta_exp
		
		grad_text = "\\Delta\\vec{w} = " * latex_form(round.(Δw, sigdigits = 3))
	end
	ΔW_text = "\$\$" * delta_text * "\\quad,\\qquad" * grad_text * "\$\$"
	
	md"""
	 $error_text
	
	 $ΔW_text
	"""
end

# ╔═╡ 68533900-3caf-11eb-26cd-dd65fe477ff6
begin
	# Plots
	gr(size = (675, 500))
	
	# Input plot
	p_inp = plot_input(sample.x, color = sample.y ? :red : :blue)
	
	# Training voltage plot
	p_V_train = plot_voltage(V_train, V_train_new, tmp.θ, sample.y ? :red : :blue, 
		title = "Training Voltage", 
		yticks = ([0, tmp.θ, V_max], 
			["\$V_{\\mathrm{rest}}\$", "\$\\theta\$", "\$V_{\\mathrm{max}}\$"]),
		xticks = ([0, t_final, t_max], 
			["\$0\$", "\$$(t_final)\$", "\$t_{\\mathrm{max}}\$"]),
		legend = :none)
	plot!(t_max*[1, 1], [ylims()[1], V_max], 
		color = :gray, linestyle = :dash,
		label = "")
	plot!([0, t_max], V_max*[1, 1], 
		color = :gray, linestyle = :dash,
		label = "")
	if spk ≠ sample.y
		if spk
			from = max_V
		else
			from = ylims()[1]
		end
		quiver!([t_max], [from], quiver = ([0], [V_max - from]), 
			color = :green)
	end
	
	# True voltage plot
	p_V_true = plot_voltage(V_true, V_true_new, tmp.θ, sample.y ? :red : :blue, 
		title = "True Voltage", 
		yticks = ([0, tmp.θ], 
			["\$V_{\\mathrm{rest}}\$", "\$\\theta\$"]),
		legend = :none)
	
	# Legend
	p_leg = plot_voltage(fill(NaN, size(t)), fill(NaN, size(t)), tmp.θ, sample.y ? :red : :blue, framestyle = :none, legend = :right, color = :white)
	
	l = @layout [a; c; d]
	plot(p_inp, p_V_train, p_V_true, p_leg, 
		layout = grid(4,1, heights = [0.27, 0.35, 0.35, 0.03]), 
		link = :x, 
		rightmargin = 5mm)
	
end

# ╔═╡ Cell order:
# ╟─a2fffdb0-3ca8-11eb-0a28-7bb88a2f2292
# ╟─16b92fa0-3cbe-11eb-129d-456f9aeebd93
# ╟─efbd08f0-3cc1-11eb-2174-af477223385e
# ╟─68533900-3caf-11eb-26cd-dd65fe477ff6
# ╟─2148f98e-3f7b-11eb-031f-9d22f9d77b47
# ╟─d9ca7af0-3ca8-11eb-1d71-b571a6089f88
# ╟─5e7bb330-3cbe-11eb-0972-35792137cf8d
# ╟─37e9d3f0-3caa-11eb-33b2-79045f34aa25
# ╟─115700f0-3caf-11eb-287d-15083df59991
# ╟─bce5d730-3cc7-11eb-338a-15172a215b22
# ╟─4ff22490-3cad-11eb-3ed5-1b2c43d4221e
# ╟─20d32260-3cbf-11eb-2fae-0f0b8823ad80
