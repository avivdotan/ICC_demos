### A Pluto.jl notebook ###
# v0.12.16

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

# ╔═╡ 502bfabe-0a90-11eb-125e-b9ef39ea15c8
begin
	
	#Imports
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add("Plots")
	using LinearAlgebra
	using PlutoUI
	using Plots
	
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
	# Binary Perceptron - Demo
	"""
end

# ╔═╡ 96778920-0687-11eb-1b84-53175daf996e
begin
	# Parameters
	
	n_runs = [0]
	do_train = [true]
	
	N = 2
	width = 2
	
	P_min = 1
	P_max = 100
	ΔP = 1
	
	η_min = -2
	η_max = 0.7
	Δη = 0.1
	
	k_max = 100
	
	md"""
	\# of samples ``P`` $(@bind P Slider(P_min:ΔP:P_max, default = 50, show_value = true))
	``\quad``
	Separable samples $(@bind separable Select(["yes", "almost", "no"], default = "yes"))
	\
	$(@bind redraw_samples Button("Redraw samples"))
	$(@bind redraw_weights Button("Redraw weights"))

	Learning rate ``\eta`` $(@bind η_log LogSlider(η_min:Δη:η_max, default = -0.523, show_value = true))
	``\quad``
	Normalize weights each step $(@bind norm_weights CheckBox(default = false))

	 $(@bind train Button("Update weights")) ``\quad`` ``k`` $(@bind k Slider(1:k_max, default = 1, show_value = true))
	
	```math
	\Delta \vec{w} = \eta \:\, \frac{1}{2}\left(1 - y^{\mu} y_0^{\mu}\right) \:\, \vec{x}^{\mu} \:\, y_0^{\mu}
	```
	"""
end

# ╔═╡ 64cf4600-0683-11eb-10c1-29673fd7bfd2
begin
	# Generate samples
	function draw_samples(p, sep)
		
		y₀ = [1.0 for μ₀ = 1:p]
		x = [width .* (2 .* rand(Float64, N) .- 1) for μ₀ = 1:p]

		while length(unique(y₀)) < 2
			if sep == "yes"
				w₀ = rand(Float64, N)
				y₀ .= dropdims(sign.(w₀'*hcat(x...)), dims = 1)
			elseif sep == "almost"
				w₀ = rand(Float64, N)
				h = w₀'*hcat(x...)
				h .+= 0.3*(2*rand(Float64, size(h)) .- 1)
				y₀ .= dropdims(sign.(h), dims = 1)
			elseif sep == "no"
				y₀ .= rand([-1.0, 1.0], p)
			end
		end
		
		return (x, y₀)
	end
	
	# Network parameters
	function draw_weights()
		w = 2*rand(Float64, N) .- 1
		w ./= norm(w)
		return w
	end
	
	# Train network
	function train!(w₀, w_old₀, η₀, p, k₀, x, y₀, norm_w)
		
		w_old₀ .= w₀
		Δw₀ = zeros(size(w₀))
		μ₀ = 1
		for m = 1:k₀
			μ₀ = rand(1:p)
			yᵐ = sign(w₀'*x[μ₀])
			Δw₀ .= η₀*0.5(y₀[μ₀] - yᵐ)*x[μ₀]
			w₀ .+= Δw₀
			if norm_w
				w₀ ./= norm(w₀)
			end
		end
		
		return Δw₀, μ₀
	end
	
end;

# ╔═╡ 921ed020-3618-11eb-3d3c-5dd20f53a21b
begin
	redraw_samples
	
	n_runs[1] = 0
	X, Y₀ = draw_samples(P, separable)
	
end;

# ╔═╡ 039ab970-0a82-11eb-36fe-d54c560f92ae
begin
	redraw_weights
	
	n_runs[1] = 0
	w = draw_weights()
	w_old = copy(w)
	
end;

# ╔═╡ a5be8f1e-361e-11eb-3fdd-7993fa47558b
begin
	η_log, norm_weights, k
	
	do_train[1] = false
end;

# ╔═╡ 6ad361e0-3625-11eb-3155-730d21df035e
begin
	η_log
	
	η = 10.0^η_log
end;

# ╔═╡ 5bb65d9e-0684-11eb-1ca2-59966bccc085
begin
	train
	
	if n_runs[1] > 0 && do_train[1]
		Δw, μ = train!(w, w_old, η, P, k, X, Y₀, norm_weights)
	end
	plot_Δw = n_runs[1] > 0 && do_train[1] && k == 1
	n_runs[1] += 1
	do_train[1] = true
		
	function set_text_pos(x)
		horz = x[1] > 0 ? :left : :right
		vert = x[2] > 0 ? :top : :bottom
		return (horz, vert)
	end
	
	# Plots
	gr(size = (525, 400))
	
	X₁ = hcat(X[filter(μ -> Y₀[μ] == 1.0, 1:P)]...)
	X₂ = hcat(X[filter(μ -> Y₀[μ] == -1.0, 1:P)]...)
	
	# Plot current sample
	if plot_Δw
		scatter([X[μ][1]], [X[μ][2]], markersize = 9, 
			color = :lightgreen, markerstrokecolor = :lightgreen, label = "",
			legend = :outertopright, legendfontsize = 12)
		quiver!([0], [0], quiver = ([0.9X[μ][1]], [0.9X[μ][2]]), linewidth = 2, 
			color = :lightgreen)
		annotate!(X[μ][1]/2, X[μ][2]/2, 
			text("\$\\vec{x}^{\\mu}\$", :lightgreen, 
				set_text_pos([X[μ][1], X[μ][2]])...))
	else
		# Empty plot
		scatter([], [], label = "", 
			legend = :outertopright, legendfontsize = 12)
	end
	
	# Plot samples
	scatter!(X₁[1, :], X₁[2, :], color = :red, label = "\$y_0=+1\$")
	scatter!(X₂[1, :], X₂[2, :], color = :blue, label = "\$y_0=-1\$")
	
	# Plot previous weights
	if n_runs[1] > 1
		quiver!([0], [0], quiver = ([w_old[1]], [w_old[2]]), 
			color = :gray, linewidth = 2)
		plot!(x -> -(w_old[1]/w_old[2])*x, color = :gray, linestyle = :dash,
			linewidth = 1.3, label = "")
	end
	
	# Plot current weights
	quiver!([0], [0], quiver = ([w[1]], [w[2]]), color = :magenta, linewidth = 3)
	annotate!(w[1]/2, w[2]/2, text("\$\\vec{w}\$", :magenta, set_text_pos(w)...))
	plot!(x -> -(w[1]/w[2])*x, color = :cyan, linewidth = 2, 
		label = "\$\\vec{w}\\cdot\\vec{x}=0\$")
	
	# Plot weights' change
	if plot_Δw && any(Δw .!= 0.0)
		quiver!([w_old[1]], [w_old[2]], quiver = ([Δw[1]], [Δw[2]]), 
			color = :orange, linewidth = 3)
		annotate!(w_old[1] + Δw[1]/2, w_old[2] + Δw[2]/2, 
			text("\$\\Delta\\vec{w}\$", :orange, set_text_pos([w[1], -w[2]])...))
	end
	
	xlims!(-width, width)
	ylims!(-width, width)
	xlabel!("\$x_1\$")
	ylabel!("\$x_2\$")
end

# ╔═╡ Cell order:
# ╟─502bfabe-0a90-11eb-125e-b9ef39ea15c8
# ╟─96778920-0687-11eb-1b84-53175daf996e
# ╟─5bb65d9e-0684-11eb-1ca2-59966bccc085
# ╟─64cf4600-0683-11eb-10c1-29673fd7bfd2
# ╟─921ed020-3618-11eb-3d3c-5dd20f53a21b
# ╟─039ab970-0a82-11eb-36fe-d54c560f92ae
# ╟─a5be8f1e-361e-11eb-3fdd-7993fa47558b
# ╟─6ad361e0-3625-11eb-3155-730d21df035e
