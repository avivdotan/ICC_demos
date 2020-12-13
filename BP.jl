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

# ╔═╡ 29fbdb70-362a-11eb-0ca4-c5c0e0621c86
begin
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add("Latexify")
	Pkg.add("Flux")
	Pkg.add("Plots")
	using PlutoUI
	using Latexify
	using Statistics
	using LinearAlgebra
	using Flux
	using Plots
	
	# Latexify default
	set_default(fmt = FancyNumberFormatter(3))
	
	# Plot colors
	plot_cols = palette(:default)
	
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
	# Back Propagation - Function Approximation
	"""
end

# ╔═╡ fcbada60-37e3-11eb-2ead-b7debfc7f291
md"""
``f\left(x\right) = `` $(@bind f⁼ TextField(default="sin(x)"))
"""

# ╔═╡ 70aba900-37f7-11eb-18ab-816ef839e2b1
begin
	
	x_range_min = -3π
	x_range_max = 3π
	Δx_range = 0.01
	
	md"""
	``x_{\mathrm{min}} = `` $(@bind x_min Slider(x_range_min:Δx_range:x_range_max, default = round(2x_range_min/3 + x_range_max/3, digits = 2), show_value = true))
	"""
end

# ╔═╡ 0a58d080-37f5-11eb-3cd5-4552c8f238ee
md"""
``x_{\mathrm{max}} = `` $(@bind x_max Slider(x_min:Δx_range:x_range_max, default = round((x_range_max + x_min)/2, digits = 2), show_value = true))
"""

# ╔═╡ 37112102-37e5-11eb-0c3c-c3d7be63f9bd
begin
	fₑₓ = Meta.parse(f⁼)
	eval(:(f(x) = $fₑₓ))
	f
	x_lims = [x_min, x_max]
	
	md"""
	``f\left(x\right) = `` $(latexify(fₑₓ))
	
	Domain: ``\left[x_{\mathrm{min}},x_{\mathrm{max}}\right] = [`` $(latexify(x_lims[1])) ``,`` $(latexify(x_lims[end])) ``]``
	
	\# of samples: $(@bind P Slider(1:2000, default = 500, show_value = true))
	``\quad``
	Noise level: $(@bind log_noise LogSlider(-3:0.1:1, default = -1, show_value = true))
	"""
end

# ╔═╡ 766d5460-3667-11eb-1312-c9dc875c5e4b
md"""
| Hyperparameter | Value |
|:--- |:--- |
| \# of layers | $(@bind n_layers Slider(1:10, default = 2, show_value = true)) |
| \# of neurons per hidden layer | $(@bind n_hidden Slider(1:30, default = 5, show_value = true)) |
| Learning rate | $(@bind log_η LogSlider(-2:0.1:0, default = -1, show_value = true)) |
| \# of training epochs | $(@bind n_epochs Slider(10:10:500, default = 50, show_value = true)) |
| Activation function | ``g\left(x\right) = `` $(@bind g⁼ TextField(default = "1/(1 + exp(-x))")) |
"""

# ╔═╡ 967c5de0-3889-11eb-174c-774e35793b36
begin
	if n_layers == 1
		n_params = 2
	else
		n_params = 2n_hidden + (n_layers - 2)*n_hidden*(n_hidden +  1) + n_hidden + 1
	end
	
	md"""
	Total \# of learned parameters: $n_params
	"""
end

# ╔═╡ e5da3880-37ee-11eb-2ecf-8b63d6992e24
begin
	gₑₓ = Meta.parse(g⁼)
	eval(:(gfunc(x) = $gₑₓ))
	gfunc
	
	latexify(:(g(x) = $gₑₓ))
end

# ╔═╡ 0dda4fe0-3665-11eb-1106-69c1457aab9c
begin
	# Functions
	
	function generate_data(xlims, func; n_samples = 1, noise_std = 1)
		x = (xlims[2] - xlims[1])*rand(Float64, n_samples) .+ xlims[1]
		y = func.(x) .+ noise_std*randn(size(x))
		return x, y
	end
	
	function model_outputs(m, x)
		return [m(x[i, :])[1] for i = 1:size(x, 1)]
	end
	
	function model_loss(cost, x, y)
		return mean([cost(x[i, :], y[i]) for i = 1:length(y)])
	end
	
	function build_train_model(L, Nₕ, g, x, y, lr, Nₑ)
	
		# Define a network
		if L == 1
			layers = [Dense(1, 1)]
		else
			layers = [Dense(1, Nₕ, g), 
					  [Dense(Nₕ, Nₕ, g) for l = 1:(L - 2)]...,
					  Dense(Nₕ, 1)]
		end
		m = Chain(layers...)
		
		# Define model parameters
		ps = Flux.params(m)

		# Define loss function
		cost(x₀, y₀) = Flux.Losses.mse(m(x₀), y₀)

		# Define the data loading process
		data = Flux.Data.DataLoader(x', y, shuffle = true)
		
		# Train
		hist = [(params	= deepcopy(ps), 
			     loss 	= model_loss(cost, x, y))]
		t_time = @elapsed for e = 1:Nₑ
			Flux.train!(cost, ps, data, Descent(lr))
			push!(hist, 
				  (params 	= deepcopy(ps), 
			       loss 	= model_loss(cost, x, y)))
		end
		
		return m, t_time, hist
	end
	
end;

# ╔═╡ 115ab2e0-36fb-11eb-31ef-ef2d29e4c9b3
begin
	# Plot functions
	
	function plot_data(x, y; 
					   color = plot_cols[1])
		# Data
		p_data = scatter(x, y, 
			color = color, 
			markersize = 2,
			xlabel = "\$x\$", ylabel = "\$y\$",
			label = "Data",
			legend = :none)
		return p_data
	end
	
	function plot_predictions!(base_p, xᵣ, m; 
						  	  color = plot_cols[1])
		# Network predictions
		p_pred = plot!(base_p, xᵣ, (x) -> (m([x])[1]),
			color = color,
			linewidth = 3,
			label = "Network predictions", 
			legend = :topright)
	end
	
	function plot_history(hist, ep; 
						  color = plot_cols[1])
		# Training history plot
		
		cost_hist = [h.loss for h ∈ hist]
		
		# Losses
		p_history = plot(cost_hist, 
			color = color,
			ylim = (0, Inf),
			xlabel = "epochs", ylabel = "training loss",
			label = "")
		# Current epoch
		vline!([ep], 
			color = :darkgray, 
			label = "")
		scatter!([ep], [cost_hist[ep]], 
			color = color, 
			label = "")
		
		return p_history
	end
	
	function plot_architecture(L, Nₕ; 
							   node_color = plot_cols[1], 
							   weight_color = plot_cols[2])
		# Plot architecture
		aspectratio = 0.6
		
		x_width = Nₕ + 1
		x_pos = [x_width/2, 
				 repeat(1:Nₕ, L - 1)..., 
				 x_width/2]
		y_pos = [L + 1, 
				 hcat([l*ones(Nₕ) for l = L:-1:2]...)..., 
				 1]
		y_pos = y_pos .- 0.4
		y_pos *= aspectratio*(x_width)/(L + 1.2)
		node_size = 20
		hidden_node_size = 5
		sizes = [node_size, 
				 repeat([hidden_node_size], Nₕ*(L - 1))..., 
				 node_size]
		nodes_text = ["\$x\$", 
					  repeat([""], Nₕ*(L - 1))..., 
					  "\$y\$"]
		nodes_text = text.(nodes_text, node_color)

		# Plot synapses
		node_fac = 0.034*(Nₕ + 1)
		hidden_node_fac = 0.01*(Nₕ + 1)
		if L == 1
			xs = [x_pos[1]]
			ys = [y_pos[1]]
			us = [x_pos[end]] .- xs
			vs = [y_pos[end]] .- ys
			source_facs = [node_fac]
			target_facs = [node_fac]
		else
			xs = vcat([x_pos[i]*ones(Nₕ) 
					   for i = 1:(length(x_pos) - Nₕ - 1)]..., 
					  x_pos[end .- (Nₕ:-1:1)]...)
			ys = vcat([y_pos[i]*ones(Nₕ) 
					   for i = 1:(length(y_pos) - Nₕ - 1)]...,
					  y_pos[end .- (Nₕ:-1:1)]...)
			us = [x_pos[1 .+ (1:Nₕ)]..., 
				  hcat([repeat(x_pos[(1 + l*Nₕ) .+ (1:Nₕ)], Nₕ)
				   for l = 1:(L - 2)]...)...,
				  x_pos[end]*ones(Nₕ)...] .- xs
			vs = [y_pos[1 .+ (1:Nₕ)]..., 
				  hcat([repeat(y_pos[(1 + l*Nₕ) .+ (1:Nₕ)], Nₕ)
				   for l = 1:(L - 2)]...)...,
				  y_pos[end]*ones(Nₕ)...] .- ys
			source_facs = [node_fac*ones(Nₕ)..., 
						   hidden_node_fac*ones(length(xs) - Nₕ)...]
			target_facs = [hidden_node_fac*ones(length(xs) - Nₕ)...,
						   node_fac*ones(Nₕ)...]
		end
		vec_norms = norm.(collect(zip(us, vs)))
		# source_facs[source_facs .== node_fac] ./= vec_norms[source_facs .== node_fac]
		# target_facs[target_facs .== node_fac] ./= vec_norms[target_facs .== node_fac]
		source_facs ./= vec_norms
		target_facs ./= vec_norms
		# xs += source_facs.*us
		# ys += source_facs.*vs
		# us = (1 .- source_facs .- target_facs).*us
		# vs = (1 .- source_facs .- target_facs).*vs
		us = (1 .- target_facs).*us
		vs = (1 .- target_facs).*vs

		# Do the actual plotting
		p_arch = quiver(xs, ys, quiver = (us, vs), 
			color = weight_color, linewidth = 0.2,
			xlim = (0, x_width), ylim = (0, aspectratio*x_width),
			framestyle = :none, aspect_ratio = :equal)
		scatter!(x_pos, y_pos, 
			markersize = sizes, 
			color = :white, markerstrokewidth = 0,
			label = "")
		scatter!(x_pos, y_pos, 
			markersize = sizes,
			series_annotations = nodes_text, 
			color = node_color, alpha = 0.3, 
			label = "")
		
		return p_arch
	end
	
end;

# ╔═╡ b2275500-36fd-11eb-0822-a7bdcf1e1e59
begin
	# Generate Data
	noise = 10.0^log_noise
	X, Y₀ = generate_data(x_lims, f; n_samples = P, noise_std = noise)
	
	# Plot data
	gr(size = (675, 250))
	
	plot_data(X, Y₀, color = plot_cols[1])
end

# ╔═╡ 93cad06e-362e-11eb-1987-e117ac0e3693
begin
	gₑₓ
	η = 10.0^log_η
	
	model, train_time, history = 
		build_train_model(n_layers, n_hidden, gfunc, 
						  X, Y₀,  
						  η, n_epochs)
	trained = true
	
	md"""
	Training time: $(round(train_time, digits = 2)) seconds. 
	"""
end

# ╔═╡ 4e9b9cf0-3665-11eb-03fb-554d5309c2bc
begin
	trained
	
	md"""
	Epoch $(@bind epoch Slider(0:n_epochs, default = n_epochs, show_value = true))
	"""
end

# ╔═╡ 81ec2090-3663-11eb-2d6f-49cf07a71df2
begin
	# Performance
	trained
	
	md"""
	Training loss: $(round(Float64(history[epoch + 1].loss), sigdigits = 3))
	"""
end

# ╔═╡ 314c5870-36f0-11eb-2638-737a6f8591f7
begin
	# Plot architecture
	gr(size = (675, 400))
	plot_architecture(n_layers, n_hidden, 
		node_color = plot_cols[3], weight_color = plot_cols[4])
	title!("Network architecture")
end

# ╔═╡ 1da0c6f0-3658-11eb-2727-596397e80a55
begin
	trained
	
	gr(size = (675, 500))
	
	# Plot spirals
	Flux.loadparams!(model, history[epoch +  1].params)
	
	p_d = plot_data(X, Y₀, color = plot_cols[1])
	plot_predictions!(p_d, x_min:0.01:x_max, model, color = plot_cols[2])
	
	# Plot training history
	p_h = plot_history(history, epoch + 1, color = plot_cols[6])
	
	plot(p_d, p_h, layout = (2, 1))
	
end

# ╔═╡ Cell order:
# ╟─29fbdb70-362a-11eb-0ca4-c5c0e0621c86
# ╟─fcbada60-37e3-11eb-2ead-b7debfc7f291
# ╟─70aba900-37f7-11eb-18ab-816ef839e2b1
# ╟─0a58d080-37f5-11eb-3cd5-4552c8f238ee
# ╟─37112102-37e5-11eb-0c3c-c3d7be63f9bd
# ╟─b2275500-36fd-11eb-0822-a7bdcf1e1e59
# ╟─314c5870-36f0-11eb-2638-737a6f8591f7
# ╟─967c5de0-3889-11eb-174c-774e35793b36
# ╟─e5da3880-37ee-11eb-2ecf-8b63d6992e24
# ╟─766d5460-3667-11eb-1312-c9dc875c5e4b
# ╟─93cad06e-362e-11eb-1987-e117ac0e3693
# ╟─4e9b9cf0-3665-11eb-03fb-554d5309c2bc
# ╟─81ec2090-3663-11eb-2d6f-49cf07a71df2
# ╟─1da0c6f0-3658-11eb-2727-596397e80a55
# ╟─0dda4fe0-3665-11eb-1106-69c1457aab9c
# ╟─115ab2e0-36fb-11eb-31ef-ef2d29e4c9b3
