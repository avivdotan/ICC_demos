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

# ╔═╡ 9653d710-269d-11eb-0499-3d429fa4086e
begin
	# Imports
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add("Latexify")
	Pkg.add("ForwardDiff")
	Pkg.add("Plots")
	# Pkg.add("PlotlyJS")
	using PlutoUI
	using Latexify
	using ForwardDiff
	using LinearAlgebra
	using Plots
	
	# Latexify defaults
	set_default(cdot = false, fmt = FancyNumberFormatter(3))
	
	# Plots default colors
	plot_cols = palette(:default)
	
	# Playground limits
	box_lim = 10
	z_lim = 10
	
	md"""
	# Gradient Descent - Demo
	"""
end

# ╔═╡ 718d6500-26b5-11eb-1d2e-27280dc0dab3
md"""
``f\left(x, y\right) = `` $(@bind fₑₓᵗ TextField((75, 1), default = "2x^4 + x^3*y + 2y^4 - 4(x + 5)^3 + 3(y + 7)^3 + 10(x + 5)^2 - 7(y - 7)^2 - 1"))
"""

# ╔═╡ 24bac9a0-34b8-11eb-1beb-cfec8e68d627
begin
	fₑₓᵗ
	
	md"""
	Select a predefined function $(@bind f_type Select(
		["ex" => "User's", 
		 "miriam3rd" => "Paraboloid", 
		 "quad" => "Another paraboloid",
		 "gmin" => "Simple", 
		 "2min" => "Two minima", 
		 "complex" => "Complex"
		], default = "ex"))
	"""
end

# ╔═╡ b0df8b00-2cf2-11eb-2a76-1d4bcfbed48a
begin
	
	# Define a function
	if f_type == "ex"
		
		fₑₓ = Meta.parse(fₑₓᵗ)
		
	elseif f_type == "miriam3rd"

		# fₑₓ = :(x^2 + x*y + 9y^2)
		fₑₓ = :(2x^2+x*y+4y^2)
		
	elseif f_type == "quad"

		fₑₓ = :(x^2 - 4x*y + 5y^2)
	
	elseif f_type == "gmin"
		
		fₑₓ = :(exp(-((x^2 + y^2)/70))*
				sin(-(x/30)) - 
				0.01x)
		
	elseif f_type == "2min"
		
		fₑₓ = :(-exp(-(((x - 4)^2 + (y + 4)^2)/50)) - 
				(3/4)*exp(-(((x + 4)^2 + (y - 4)^2)/30)))
		
	elseif f_type == "complex"
		
		fₑₓ = :(exp(-((x^2 + y^2)/70))*
				cos((x + y)/3)*
				cos((x - y)/10)*
				sin((y^2 - x^2)/300))
		
	end
	
	# Convert to f(x⃗)
	eval(:(fˣʸ(x, y) = $fₑₓ))
	fᵗ(x) = fˣʸ(x[1], x[2])

	# renormalize for better display
	mx = -box_lim:0.1:box_lim
	my = -box_lim:0.1:box_lim
	f_vals = [fᵗ([x, y]) for x ∈ mx, y ∈ my]
	f_min, f_max = minimum(f_vals), maximum(f_vals)
	f(x) = 5(1/2 + (fᵗ(x) - f_min)/(f_max - f_min))
	
	# Define the gradient
	∇f(x) = ForwardDiff.gradient(f, x)
	
	md"""
	``f\left(x,y\right) = `` $(latexify(fₑₓ))
	
	Update rules:
	
	```math
	\begin{eqnarray}
	\\
	&&\Delta x_n &= -\eta \left. \frac{\partial f}{\partial x} \right\vert_{\left(x_n, y_n\right)}&&
	,\quad
	\Delta y_n &= -\eta \left. \frac{\partial f}{\partial y} \right\vert_{\left(x_n, y_n\right)}&&
	\\
	&&x_{n+1} &= x_n + \Delta x_n&&
	,\quad
	y_{n+1} &= y_n + \Delta y_n&&
	\\
	\end{eqnarray}
	```
	"""
end

# ╔═╡ b7f850a0-34d3-11eb-16f9-e3a421c96735
md"""
``\:`` $(@bind redraw_w₀ Button("Draw initial guess"))
``\qquad``
``\eta`` $(@bind η Slider(0.1:0.1:10, default = 1, show_value = true))
``\qquad``
\# of steps $(@bind T Slider(0:1000, default = 0, show_value = true))
"""

# ╔═╡ 5e198190-2cfd-11eb-2b6e-93daedb94e98
begin
	redraw_w₀
	
	dw = 0.1
	
	md"""
	``x_0`` $(@bind x₀ Slider((-box_lim + dw):dw:(box_lim - dw), default = rand((-box_lim + dw):dw:(box_lim - dw)), show_value = true))
	\
	``y_0`` $(@bind y₀ Slider((-box_lim + dw):dw:(box_lim - dw), default = rand((-box_lim + dw):dw:(box_lim - dw)), show_value = true))
	"""
end

# ╔═╡ a2c8c110-26af-11eb-0577-25ecbeb860e8
begin
	w₀ = [x₀, y₀]
	w_hist = [w₀]
	w = copy(w₀)
	for t = 1:T
		w .-= η*∇f(w)
		push!(w_hist, copy(w))
	end
	
	curr_∇f = ∇f(w)
	next_Δw = -η*curr_∇f
	
	x_names = ["x", "y"]
	x_expr(i) = "$(x_names[i])_{$T} &= $(latexify(w[i], env = :raw))&"
	pos_str = "&&\\text{Current position:}\\quad&" * x_expr(1) * ",\\quad " * x_expr(2)
	∇x_expr(i) = "\\left. \\frac{\\partial f}{\\partial $(x_names[i])} \\right\\vert_{\\left($(x_names[1])_{$T}, $(x_names[2])_{$T}\\right)} &= $(latexify(curr_∇f[i], env = :raw))&"
	grad_str = "&&\\text{Current gradient:}\\quad&" * ∇x_expr(1) * ",\\quad " * ∇x_expr(2)
	Δx_expr(i) = "\\Delta $(x_names[i])_{$T} &= $(latexify(next_Δw[i], env = :raw))&"
	delta_str = "&&\\text{Next update:}\\quad&" * Δx_expr(1) * ",\\quad " * Δx_expr(2)
	status_str = "\$ \\begin{eqnarray} \\\\ " * pos_str * "\\\\" * grad_str * "\\\\" * delta_str * " \\\\ \\end{eqnarray} \$"
	
	md"""
	 $status_str
	
	**Note:** ``z``-axiz is automatically  *shifted* and *rescaled* for visualization purposes.
	"""
end

# ╔═╡ 30ee8b80-269e-11eb-1d38-6928957188f3
begin
	#Plots
	gr(size = (675, 300))
	
	# Wireframe grid
	dx = box_lim/20
	x = -box_lim:dx:box_lim
	y = -box_lim:dx:box_lim
	
	# Split vectors history by dimension
	w₁ = [w[1] for w ∈ w_hist[1:end - 1]]
	w₂ = [w[2] for w ∈ w_hist[1:end - 1]]
	
	# Prepare gradients for display
	Δw₂ = [(-η*∇f(w)..., ) for w ∈ w_hist[1:end - 1]]
	Δw₃ = [(Δw..., 0) for Δw ∈ Δw₂]
	
	# 3D Plot
	p3d = wireframe(x, y, (x, y) -> f([x, y]), 
		camera = (30, 40), 
		linecolor = :black,
		xlim = (-box_lim, box_lim), ylim = (-box_lim, box_lim), zlim = (0, z_lim), 
		xlabel = "x", ylabel = "y", zlabel = "f(x,y)")
	quiver!(p3d, w₁, w₂, zeros(size(w₁)), quiver = Δw₃, 
		linewidth = 3, linecolor = plot_cols[2])
	plot!(p3d, w_hist[end][1]*[1, 1], w_hist[end][2]*[1, 1], [0, f(w_hist[end])], 
		linestyle = :dash, linewidth = 3, linecolor = plot_cols[3],
		label = "")
	scatter!(p3d, [w_hist[end][1]], [w_hist[end][2]], [f(w_hist[end])], 
		markersize = 5, markercolor = plot_cols[4], 
		label = "")
	
	# 2D Plot
	p2d = contour(x, y, (x, y) -> f([x, y]),
		color = :rainbow, levels = 40, aspect_ratio = :equal,
		xlim = (-box_lim, box_lim), ylim = (-box_lim, box_lim), 
		xlabel = "x", ylabel = "y")
	quiver!(p2d, w₁, w₂, quiver = Δw₂, 
		linecolor = plot_cols[2])
	plot!(p2d, [], [], label = "")
	scatter!(p2d, [w_hist[end][1]], [w_hist[end][2]], 
		markersize = 5, markercolor = plot_cols[4], 
		label = "")
	
	plot(p3d, p2d, layout = grid(1, 2, widths = [0.45 ,0.65]))
end

# ╔═╡ Cell order:
# ╟─9653d710-269d-11eb-0499-3d429fa4086e
# ╟─718d6500-26b5-11eb-1d2e-27280dc0dab3
# ╟─24bac9a0-34b8-11eb-1beb-cfec8e68d627
# ╟─b0df8b00-2cf2-11eb-2a76-1d4bcfbed48a
# ╟─b7f850a0-34d3-11eb-16f9-e3a421c96735
# ╟─5e198190-2cfd-11eb-2b6e-93daedb94e98
# ╟─a2c8c110-26af-11eb-0577-25ecbeb860e8
# ╟─30ee8b80-269e-11eb-1d38-6928957188f3
