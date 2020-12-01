### A Pluto.jl notebook ###
# v0.12.15

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

# ╔═╡ 08778ee0-2a77-11eb-296c-453f1932bf58
md"""
Now we can have a look at the gradients of ``f(x,y)`` at many different points simultaneously: 

**Note:** The displayed *gradient* vectors are *rescaled* for better visualization. 
"""

# ╔═╡ 6fa8a130-29af-11eb-0f35-776973ddd076
begin
	# Imports
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add("Latexify")
	Pkg.add("ForwardDiff")
	Pkg.add("Plots")
	Pkg.add("PlotlyJS")
	using PlutoUI
	using Latexify
	using ForwardDiff
	using LinearAlgebra
	using Plots
	using Plots.PlotMeasures
	
	# Latexify defaults
	set_default(cdot = false, fmt = FancyNumberFormatter(3))
	
	# Plots default colors
	plot_cols = palette(:default)
	cmap = :rainbow
	
	# Playground limits
	box_lim = 10
	z_lim = 10
	
	md"""
	# Partial Derivatives, Directional Derivatives and Gradients
	
	## Multivariable functions
	
	**Multivariable functions**  are mappings of inputs from a high-dimensional domain to the real line, *i.e.* ``f:\mathbb{R}^n\to\mathbb{R}``. For simplicity, we will focus here on functions from a 2-dimensional domain, *i.e.*  ``f:\mathbb{R}^2\to\mathbb{R}``, though the definitions and results are quite general. 
	
	For example:
	\
	``f\left(x, y\right) = `` $(@bind fₑₓᵗ TextField(default = "x^2 + 3x*y + y - 1"))
	
	Select a predefined function $(@bind f_type Select(["ex" => "User's", "quad" => "Paraboloid", "gmin" => "Simple", "complex" => "Complex"], default = "ex"))
	"""
end

# ╔═╡ 75a0a510-29af-11eb-2005-dfb6db948065
begin
	
	# Define a function
	if f_type == "ex"
		
		fₑₓ = Meta.parse(fₑₓᵗ)
		
	elseif f_type == "quad"

		# C = [ 1  -2 
		# 	 -2  5]
		# f(x) = (0.02x'*C*x + 10)/4
		fₑₓ = :(x^2 - 4x*y + 5y^2)
		
	elseif f_type == "gmin"
		
		fₑₓ = :(exp(-((x^2 + y^2)/70))*
				sin(-(x/30)) - 
				0.01x)
		
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
	``f(x,y) = `` $(latexify(fₑₓ))
	
	**Note:** ``z``-axiz is automatically  *shifted* and *rescaled* for visualization purposes.
	"""
end

# ╔═╡ 3fab8fb0-29f0-11eb-058e-b1717bab59bd
begin
	# Plots Functions
	
	function xy_grid()
		# xy grid
		dx = box_lim/20
		gx = -box_lim:dx:box_lim
		gy = -box_lim:dx:box_lim
		return (gx, gy)
	end
	
	function plot_3d_base(gx, gy)

		# 3D plot
		p3d = surface(gx, gy, (x, y) -> f([x, y]), 
			# camera = (30, 40), 
			color = cmap,
			seriesalpha = 0.8,
			xlim = (-box_lim, box_lim), ylim = (-box_lim, box_lim), zlim = (0, z_lim), 
			xlabel = "x", ylabel = "y", zlabel = "f(x,y)", 
			colorbar = :none)
		
		return p3d
		
	end
	
	function plot_2d_base(gx, gy)
		
		# 2D plot
		p2d = contour(gx, gy, (x, y) -> f([x, y]),
			color = cmap, levels = 40,
			aspect_ratio = :equal,
			xlim = (-box_lim, box_lim), ylim = (-box_lim, box_lim), 
			xlabel = "x", ylabel = "y", 
			colorbar = :none)
		
		return p2d
	end
	
end;

# ╔═╡ ed15ed70-29ec-11eb-380f-139b9e82f0f6
begin
	
	plotly(size = (675, 750))
	
	gx, gy = xy_grid()
	
	# 3D plot
	p3d = plot_3d_base(gx, gy)
	title!("f(x,y)")
	
	# 2D plot
	p2d = plot_2d_base(gx, gy)
	title!("Contour plot")
	
	plot(p3d, p2d, layout = grid(2, 1, heights = (0.7, 0.3)))
	
end

# ╔═╡ 868e4b20-29af-11eb-289e-ed49c852493d
begin
	# Functions definitions
	
	function s_grid(x₀, y₀, u)
		# 1D axis
		ds = 0.01
		s = √(2)*(-2box_lim:ds:2box_lim)
		s = s[(-box_lim .< (x₀ .+ u[1]*s) .< box_lim) .& 
			  (-box_lim .< (y₀ .+ u[2]*s) .< box_lim)]
		return s
	end
	
	function t_grid(x, y)
		# 1D axis
		dt = 0.01
		t = -5box_lim:dt:5box_lim
		t_in(τ) = (-box_lim < x(τ) < box_lim) & 
			  	  (-box_lim < y(τ) < box_lim)
		t_ind = t_in.(t)
		
		t_lims_ind = []
		if t_in(t[1])
			push!(t_lims_ind, 1)
		end
		for i = 1:(length(t) - 1)
			if t_ind[i] ≠ t_ind[i + 1]
				push!(t_lims_ind, t_in(t[i]) ? i : (i + 1))
			end
		end
		if t_in(t[end])
			push!(t_lims_ind, length(t))
		end
		
		t_lims = t[t_lims_ind]
		t = t[t_ind]
		
		let i = 1
			while i < length(t)
				if t[i] ∈ t_lims && t[i + 1] ∈ t_lims
					insert!(t, i + 1, NaN)
				end
				i += 1
			end
		end
		
		return t, t_lims
	end
	
	function plot_3d_plane!(p3d, x₀, y₀, u, s, x₁, x₂, y₁, y₂, col)

		# 3D plot
		plot!(p3d, x₀ .+ s*u[1], y₀ .+ s*u[2], 
			f.([[x₀ .+ sᵢ*u[1], y₀ .+ sᵢ*u[2]] for sᵢ ∈ s]), 
			linetype = :path3d, linewidth = 5,
			color = col, 
			label = "")
		lin_app = f([x₀, y₀]) .+ ∇f([x₀, y₀])'*u*s
		ind = (0 .<= lin_app .<= z_lim)
		plot!(p3d, x₀ .+ s[ind]*u[1], y₀ .+ s[ind]*u[2], 
			lin_app[ind],
			seriescolor = col, linewidth = 5, linestyle = :dash,
			label = "")
		plot!([x₁, x₁, x₂, x₂, x₁],
			  [y₁, y₁, y₂, y₂, y₁],
			  [0, z_lim, z_lim, 0, 0], 
			seriestype = :path3d,
			seriescolor = col, linewidth = 5, fill = 0,
			label = "")
		
	end
	
	function plot_3d_trajectory!(p3d, x, y, t, t_lims, col)
	
		# 3D plot
		plot!(p3d, x.(t), y.(t), 
			f.([[x(tᵢ), y(tᵢ)] for tᵢ ∈ t]), 
			linetype = :path3d, linewidth = 5,
			color = col, 
			label = "")
		plot!(p3d, x.(t), y.(t), 
			zeros(size(t)), 
			linetype = :path3d, linewidth = 5,
			color = col, 
			label = "")
		plot!(p3d, x.(t), y.(t), 
			z_lim .* ones(size(t)), 
			linetype = :path3d, linewidth = 5,
			color = col, 
			label = "")
		for tl ∈ t_lims
			plot!(p3d, x(tl) .* [1, 1], y(tl) .* [1, 1], 
				[0, z_lim], 
				linetype = :path3d, linewidth = 5,
				color = col, 
				label = "")
		end
		
	end
	
	function plot_3d_point!(p3d, x₀, y₀, col)

		# 3D plot
		scatter!(p3d, [x₀], [y₀], [f([x₀, y₀])], 
			markersize = 3, markercolor = col, 
			label = "")
		
	end
	
	function plot_2d_plane!(p2d, x₀, y₀, u, x₁, x₂, y₁, y₂, col)
	
		# 2D plot
		plot!(p2d, [x₁, x₂], [y₁, y₂], 
			linecolor = col, linewidth = 3, linestyle = :dash,
			label = "")
		quiver!(p2d, [x₀], [y₀], quiver = [(3u..., )],
			linecolor = col, linewidth = 5)
		
	end
	
	function plot_2d_trajectory!(p2d, x, y, t, col)
	
		# 2D plot
		plot!(p2d, x.(t), y.(t), 
			linecolor = col, linewidth = 3,
			label = "")
		
	end
	
	function plot_2d_point!(p2d, x₀, y₀, col)
		
		# 2D plot
		scatter!(p2d, [x₀], [y₀], 
			markersize = 5, markercolor = col, 
			label = "")
		
	end
	
	function plot_1d(x₀, y₀, u, s, col, x₀_col; s₀ = 0)
		
		# 1D plot
		p1d = plot(s, [f([x₀ .+ (sᵢ - s₀)*u[1], y₀ .+ (sᵢ - s₀)*u[2]]) for sᵢ ∈ s],
			seriescolor = col, linewidth = 3,
			ylim = (0, z_lim),
			xlabel = "s", ylabel = "f(x₀ + su)",
			label = "")
		plot!(p1d, s, f([x₀, y₀]) .+ ∇f([x₀, y₀])'*u*(s .- s₀),
			seriescolor = col, linewidth = 2, linestyle = :dash,
			label = "")
		scatter!(p1d, [s₀], [f([x₀, y₀])], 
			markersize = 5, markercolor = x₀_col, 
			label = "")
		
		return p1d
	end
	
	function plot_partials(x₀, y₀; col_x = plot_cols[1], 
						   col_y = plot_cols[2], 
						   x₀_col = plot_cols[4])

		# 1D axis
		sˣ = s_grid(x₀, y₀, [1, 0])
		sʸ = s_grid(x₀, y₀, [0, 1])
		
		# 2D axes
		(gx, gy) = xy_grid()

		# 3D Plot
		p3d = plot_3d_base(gx, gy)
		plot_3d_plane!(p3d, x₀, y₀, [1, 0], sˣ, -box_lim, box_lim, y₀, y₀, col_x)
		plot_3d_plane!(p3d, x₀, y₀, [0, 1], sʸ, x₀, x₀, -box_lim, box_lim, col_y)
		plot_3d_point!(p3d, x₀, y₀, x₀_col)

		# 2D Plot
		p2d = plot_2d_base(gx, gy)
		plot_2d_plane!(p2d, x₀, y₀, [1, 0], -box_lim, box_lim, y₀, y₀, col_x)
		plot_2d_plane!(p2d, x₀, y₀, [0, 1], x₀, x₀, -box_lim, box_lim, col_y)
		plot_2d_point!(p2d, x₀, y₀, x₀_col)

		# 1D Plot
		p1d_x = plot_1d(x₀, y₀, [1, 0], x₀ .+ sˣ, col_x, x₀_col, s₀ = x₀)
		xlabel!("x")
		ylabel!("f(x,y₀)")
		title!("∂f/∂x = $(round(∇f([x₀, y₀])[1], digits = 4))")
		p1d_y = plot_1d(x₀, y₀, [0, 1], y₀ .+ sʸ, col_y, x₀_col, s₀ = y₀)
		xlabel!("y")
		ylabel!("f(x₀,y)")
		title!("∂f/∂y = $(round(∇f([x₀, y₀])[2], digits = 4))")

		l = @layout([a{0.7h}; b{0.4w} c d])
		return plot(p3d, p2d, p1d_x, p1d_y, layout = l, 
			leftmargin = 3mm, rightmargin = 3mm)
		
	end
	
	function plot_directional(x₀, y₀, α; col = plot_cols[1], x₀_col = plot_cols[4])

		u = [cos(α), sin(α)]

		# 1D axis
		s = s_grid(x₀, y₀, u)
		x₁, x₂ = x₀ .+ s[1]*u[1], x₀ .+ s[end]*u[1]
		y₁, y₂ = y₀ .+ s[1]*u[2], y₀ .+ s[end]*u[2]
		
		# 2D axes
		(gx, gy) = xy_grid()

		# 3D Plot
		p3d = plot_3d_base(gx, gy)
		plot_3d_plane!(p3d, x₀, y₀, u, s, x₁, x₂, y₁, y₂, col)
		plot_3d_point!(p3d, x₀, y₀, x₀_col)

		# 2D Plot
		p2d = plot_2d_base(gx, gy)
		plot_2d_plane!(p2d, x₀, y₀, u, x₁, x₂, y₁, y₂, col)
		plot_2d_point!(p2d, x₀, y₀, x₀_col)
		annotate!(x₀ + 3u[1], y₀ + 3u[2], 
			text("u", :top, :left, col))

		# 1D Plot
		p1d = plot_1d(x₀, y₀, u, s, col, x₀_col)
		title!("Dᵤf = $(round(∇f([x₀, y₀])'*u, digits = 4))")

		l = @layout([a{0.7h}; b{0.4w} c])
		return plot(p3d, p2d, p1d, layout = l)
		
	end
	
	
	function plot_gradient(x₀, y₀, α; col_u = plot_cols[1], 
						   col_∇ = plot_cols[2], x₀_col = plot_cols[4])

		u = [cos(α), sin(α)]
		∇ᶠ = ∇f([x₀, y₀])
		∇ = ∇ᶠ./norm(∇ᶠ)

		# 1D axis for u
		su = s_grid(x₀, y₀, u)
		x₁u, x₂u = x₀ .+ su[1]*u[1], x₀ .+ su[end]*u[1]
		y₁u, y₂u = y₀ .+ su[1]*u[2], y₀ .+ su[end]*u[2]
		
		# 1D axis for ∇f
		s∇ = s_grid(x₀, y₀, ∇)
		x₁∇, x₂∇ = x₀ .+ s∇[1]*∇[1], x₀ .+ s∇[end]*∇[1]
		y₁∇, y₂∇ = y₀ .+ s∇[1]*∇[2], y₀ .+ s∇[end]*∇[2]
		
		# 2D axes
		(gx, gy) = xy_grid()

		# 3D Plot
		p3d = plot_3d_base(gx, gy)
		plot_3d_plane!(p3d, x₀, y₀, u, su, x₁u, x₂u, y₁u, y₂u, col_u)
		plot_3d_plane!(p3d, x₀, y₀, ∇, s∇, x₁∇, x₂∇, y₁∇, y₂∇, col_∇)
		plot_3d_point!(p3d, x₀, y₀, x₀_col)

		# 2D Plot
		p2d = plot_2d_base(gx, gy)
		plot_2d_plane!(p2d, x₀, y₀, u, x₁u, x₂u, y₁u, y₂u, col_u)
		plot_2d_plane!(p2d, x₀, y₀, ∇, x₁∇, x₂∇, y₁∇, y₂∇, col_∇)
		plot_2d_point!(p2d, x₀, y₀, x₀_col)
		annotate!(x₀ + 3u[1], y₀ + 3u[2], 
			text("u", :top, :left, col_u))
		annotate!(x₀ + 3∇[1], y₀ + 3∇[2], 
			text("∇f", :bottom, :right, col_∇))

		# 1D Plot for u
		p1du = plot_1d(x₀, y₀, u, su, col_u, x₀_col)
		title!("Dᵤf = $(round(∇f([x₀, y₀])'*u, digits = 4))")
		p1d∇ = plot_1d(x₀, y₀, ∇, s∇, col_∇, x₀_col)
		ylabel!("f(x₀ + s∇f)")
		title!("‖∇f‖ = $(round(norm(∇f([x₀, y₀])), digits = 4))")

		l = @layout([a{0.7h}; b{0.4w} c d])
		return plot(p3d, p2d, p1du, p1d∇, layout = l, 
			leftmargin = 3mm, rightmargin = 3mm)
		
	end
	
	function plot_chain_rule(x, y, dx, dy, t₀; 
							 col = plot_cols[1], x₀_col = plot_cols[4])

		# 1D axis
		t, t_lims = t_grid(x, y)
		t_in(τ) = (-box_lim < x(τ) < box_lim) & 
			  	  (-box_lim < y(τ) < box_lim)
		
		# 2D axes
		(gx, gy) = xy_grid()

		# 3D Plot
		p3d = plot_3d_base(gx, gy)
		plot_3d_trajectory!(p3d, x, y, t, t_lims, col)
		if t_in(t₀)
			plot_3d_point!(p3d, x(t₀), y(t₀), x₀_col)
		end

		# 2D Plot
		p2d = plot_2d_base(gx, gy)
		plot_2d_trajectory!(p2d, x, y, t, col)
		if t_in(t₀)
			plot_2d_point!(p2d, x(t₀), y(t₀), x₀_col)
		end

		# 1D Plot
		p1d = plot(t, [f([x(tᵢ), y(tᵢ)]) for tᵢ ∈ t],
			seriescolor = col, linewidth = 3,
			ylim = (0, z_lim),
			xlabel = "t", ylabel = "f",
			label = "")
		if t_in(t₀)
			scatter!(p1d, [t₀], [f([x(t₀), y(t₀)])], 
				markersize = 5, markercolor = x₀_col, 
				label = "")
		end
		title!("f(x(t), y(t))")
		
		# 1D Plot - derivative
		df(τ) = ([dx(τ) dy(τ)] * ∇f([x(τ), y(τ)]))[1]
		p1d_d = plot(t, df.(t),
			seriescolor = col, linewidth = 3, linestyle = :dash,
			xlabel = "t", ylabel = "df/dt",
			label = "")
		if t_in(t₀)
			scatter!(p1d_d, [t₀], [df(t₀)], 
				markersize = 5, markercolor = x₀_col, 
				label = "")
		end
		title!("d/dt f(x(t), y(t))")

		l = @layout([a{0.7h}; b{0.4w} c d])
		return plot(p3d, p2d, p1d, p1d_d, layout = l, 
			leftmargin = 3mm, rightmargin = 3mm)
		
	end
	
end;

# ╔═╡ 6990bc60-2a77-11eb-2415-b38966242435
begin
	plotlyjs(size = (675, 675))
	
	dg = 0.05
	gx₁ = -box_lim:dg:box_lim
	gy₁ = -box_lim:dg:box_lim
	plot_2d_base(gx₁, gy₁)
	
	gs = 1/10
	xg = (-box_lim + gs/2):(gs*box_lim):(box_lim - gs/2)
	yg = (-box_lim + gs/2):(gs*box_lim):(box_lim - gs/2)
	∇ᶠ = [∇f([x, y]) for x ∈ xg, y ∈ yg]
	∇ᶠ ./= maximum(norm.(∇ᶠ))
	∇ᶠ = [(gf..., ) for gf ∈ ∇ᶠ]
	quiver!(hcat([xg for y ∈ yg]...), vcat([yg' for x ∈ xg]...), quiver = ∇ᶠ, 
		linecolor = plot_cols[7], linewidth = 1.5)
end

# ╔═╡ 859be820-29ec-11eb-2053-97a2830c8a8d
md"""
## Partial derivatives

The **partial derivative** of the multivariable function ``f\left(x,y\right)`` with respect to ``x`` at the point ``\left(x_0, y_0\right)`` is defined by keeping ``y`` constant at ``y_0``, thus effectively turning ``f\left(x, y\right)`` into a *single variable* function ``f\left(x, y_0\right)``: 

```math
\left. \frac{\partial f}{\partial x}\right\vert_{\left(x_0, y_0\right)} = \left. \frac{\mathrm{d}}{\mathrm{d}x} f\left(x, y_0\right) \right|_{x = x_0}= \lim_{h\to 0} {\frac{f\left(x_0 + h, y_0\right) - f\left(x_0, y_0\right)}{h}}.
```

**Note:** The notation ``\partial`` is pronounced *"dee"* (same as a regular ``\mathrm{d}``), but in mathematics it is used for different purposes.  We will use ``\mathrm{d}`` for regular derivatives and ``\partial`` for partial derivatives. In MS Word or ``\LaTeX`` use `\partial` to produce it. 

Similarly, the partial derivative of ``f\left(x,y\right)`` with respect to ``y`` at the point ``\left(x_0, y_0\right)`` is defined by keeping ``x`` constant at ``x_0``: 

```math
\left. \frac{\partial f}{\partial y}\right\vert_{\left(x_0, y_0\right)} = \left. \frac{\mathrm{d}}{\mathrm{d}y} f\left(x_0, y\right) \right|_{y = y_0}= \lim_{h\to 0} {\frac{f\left(x_0, y_0 + h\right) - f\left(x_0, y_0\right)}{h}}.
```

Let's have a look at these derivatives! 
Choose a point: 
\
``x_0`` $(@bind x₀ Slider(-box_lim:0.1:box_lim, default = round(rand((-box_lim/2):0.1:(box_lim/2)), digits = 1), show_value = true))
\
``y_0`` $(@bind y₀ Slider(-box_lim:0.1:box_lim, default = round(rand((-box_lim/2):0.1:(box_lim/2)), digits = 1), show_value = true))
"""

# ╔═╡ 33e9bc80-29f3-11eb-11a7-b1cd60c8ae03
begin
	# Plot directional derivative
	plotlyjs(size = (675, 675))
	plot_partials(x₀, y₀, col_x = plot_cols[1], col_y = plot_cols[2])
end

# ╔═╡ cb617850-29ee-11eb-2e28-a99ff45f2ed1
md"""
## Directional derivative

The **directional derivative** of ``f\left(x,y\right)`` with respect to a *unit vector* ``\vec{u}=\left(u_1, u_2\right)`` at a point ``\left(x_0, y_0\right)`` is defined by keeping the function constant in any direction orthogonal to ``\vec{u}``, thus essentially turning ``f\left(x, y\right)`` into a *single variable* function (similarly to *partial derivatives*): 

```math
\left. D_{\vec{u}}f \right\vert_{\left(x_0, y_0\right)} = \left. \frac{\mathrm{d} f}{\mathrm{d} s}\right\vert_{\vec{u}, \left(x_0, y_0\right)} = \lim_{s\to 0} {\frac{f\left(x_0 + s\cdot u_1, y_0 + s\cdot u_2\right) - f\left(x_0, y_0\right)}{s}}\qquad \text{where } \left\Vert \vec{u} \right\Vert = 1.
```
	
**Note:** the *directional derivative* ``\left. D_{\left(1, 0\right)}f  \right\vert_{\left(x_0, y_0\right)}`` actually coincides with the definition of the *partial derivative* ``\left. \frac{\partial f}{\partial x} \right\vert_{\left(x_0, y_0\right)}``. Similarly, ``\left. D_{\left(0, 1\right)}f  \right\vert_{\left(x_0, y_0\right)} = \left. \frac{\partial f}{\partial y} \right\vert_{\left(x_0, y_0\right)}``.

Let's have a look at such a derivative! 
Change the direction: 
\
``θ`` $(@bind θₓ Slider(-1:(1/16):1, default = rand(-1:(1/16):1), show_value = true))π
"""

# ╔═╡ a1504642-29ee-11eb-25e9-97312b012e13
begin
	# Set θ units
	θ = θₓ*π
	
	md"""
	``\vec{u} = (`` $(latexify(cos(θ))) ``, `` $(latexify(sin(θ))) ``)``
	"""
end

# ╔═╡ b35b9120-29c9-11eb-1e50-8f2c7d0a635d
begin
	# Plot directional derivative
	plotlyjs(size = (675, 675))
	plot_directional(x₀, y₀, θ, col = plot_cols[3])
end

# ╔═╡ ccd642b0-2d8b-11eb-150a-f1ce567ecee1
md"""
## Chain rule

The chain rule for *univariate* functions ``f:\mathbb{R}\to\mathbb{R}`` and ``x:\mathbb{R}\to\mathbb{R}`` states that the derivative of ``f\left(x\left(t\right)\right)`` is:

```math
f'\left(x\left(t\right)\right)\dot{x}\left(t\right) \qquad \text{or} \qquad \frac{\mathrm{d} f}{\mathrm{d} t} = \frac{\mathrm{d} f}{\mathrm{d} x} \frac{\mathrm{d} x}{\mathrm{d} t}
```

The chain rule for a *multivariable* function ``f:\mathbb{R}^2\to\mathbb{R}`` and *univariate* functions ``x:\mathbb{R}\to\mathbb{R}`` and ``y:\mathbb{R}\to\mathbb{R}`` states that the derivative of ``f\left(x\left(t\right), y\left(t\right)\right)`` is:

```math
\frac{\mathrm{d} f}{\mathrm{d} t} = \frac{\partial f}{\partial x} \frac{\mathrm{d} x}{\mathrm{d} t} + \frac{\partial f}{\partial y} \frac{\mathrm{d} y}{\mathrm{d} t}
```

The functions ``x(t)`` and ``y(t)`` define a 2-dimentional *trajectory*. The derivative ``\frac{\mathrm{d} f}{\mathrm{d} t}`` is the derivative along this trajectory. 

Let's have a look at ``f`` and its derivative along the trajectory! Choose a trajectory:
\
``x(t) = `` $(@bind xₑₓᵗ TextField(default = "5t"))
``\qquad``
``y(t) = `` $(@bind yₑₓᵗ TextField(default = "t^3 - 3t"))
"""

# ╔═╡ e7612a0e-2d8f-11eb-3e0b-8da86e8dc20d
begin
	
	xₑₓ = Meta.parse(xₑₓᵗ)
	yₑₓ = Meta.parse(yₑₓᵗ)
	
	# Convert to x(t), y(t)
	eval(:(x(t) = $xₑₓ))
	eval(:(y(t) = $yₑₓ))
	
	# Define the derivatives
	dx(t) = ForwardDiff.derivative(x, t)
	dy(t) = ForwardDiff.derivative(y, t)
	
	md"""
	``x(t) = `` $(latexify(xₑₓ)) ``\qquad\qquad\qquad\qquad\qquad\quad\, y(t) = `` $(latexify(yₑₓ))
	
	``t`` $(@bind t Slider(-10:0.1:10, default = 0, show_value = true))
	"""
end

# ╔═╡ e533e6a0-2d8b-11eb-2a42-0fc2f7c7c276
begin
	plotlyjs(size = (675, 675))
	plot_chain_rule(x, y, dx, dy, t; col = plot_cols[6])
end

# ╔═╡ 744d6c50-2a0b-11eb-02c8-bd9f527a8d42
md"""
## Gradients

The **gradient** of a multivariate function ``f\left(x,y\right)`` at a point ``\left(x_0, y_0\right)`` is defined as a **vector** of the *partial derivatives* of ``f``:

```math
\left. \vec{\nabla} f \right\vert_{\left(x_0, y_0\right)} = \left({\left. \frac{\partial f}{\partial x} \right\vert_{\left(x_0, y_0\right)}, \left. \frac{\partial f}{\partial y} \right\vert_{\left(x_0, y_0\right)}}\right).
```

**Note:** The notation ``\nabla`` is pronounced *"nabla"* (no, it is not a greek letter). In MS Word or ``\LaTeX`` use `\nabla` to produce it. 

The *directional derivative* of ``f\left(x,y\right)`` with respect to a unit vector ``u`` at a point ``\left(x_0, y_0\right)`` can be expressed using the *gradient* of ``f`` at ``\left(x_0, y_0\right)``: 

```math
\left. D_{\vec{u}}f \right\vert_{\left(x_0, y_0\right)} = \left. \vec{\nabla} f \right\vert_{\left(x_0, y_0\right)} \cdot \vec{u},
```

where ``\cdot`` is the standard inner product over ``\mathbb{R}^2`` (*i.e.* ``\vec{a\vphantom{b}}\cdot\vec{b}=a_1 b_1 + a_2 b_2``). 
\
This result allows us to gain a geometrical perspective about the *gradient*. Note that ``\vec{a\vphantom{b}}\cdot\vec{b} = \left\Vert \vec{a\vphantom{b}} \right\Vert \left\Vert \vec{b} \right\Vert \cos{\alpha}`` where ``\alpha`` is the angle between the vectors ``\vec{a\vphantom{b}}`` and ``\vec{b}``. We defined ``\vec{u}`` to be a unit vector, so: 

```math
\left. D_{\vec{u}}f \right\vert_{\left(x_0, y_0\right)} = \left. \vec{\nabla} f \right\vert_{\left(x_0, y_0\right)} \cdot \vec{u} = \left\Vert \left. \vec{\nabla} f \right\vert_{\left(x_0, y_0\right)} \right\Vert \left\Vert \vec{u} \vphantom{\left. \vec{\nabla} f \right\vert_{\left(x_0, y_0\right)}}\right\Vert \cos{\alpha} = \left\Vert \left. \vec{\nabla} f \right\vert_{\left(x_0, y_0\right)} \right\Vert \cos{\alpha}.
```

From this we can gain several insights: 

- The **maximal value** of the *directional derivative* in the direction of ``\vec{u}`` is obtained *if and only if* ``\vec{u}`` is at **the same direction as the _gradient_**, since then (and only then) we will get ``\cos{\alpha} = 1``. 

    - In this case the *directional derivative* will be *nonnegative*, thus the function will usualy **increase in the direction of the gradient**. 

- The **minimal value** of the *directional derivative* in the direction of ``\vec{u}`` is obtained *if and only if* ``\vec{u}`` is at **the opposite direction of the _gradient_**, since then (and only then) we will get ``\cos{\alpha} = -1``. 

    - In this case the *directional derivative* will be *nonpositive*, thus the function will usualy **decrease in the direction opposite to that of the gradient**.

- The *gradient* is **perpendicular** the the level curves (curves of constant height ``f\left(x,y\right) = C``, see contour plots), because the *directional derivative* in their direction must be ``0``, and that is possible *if and only if* ``\cos{\alpha} = 0``. 

Let's have a look at the gradient! Choose a point and a direction:
\
``x_0`` $(@bind x⁰ Slider(-box_lim:0.1:box_lim, default = x₀, show_value = true))
``\qquad\qquad\quad``
``\alpha`` $(@bind αₓ Slider(-1:(1/16):1, default = rand(-(1 - (1/16)):(1/16):(1 - (1/16))), show_value = true))π
\
``y_0`` $(@bind y⁰ Slider(-box_lim:0.1:box_lim, default = y₀, show_value = true))
"""

# ╔═╡ 3c6c8980-2a6d-11eb-2257-bb57b19a9d1f
begin 
	α = αₓ*π
	∇ = ∇f([x⁰, y⁰])
	θ⁰ = α + atan(∇[2], ∇[1])
	
	md"""
	``\left .\vec{∇}f \right\vert_{\left(x_0, y_0\right)} = (`` $(latexify(∇[1])) ``, `` $(latexify(∇[2])) ``)``
	``\qquad\quad\!\!``
	``\vec{u} = (`` $(latexify(cos(θ⁰))) ``, `` $(latexify(sin(θ⁰))) ``)``
	
	**Note:** The displayed *gradient* vector in the contour plot is *rescaled to a unit vector* for better visualization. 
	"""
end

# ╔═╡ 32e938d0-2a69-11eb-017c-6b53e7d78c2e
begin
	# Plot directional derivative
	plotlyjs(size = (675, 675))
	plot_gradient(x⁰, y⁰, θ⁰, col_u = plot_cols[3], col_∇ = plot_cols[7])
end

# ╔═╡ a89bad10-341e-11eb-31b8-2d7460800706
begin
	# FUTURE: TOC should be a part of PlutoUI in the near future. Once that happens,  remove all but last line
	
	import Markdown: withtag
	
	"""Generate Table of Contents using Markdown cells. Headers h1-h6 are used. 
	`title` header to this element, defaults to "Table of Contents"
	`indent` flag indicating whether to vertically align elements by hierarchy
	`depth` value to limit the header elements, should be in range 1 to 6 (default = 3)
	`aside` fix the element to right of page, defaults to true
	# Examples:
	`TableOfContents()`
	`TableOfContents("Experiments 🔬")`
	`TableOfContents("📚 Table of Contents", true, 4, true)`
	"""
	struct TableOfContents
		title::AbstractString
		indent::Bool
		depth::Int
		aside::Bool
	end
	TableOfContents(title::AbstractString; indent::Bool=true, depth::Int=3, aside::Bool=true) = TableOfContents(title, indent, depth, aside)
	TableOfContents() = TableOfContents("Table of Contents", true, 3, true)

	function Base.show(io::IO, ::MIME"text/html", toc::TableOfContents)

		if toc.title === nothing || toc.title === missing 
			toc.title = ""
		end

		withtag(io, :script) do
			print(io, """
				if (document.getElementById("toc") !== null){
					return html`<div>TableOfContents already added. Cannot add another.</div>`
				}
				const getParentCellId = el => {
					// Traverse up the DOM tree until you reach a pluto-cell
					while (el.nodeName != 'PLUTO-CELL') {
						el = el.parentNode;
						if (!el) return null;
					}
					return el.id;
				}     
				const getElementsByNodename = nodeName => Array.from(
					document.querySelectorAll(
						"pluto-notebook pluto-output " + nodeName
					)
				).map(el => {
					return {
						"nodeName" : el.nodeName,
						"parentCellId": getParentCellId(el),
						"innerText": el.innerText
					}
				})

				const getPlutoCellIds = () => Array.from(
					document.querySelectorAll(
						"pluto-notebook pluto-cell"
					)
				).map(el => el.id)

				const isSelf = el => {
					try {
						return el.childNodes[1].id === "toc"
					} catch {                    
					}
					return false
				}            
				const getHeaders = () => {
					const depth = Math.max(1, Math.min(6, $(toc.depth))) // should be in range 1:6
					const range = Array.from({length: depth}, (x, i) => i+1) // [1, ... depth]
					let headers = [].concat.apply([], range.map(i => getElementsByNodename("h"+i))); // flatten [[h1s...], [h2s...], ...]
					const plutoCellIds = getPlutoCellIds()
					headers.sort((a,b) => plutoCellIds.indexOf(a.parentCellId) - plutoCellIds.indexOf(b.parentCellId)); // sort in the order of appearance
					return headers
				}
				const tocIndentClass = '$(toc.indent ? "-indent" : "")'
				const render = (el) => `\${el.map(h => `<div class="toc-row">
												<a class="\${h.nodeName}\${tocIndentClass}" 
													href="#\${h.parentCellId}" 
													onmouseover="(()=>{document.getElementById('\${h.parentCellId}').firstElementChild.classList.add('highlight-pluto-cell-shoulder')})()" 
													onmouseout="(()=>{document.getElementById('\${h.parentCellId}').firstElementChild.classList.remove('highlight-pluto-cell-shoulder')})()"
													onclick="((e)=>{
														e.preventDefault();
														document.getElementById('\${h.parentCellId}').scrollIntoView({
															behavior: 'smooth', 
															block: 'center'
														});
													})(event)"
													> \${h.innerText}</a>
											</div>`).join('')}`
				const updateCallback = e => {
					if (isSelf(e.detail.cell_id)) return
					document.getElementById('toc-content').innerHTML = render(getHeaders())                
				}
				window.addEventListener('cell_output_changed', updateCallback)
				const tocClass = '$(toc.aside ? "toc-aside" : "")'
				return html`<div class=\${tocClass} id="toc">
								<div class="markdown">
									<p class="toc-title">$(toc.title)</p>
									<div class="toc-content" id="toc-content">
											\${render(getHeaders())}
									</div>
								</div>
							</div>`
			""")
		end

		withtag(io, :style) do        
			print(io, """
				@media screen and (min-width: 1081px) {
					.toc-aside {
						position:fixed; 
						right: 1rem;
						top: 5rem; 
						width:25%; 
						padding: 10px;
						border: 3px solid rgba(0, 0, 0, 0.15);
						border-radius: 10px;
						box-shadow: 0 0 11px 0px #00000010;
						max-height: 500px;
						overflow: auto;
					}
				}    
				.toc-title{
					display: block;
					font-size: 1.5em;
					margin-top: 0.67em;
					margin-bottom: 0.67em;
					margin-left: 0;
					margin-right: 0;
					font-weight: bold;
					border-bottom: 2px solid rgba(0, 0, 0, 0.15);
				}
				.toc-row {
					white-space: nowrap;
					overflow: hidden;
					text-overflow: ellipsis;
					padding-bottom: 2px;
				}
				.highlight-pluto-cell-shoulder {
					background: rgba(0, 0, 0, 0.05);
					background-clip: padding-box;
				}
				a {
					text-decoration: none;
					font-weight: normal;
					color: gray;
				}
				a:hover {
					color: black;
				}
				a.H1-indent {
					padding: 0px 0px;
				}
				a.H2-indent {
					padding: 0px 10px;
				}
				a.H3-indent {
					padding: 0px 20px;
				}
				a.H4-indent {
					padding: 0px 30px;
				}
				a.H5-indent {
					padding: 0px 40px;
				}
				a.H6-indent {
					padding: 0px 50px;
				}
				""")
		end
	end

	get(toc::TableOfContents) = toc.default
	
	TableOfContents("Table of Contents")
end

# ╔═╡ Cell order:
# ╟─6fa8a130-29af-11eb-0f35-776973ddd076
# ╟─75a0a510-29af-11eb-2005-dfb6db948065
# ╟─3fab8fb0-29f0-11eb-058e-b1717bab59bd
# ╟─ed15ed70-29ec-11eb-380f-139b9e82f0f6
# ╟─859be820-29ec-11eb-2053-97a2830c8a8d
# ╟─868e4b20-29af-11eb-289e-ed49c852493d
# ╟─33e9bc80-29f3-11eb-11a7-b1cd60c8ae03
# ╟─cb617850-29ee-11eb-2e28-a99ff45f2ed1
# ╟─a1504642-29ee-11eb-25e9-97312b012e13
# ╟─b35b9120-29c9-11eb-1e50-8f2c7d0a635d
# ╟─ccd642b0-2d8b-11eb-150a-f1ce567ecee1
# ╟─e7612a0e-2d8f-11eb-3e0b-8da86e8dc20d
# ╟─e533e6a0-2d8b-11eb-2a42-0fc2f7c7c276
# ╟─744d6c50-2a0b-11eb-02c8-bd9f527a8d42
# ╟─3c6c8980-2a6d-11eb-2257-bb57b19a9d1f
# ╟─32e938d0-2a69-11eb-017c-6b53e7d78c2e
# ╟─08778ee0-2a77-11eb-296c-453f1932bf58
# ╟─6990bc60-2a77-11eb-2415-b38966242435
# ╟─a89bad10-341e-11eb-31b8-2d7460800706
