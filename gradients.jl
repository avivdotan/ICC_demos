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
	
	# Toy example
	# fₑₓ = :(x^2 + 3x*y + y - 1)
	
	md"""
	## Partial Derivatives, Directional Derivatives and Gradients
	
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
	
	**Note:** ``z``-axiz is automatically  *shifted* and *rescaled* for display purposes.
	
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

# ╔═╡ 859be820-29ec-11eb-2053-97a2830c8a8d
md"""
## Partial derivatives

The **partial derivative** of the multivariable function ``f\left(x,y\right)`` with respect to ``x`` at the point ``\left(x_0, y_0\right)`` is defined by keeping ``y`` constant at ``y_0``, thus effectively turning ``f\left(x, y\right)`` into a *single variable* function ``f\left(x, y_0\right)``: 

```math
\left. \frac{\partial f}{\partial x}\right\vert_{\left(x_0, y_0\right)} = \left. \frac{d}{dx} f\left(x, y_0\right) \right|_{x = x_0}= \lim_{h\to 0} {\frac{f\left(x_0 + h, y_0\right) - f\left(x_0, y_0\right)}{h}}.
```

**Note:** The notation ``\partial`` is pronounced *"dee"* (same as a regular ``\mathrm{d}``), but in mathematics it is used for different purposes.  We will use ``\mathrm{d}`` for regular derivatives and ``\partial`` for partial derivatives. In MS Word or ``\LaTeX`` use `\partial` to produce it. 

Similarly, the partial derivative of ``f\left(x,y\right)`` with respect to ``y`` at the point ``\left(x_0, y_0\right)`` is defined by keeping ``x`` constant at ``x_0``: 

```math
\left. \frac{\partial f}{\partial y}\right\vert_{\left(x_0, y_0\right)} = \left. \frac{d}{dy} f\left(x_0, y\right) \right|_{y = y_0}= \lim_{h\to 0} {\frac{f\left(x_0, y_0 + h\right) - f\left(x_0, y_0\right)}{h}}.
```

Let's have a look at these derivatives! 
Choose a point: 
\
``x_0`` $(@bind x₀ Slider(-box_lim:0.1:box_lim, default = round(rand((-box_lim/2):0.1:(box_lim/2)), digits = 1), show_value = true))
\
``y_0`` $(@bind y₀ Slider(-box_lim:0.1:box_lim, default = round(rand((-box_lim/2):0.1:(box_lim/2)), digits = 1), show_value = true))
"""

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
	
	function plot_3d_plane!(p3d, x₀, y₀, u, s, x₁, x₂, y₁, y₂, col)

		# 3D plot
		plot!(x₀ .+ s*u[1], y₀ .+ s*u[2], 
			f.([[x₀ .+ sᵢ*u[1], y₀ .+ sᵢ*u[2]] for sᵢ ∈ s]), 
			linetype = :path3d, linewidth = 3,
			color = col, 
			label = "")
		lin_app = f([x₀, y₀]) .+ ∇f([x₀, y₀])'*u*s
		ind = (0 .<= lin_app .<= z_lim)
		plot!(p3d, x₀ .+ s[ind]*u[1], y₀ .+ s[ind]*u[2], 
			lin_app[ind],
			seriescolor = col, linewidth = 2, linestyle = :dash,
			label = "")
		plot!([x₁, x₁, x₂, x₂, x₁],
			  [y₁, y₁, y₂, y₂, y₁],
			  [0, z_lim, z_lim, 0, 0], 
			seriestype = :path3d,
			seriescolor = col, linewidth = 3, fill = 0,
			label = "")
		
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
			linecolor = col, linewidth = 1.5, linestyle = :dash,
			label = "")
		quiver!(p2d, [x₀], [y₀], quiver = [(3u..., )],
			linecolor = col, linewidth = 5)
		
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
			seriescolor = col, linewidth = 1.5,
			ylim = (0, z_lim),
			xlabel = "s", ylabel = "f(x₀ + su)",
			label = "")
		plot!(p1d, s, f([x₀, y₀]) .+ ∇f([x₀, y₀])'*u*(s .- s₀),
			seriescolor = col, linewidth = 1, linestyle = :dash,
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
	
end;

# ╔═╡ 33e9bc80-29f3-11eb-11a7-b1cd60c8ae03
begin
	# Plot directional derivative
	plotlyjs(size = (675, 675))
	plot_partials(x₀, y₀, col_x = plot_cols[1], col_y = plot_cols[2])
end

# ╔═╡ cb617850-29ee-11eb-2e28-a99ff45f2ed1
md"""
## Directional derivative

The **directional derivative** of ``f\left(x,y\right)`` with respect to a *unit vector* ``u`` at a point ``\left(x_0, y_0\right)`` is defined by keeping the function constant in any direction orthogonal to ``u``, thus essentially turning ``f\left(x, y\right)`` into a *single variable* function (similarly to *partial derivatives*): 

```math
\left. D_{\vec{u}}f \right\vert_{\left(x_0, y_0\right)} = \left. \frac{\partial f}{\partial s}\right\vert_{\vec{u}, \left(x_0, y_0\right)} = \lim_{s\to 0} {\frac{f\left(x_0 + s\cdot u_1, y_0 + s\cdot u_2\right) - f\left(x_0, y_0\right)}{s}}\qquad \text{where } \left\Vert u\right\Vert = 1.
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

- The *gradient* is **perpendicular** the the constant height lines (see contour plots), because the *directional derivative* in their direction must be ``0``, and that is possible *if and only if* ``\cos{\alpha} = 0``. 

Let's have a look at the gradient! Choose a point and a direction:
\
``x_0`` $(@bind x⁰ Slider(-box_lim:0.1:box_lim, default = x₀, show_value = true))
``\qquad\qquad\quad``
``\alpha`` $(@bind αₓ Slider(-1:(1/16):1, default = rand(-1:(1/16):1), show_value = true))π
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
	``\vec{u} = (`` $(latexify(cos(θ))) ``, `` $(latexify(sin(θ))) ``)``
	
	**Note:** The displayed *gradient* vector in the contour plot is *rescaled to a unit vector* for better visualization. 
	"""
end

# ╔═╡ 32e938d0-2a69-11eb-017c-6b53e7d78c2e
begin
	# Plot directional derivative
	plotlyjs(size = (675, 675))
	plot_gradient(x⁰, y⁰, θ⁰, col_u = plot_cols[3], col_∇ = plot_cols[7])
end

# ╔═╡ 08778ee0-2a77-11eb-296c-453f1932bf58
md"""
Now we can have a look at the gradients of ``f(x,y)`` at many different points simultaneously: 

**Note:** The displayed *gradient* vectors are *rescaled* for better visualization. 
"""

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
		linecolor = plot_cols[7])
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
# ╟─744d6c50-2a0b-11eb-02c8-bd9f527a8d42
# ╟─3c6c8980-2a6d-11eb-2257-bb57b19a9d1f
# ╟─32e938d0-2a69-11eb-017c-6b53e7d78c2e
# ╟─08778ee0-2a77-11eb-296c-453f1932bf58
# ╟─6990bc60-2a77-11eb-2415-b38966242435
