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

# ╔═╡ 2824ff00-308c-11eb-3ea2-857fe96ce19d
begin
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add("Latexify")
	Pkg.add("LatexPrint")
	Pkg.add("Distributions")
	Pkg.add("Plots")
	Pkg.add("PlotlyJS")
	using PlutoUI
	using Latexify
	using LatexPrint
	using Distributions
	using Random
	using LinearAlgebra
	using Plots
	
	# Latexify default
	set_default(cdot = false, fmt = FancyNumberFormatter(3))
	
	# Plots default palette
	plot_cols = palette(:default)
	
	md"""
	# Linear Perceptron - Demo
	
	$(@bind redraw_all Button("Redraw all samples"))
	"""
end

# ╔═╡ 555ec140-308c-11eb-05ca-878451b80eb5
begin
	redraw_all
	
	P_max = 250
	N = 2
	box_lim = 10
	noise_std = 6
	
	w₀ = [3, 2]
	Xᶠ = rand(Uniform(-box_lim, box_lim), N, P_max)
	Y₀ᶠ = w₀'*Xᶠ .+ rand(Normal(0, noise_std), 1, P_max)
	
	md"""
	``P`` $(@bind P Slider(N:P_max, default = N, show_value = true))
	``\qquad``
	$(@bind redraw_chosen Button("Redraw chosen samples"))
	"""
end

# ╔═╡ a6ad7680-308d-11eb-3eeb-2506dc6c0d1e
begin
	redraw_chosen
	
	ind = randperm(P_max)[1:P]
	X = Xᶠ[:, ind]
	Y₀ = Y₀ᶠ[:, ind]
	
	C = X*X'./P
	u = X*Y₀'./P
	w = C\u
	
# 	md"""
# 	``\qquad\qquad\qquad\qquad\qquad\qquad`` $(latexify(:( C = $C, u = $u)))
	
# 	``\qquad\qquad\qquad\qquad\qquad\qquad\qquad\quad`` $(latexify(:(w^* = $w)))
# 	"""
	
	# eq_str_C_u = latexify(:(C = $C)) * " " * L"\quad, \qquad \vec{u} = %$(round.(u, sigdigits = 3))"
	# eq_str_w = L"\vec{w}^* = C^{-1}\vec{u} = %$(round.(w, sigdigits = 3))"
	
	w_str = "C \\equiv \\left\\langle\\vec{x}\\vec{x}^T\\right\\rangle = " * 
		latex_form(round.(C, sigdigits = 3)) * 
		"\\;,\\quad\\vec{u} \\equiv \\left\\langle y_0 \\vec{x}\\right\\rangle = " * 
		latex_form(round.(u, sigdigits = 3)) * 
		"\\;,\\quad\\vec{w}^* = C^{-1}\\vec{u}" 
		# latex_form(round.(w, sigdigits = 3))
	
	latexify(:($w_str = $w))
	
end

# ╔═╡ c8adf350-308c-11eb-32dd-ef8accbca31b
begin
	plotlyjs(size = (670, 500), extra_plot_kwargs = KW(:include_mathjax => "cdn"))
	
	# 3D plot of samples and hyperplane
	gx₁ = -box_lim:0.1:box_lim
	gx₂ = -box_lim:0.1:box_lim
	p3d = surface(gx₁, gx₂, (x₁, x₂) -> (w'*[x₁, x₂])[1], 
		# camera = (30, 40), 
		color = plot_cols[1],
		seriesalpha = 0.5,
		xlim = (-box_lim, box_lim), ylim = (-box_lim, box_lim), 
		xlabel = "x₁", ylabel = "x₂", zlabel = "y", 
		colorbar = :none, 
		legendfontsize = 12,
		legend = :topleft)
	scatter!(p3d, [Xᶠ[1, :]...], [Xᶠ[2, :]...], [Y₀ᶠ...], 
		markersize = 2, color = plot_cols[2], 
		label = "unchosen samples")
	scatter!(p3d, [X[1, :]...], [X[2, :]...], [Y₀...], 
		markersize = 2, color = plot_cols[1], 
		label = "chosen sampels")
	
	# Network illustration
	x_pos = [0, 2, 1]
	y_pos = [6, 6, 2]
	labels = ["\$x_1\$", "\$x_2\$", "\$y\$"]
	x_arrows = x_pos[3]*[1, 1] .- x_pos[1:2]
	y_arrows = y_pos[3]*[1, 1] .- y_pos[1:2]
	arrows_margins = 0.12
	arrow_text_pos = 0.4
	pni = scatter(x_pos, y_pos,
		series_annotations = text.(labels, 18, plot_cols[3]), 
		marker = (20, 0.2, plot_cols[3]), 
		xlim = (-1.1, 3.1), ylim = (-1, 7),
		aspect_ratio = :equal,
		framestyle = :none, 
		legendfontsize = 14,
		legend = :topleft,
		label = "")
	quiver!(pni, x_pos[1:2] .+ arrows_margins*x_arrows, 
		y_pos[1:2] .+ arrows_margins*y_arrows,
		quiver = ((1 - 2arrows_margins)*(x_pos[3]*[1, 1] .- x_pos[1:2]), 
				  (1 - 2arrows_margins)*(y_pos[3]*[1, 1] .- y_pos[1:2])), 
		seriescolor = plot_cols[4])
	annotate!(arrow_text_pos*x_pos[3] + (1 - arrow_text_pos)*x_pos[1], 
		arrow_text_pos*y_pos[3] + (1 - arrow_text_pos)*y_pos[1], 
		text("\$w_1\$", 18, :top, :right, plot_cols[4]))
	annotate!(arrow_text_pos*x_pos[3] + (1 - arrow_text_pos)*x_pos[2], 
		arrow_text_pos*y_pos[3] + (1 - arrow_text_pos)*y_pos[2], 
		text("\$w_2\$", 18, :top, :left, plot_cols[4]))
	# annotate!(-1.1, 0, 
	# 	text("\$y = \\vec{w}\\cdot\\vec{x}\$", 
	# 		12, :bottom, :left, :black))
	# annotate!(-0.78, -0.5, 
	# 	text("\$ = w_1x_1+w_2x_2\$", 
	# 		12, :middle, :left, :black))
	# annotate!(-0.997, -1, 
	# 	text("\$ \\underbrace{=}_{\\vec{w}=\\vec{w}^*} " * 
	# 		latexify(:($(w[1])*x₁ + $(w[2])*x₂), env = :raw) * "\$", 
	# 		# latex_form(round.(w[1], sigdigits = 3)) * "x_1 + " *
	# 		# latex_form(round.(w[2], sigdigits = 3)) * "x_2\$", 
	# 		12, :top, :left, :black))
	y_str = "\$ \\begin{array} \\\\" * 
			"y &= \\vec{w}\\cdot\\vec{x} \\\\" * 
			" &= w_1x_1+w_2x_2 \\\\" *
			" &\\!\\!\\!\\!\\!\\; \\underbrace{=}_{\\vec{w}=\\vec{w}^*} " * 
			latexify(:($(w[1])*x₁ + $(w[2])*x₂), env = :raw) * 
			"\\\\ \\end{array} \$"
	annotate!(3, -0.5, 
		text(y_str, 
			12, :middle, :right, :black))
	
	plot(p3d, pni, layout = grid(1, 2, widths = [0.7, 0.3]))
end

# ╔═╡ e8f51ad0-30a0-11eb-2cc1-5d3a4ed9ddf4
begin
	Y = w'*X
	εₜᵣ = mean((Y₀ - Y).^2)/2
	if abs(εₜᵣ - 0) < eps()
		εₜᵣ = 0
	end
	
	Cₐ = box_lim^2/3*Diagonal(ones(N))
	uₐ = (box_lim^2/3)*w₀
	#FUTURE:  εᵍ should be just ((w - w₀)'*Cₐ*(w - w₀) .+ noise_std^2)[1]/2
	εᵍ = ((w - w₀)'*Cₐ*(w - w₀) .- uₐ'*inv(Cₐ)*uₐ .+ w₀'*Cₐ*w₀ .+ noise_std^2)[1]/2
	
	md"""
	## Training Error vs. Generalization Error
	
	Training error: ``\qquad\qquad\varepsilon_{\mathrm{tr}}\left(\vec{w}\right) = \frac{1}{2P}\sum_{\mu = 1}^P \left(y_0^{\mu} - \vec{w}\cdot \vec{x}^{\mu}\right)^2 \qquad\ \: \!\! \underbrace{=}_{\vec{w}=\vec{w}^*} `` $(latexify(εₜᵣ))
	\
	Generalization error: ``\quad\:\:\,\varepsilon_{\mathrm{g}}\left(\vec{w}\right) \,= \frac{1}{2}\int_{\mathbb{R}^N} \left(y_0\left(\vec{x}\right) - \vec{w}\cdot \vec{x}\right)^2 p\left[\vec{x}\right]\mathrm{d}\vec{x} \,\underbrace{=}_{\vec{w}=\vec{w}^*} `` $(latexify(εᵍ))
	"""
end

# ╔═╡ 3ea0f380-31b1-11eb-2d80-598995e97f75
begin
	gr(size = (670, 300))
	
	εᵗʳ_ave(p) = (noise_std^2/2)*max(1 - N/p, 0)
	εᵍ_ave(p) = (noise_std^2/2)*(1 + N/p)
	
	psᵍ = 1:P_max
	psᵗʳ = 1:P_max
	plot(psᵍ, εᵍ_ave, 
		linecolor = plot_cols[1],
		legendfontsize = 14,
		label = "\$\\left\\langle\\varepsilon_{\\mathrm{g}}\\right\\rangle\$")
	plot!(psᵗʳ, εᵗʳ_ave, 
		linecolor = plot_cols[2],
		label = "\$\\left\\langle\\varepsilon_{\\mathrm{tr}}\\right\\rangle\$")
	scatter!([P], [εᵍ], 
		markercolor = plot_cols[1], 
		label = "\$\\varepsilon_{\\mathrm{g}}\$")
	scatter!([P], [εₜᵣ], 
		markercolor = plot_cols[2], 
		label = "\$\\varepsilon_{\\mathrm{tr}}\$")
	xlabel!("\$P\$")
	ylabel!("Error")
	
end

# ╔═╡ Cell order:
# ╟─2824ff00-308c-11eb-3ea2-857fe96ce19d
# ╟─555ec140-308c-11eb-05ca-878451b80eb5
# ╟─a6ad7680-308d-11eb-3eeb-2506dc6c0d1e
# ╟─c8adf350-308c-11eb-32dd-ef8accbca31b
# ╟─e8f51ad0-30a0-11eb-2cc1-5d3a4ed9ddf4
# ╟─3ea0f380-31b1-11eb-2d80-598995e97f75
