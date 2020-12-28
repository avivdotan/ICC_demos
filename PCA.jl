### A Pluto.jl notebook ###
# v0.12.18

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

# ╔═╡ a05e9050-4426-11eb-0c69-79b297463ea2
begin
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add("Distributions")
	Pkg.add("Latexify")
	Pkg.add("Plots")
	Pkg.add("PlotlyJS")
	using PlutoUI
	using Distributions
	using LinearAlgebra
	using Latexify
	using Plots
	
	# Latexify default
	set_default(cdot = false, fmt = FancyNumberFormatter(3))
	
	# Plots default palette
	plot_cols = palette(:default)
	
	md"""
	# Principal Component Analysis (PCA) - Demo
	"""
end

# ╔═╡ c461c17e-4452-11eb-29dc-1399eb26a0c0
md"""
Input dimension: $(@bind N_str Select(["2" => "2D", "3" => "3D"]))
"""

# ╔═╡ f3540d30-4426-11eb-0983-d74a034c19fa
md"""
\# of samples: $(@bind P Slider(10:10:300, default = 100, show_value = true))
"""

# ╔═╡ 4a29cd80-4453-11eb-0786-3790d6ae767e
md"""
Input distribution: $(@bind dist Select(["gauss" => "Normal", "ellipse" => "Ellipsoid"], default = "ellipse"))
$(@bind redraw_samples Button("Redraw samples"))
"""

# ╔═╡ 6bfef7ee-447b-11eb-1001-e3ba0a3683c9
md"""
$(@bind redraw_sample Button("Redraw sample point"))
"""

# ╔═╡ fbca5380-4452-11eb-1542-956601ca296c
N = Meta.parse(N_str);

# ╔═╡ ee549620-4452-11eb-3204-0f0713dbf55e
if N > 2
	md"""
	Output dimension: $(@bind M_str Select(["1" => "1D", "2" => "2D"]))
	"""
else
	M_str = "1"
	md""
end

# ╔═╡ 31ec2290-4430-11eb-1a19-1fd62361adef
begin
	redraw_samples
	
	μ⁰ = 2*rand(N) .- 1
	
	if dist == "gauss"
		Σ⁰ = randn(N, N)*(2/π)
		Σ⁰ = Σ⁰*Σ⁰'
		X = rand(MvNormal(μ⁰, Σ⁰), P)
	
	elseif dist == "ellipse"
		# Ellipsoid (approximately uniform)
		v₁ = rand(N)
		v₁ /= norm(v₁)
		v₂ = [-v₁[2], v₁[1], zeros(N - 2)...]
		a, b = 1 + rand(), 0.3 + 1.3rand()
		r = rand(P).^(1/N)
		θ = 2π*rand(P)
		if N == 2
			X = a*v₁*(r.*cos.(θ))' .+ b*v₂*(r.*sin.(θ))'
		elseif N == 3
			v₃ = v₁ × v₂
			c = 0.3 + 1rand()
			φ = π*(2rand(P) .- 1)
			X = a*v₁*(r.*cos.(θ).*cos.(φ))' .+ b*v₂*(r.*sin.(θ).*cos.(φ))' .+ c*v₃*(r.*sin.(φ))'
		end
		X .+= μ⁰
	end
	
end;

# ╔═╡ 74e7de00-442a-11eb-1083-593a38411b4a
begin
	
	function plot1d(X; kwargs...)
		D = maximum(X) - minimum(X)
		scatter(X', zeros(size(X'));
			markersize = 3,
			markerstrokewidth = 0,
			yticks = [], ylim = (-D/2, D/2),
			aspect_ratio = 0.5,
			label = "", kwargs...)
	end
	function plot1d!(plt, X; kwargs...)
		D = maximum(X) - minimum(X)
		scatter!(plt, X', zeros(size(X'));
			ylim = (-D/2, D/2),
			markersize = 3,
			markerstrokewidth = 0,
			label = "", kwargs...)
	end
	function plot2d(X; kwargs...)
		scatter(X[1, :], X[2, :];
			markersize = 3,
			markerstrokewidth = 0,
			aspect_ratio = :equal,
			label = "", kwargs...)
	end
	function plot2d!(plt, X; kwargs...)
		scatter!(plt, X[1, :], X[2, :];
			markersize = 3,
			markerstrokewidth = 0,
			aspect_ratio = :equal,
			label = "", kwargs...)
	end
	function plot3d(X; kwargs...)
		scatter(X[1, :], X[2, :], X[3, :];
			markersize = 1,
			markerstrokewidth = 0,
			aspect_ratio = :equal,
			label = "", kwargs...)
	end
	function plot3d!(plt, X; kwargs...)
		scatter!(plt, X[1, :], X[2, :], X[3, :];
			markersize = 1,
			markerstrokewidth = 0,
			aspect_ratio = :equal,
			label = "", kwargs...)
	end
	function add_line!(plt, pnt, vec; kwargs...)
		fʸˣ(x) = pnt[2] + (vec[2]/vec[1])*(x - pnt[1])
		fˣʸ(y) = pnt[1] + (vec[1]/vec[2])*(y - pnt[2])
		if length(pnt) == 3
			fᶻˣ(x) = pnt[3] + (vec[3]/vec[1])*(x - pnt[1])
			fᶻʸ(y) = pnt[3] + (vec[3]/vec[2])*(y - pnt[2])
			fˣᶻ(z) = pnt[1] + (vec[1]/vec[3])*(z - pnt[3])
			fʸᶻ(z) = pnt[2] + (vec[2]/vec[3])*(z - pnt[3])
		end
		xl = xlims(plt)
		yl = ylims(plt)
		if length(pnt) == 3
			zl = zlims(plt)
		end
		if length(pnt) == 2
			pnts = [(xl[1], fʸˣ(xl[1])),
					(xl[2], fʸˣ(xl[2])),
					(fˣʸ(yl[1]), yl[1]),
					(fˣʸ(yl[2]), yl[2])]
		elseif length(pnt) == 3
			pnts = [(xl[1], fʸˣ(xl[1]), fᶻˣ(xl[1])),
					(xl[2], fʸˣ(xl[2]), fᶻˣ(xl[2])),
					(fˣʸ(yl[1]), yl[1], fᶻʸ(yl[1])),
					(fˣʸ(yl[2]), yl[2], fᶻʸ(yl[2])),
					(fˣᶻ(zl[1]), fʸᶻ(zl[1]), zl[1]),
					(fˣᶻ(zl[2]), fʸᶻ(zl[2]), zl[2])]
		end
		filter!(p -> xl[1] ≤ p[1] ≤ xl[2], pnts)
		filter!(p -> yl[1] ≤ p[2] ≤ yl[2], pnts)
		if length(pnt) == 3
			filter!(p -> zl[1] ≤ p[3] ≤ zl[2], pnts)
		end
		plot!(plt, pnts; kwargs...)
	end
	function add_plane!(plt, pnt, vec1, vec2; kwargs...)
		# f(x, y) = pnt[3] + ((vec2[2]*(x - pnt[1]) - vec2[1]*(y - pnt[2]))/(vec1[1]*vec2[2] - vec1[2]*vec2[1]))*vec1[3] - ((vec1[2]*(x - pnt[1]) - vec1[1]*(y - pnt[2]))/(vec1[1]*vec2[2] - vec1[2]*vec2[1]))*vec2[3]
		# gx = [xl[1], xl[2]]
		# gy = [yl[1], yl[2]]
		# surface!(plt, gx, gy, f; 
		# 	seriesalpha = 0.5, 
		# 	kwargs...)
		fᶻˣʸ(x, y) = pnt[3] + ((vec2[2]*(x - pnt[1]) - vec2[1]*(y - pnt[2]))/(vec1[1]*vec2[2] - vec1[2]*vec2[1]))*vec1[3] - ((vec1[2]*(x - pnt[1]) - vec1[1]*(y - pnt[2]))/(vec1[1]*vec2[2] - vec1[2]*vec2[1]))*vec2[3]
		fʸᶻˣ(z, x) = pnt[2] + ((vec2[1]*(z - pnt[3]) - vec2[3]*(x - pnt[1]))/(vec1[3]*vec2[1] - vec1[1]*vec2[3]))*vec1[2] - ((vec1[1]*(z - pnt[3]) - vec1[3]*(x - pnt[1]))/(vec1[3]*vec2[1] - vec1[1]*vec2[3]))*vec2[2]
		fˣʸᶻ(y, z) = pnt[1] + ((vec2[3]*(y - pnt[2]) - vec2[2]*(z - pnt[3]))/(vec1[2]*vec2[3] - vec1[3]*vec2[2]))*vec1[1] - ((vec1[3]*(y - pnt[2]) - vec1[2]*(z - pnt[3]))/(vec1[2]*vec2[3] - vec1[3]*vec2[2]))*vec2[1]
		xl = xlims(plt)
		yl = ylims(plt)
		zl = zlims(plt)
		pnts = [(xl[1], yl[1], fᶻˣʸ(xl[1], yl[1])),
				(xl[1], yl[2], fᶻˣʸ(xl[1], yl[2])),
				(xl[2], yl[1], fᶻˣʸ(xl[2], yl[1])),
				(xl[2], yl[2], fᶻˣʸ(xl[2], yl[2])),
				(xl[1], fʸᶻˣ(zl[1], xl[1]), zl[1]),
				(xl[1], fʸᶻˣ(zl[2], xl[1]), zl[2]),
				(xl[2], fʸᶻˣ(zl[1], xl[2]), zl[1]),
				(xl[2], fʸᶻˣ(zl[2], xl[2]), zl[2]),
				(fˣʸᶻ(yl[1], zl[1]), yl[1], zl[1]),
				(fˣʸᶻ(yl[1], zl[2]), yl[1], zl[2]),
				(fˣʸᶻ(yl[2], zl[1]), yl[2], zl[1]),
				(fˣʸᶻ(yl[2], zl[2]), yl[2], zl[2])]
		filter!(p -> xl[1] ≤ p[1] ≤ xl[2], pnts)
		filter!(p -> yl[1] ≤ p[2] ≤ yl[2], pnts)
		filter!(p -> zl[1] ≤ p[3] ≤ zl[2], pnts)
		θ(p) = acos([p...]⋅[1, 0, 0]/norm([p...]))
		sort!(pnts, by = θ)
		push!(pnts, pnts[1])
		mesh3d!(plt, pnts;
			alpha = 0.5, 
			label = "",
			kwargs...)
		xlims!(plt, xl...)
		ylims!(plt, yl...)
		zlims!(plt, zl...)
	end
	function set_2d_lims!(plt)
		xl = xlims(plt)
		yl = ylims(plt)
		lims = [min(xl[1], yl[1]), max(xl[2], yl[2])]
		xlims!(plt, lims...)
		ylims!(plt, lims...)
	end
	function set_3d_lims!(plt)
		xl = xlims(plt)
		yl = ylims(plt)
		zl = ylims(plt)
		lims = [min(xl[1], yl[1], zl[1]), max(xl[2], yl[2], zl[2])]
		xlims!(plt, lims...)
		ylims!(plt, lims...)
		zlims!(plt, lims...)
	end
	function plot_legend(n; kwargs...)
		p_legend = scatter([-1], [-1];
			xlim = (0, 1), ylim = (0, 1),
			aspect_ratio = 0.1,
			framestyle = :none,
			legend = (n == 2 ? :bottomleft : :bottomright),
			legendfontsize = 12, 
			label = "", kwargs...)
		scatter!([-1], [-1],
			markersize = 3,
			markerstrokewidth = 0,
			color = plot_cols[1],
			label = "\$X\$")
		scatter!([-1], [-1],
			markersize = 3,
			markerstrokewidth = 0,
			color = plot_cols[6],
			label = "\$Y\$")
		scatter!([-1], [-1],
			markersize = 3,
			markerstrokewidth = 0,
			color = plot_cols[2],
			label = "\$\\hat{X}\$")
		if n == 3
			scatter!([-1], [-1],
				markersize = 4,
				markerstrokewidth = 0,
				color = plot_cols[5],
				label = "\$\\left\\langle\\vec{x}\\right\\rangle\$")
			plot!([-1], [-1],
				linewidth = 3,
				color = plot_cols[3],
				label = "\$\\vec{u}^{(1)}\$")
			plot!([-1], [-1],
				linewidth = 3,
				color = plot_cols[4],
				label = "\$\\vec{u}^{(2)}\$")
			plot!([-1], [-1],
				linewidth = 3,
				color = plot_cols[7],
				label = "\$\\vec{u}^{(3)}\$")
		end
		return p_legend
	end
	function myquiver3d!(plt, from, vec; 
						 arrowhead_scale = 0.15, arrowhead_ratio = 0.4, 
						 kwargs...)

		# Get arrowhead
		to = from .+ vec
		arrow_base = from .+ (1 - arrowhead_scale)*vec
		perp_plane₁ = [vec[2], -vec[1], 0]
		perp_plane₂ = perp_plane₁ × vec
		perp_plane₁ /= norm(perp_plane₁)
		perp_plane₂ /= norm(perp_plane₂)
		θ = 0:(π/3):2π
		circ = perp_plane₁*cos.(θ)' .+ perp_plane₂*sin.(θ)'
		arrow_circ = arrow_base .+ arrowhead_scale*arrowhead_ratio*circ

		# Plot arrow line
		plot!(plt, [from[1], to[1]], [from[2], to[2]], [from[3], to[3]]; 
			label = "", kwargs...)

		# Plot arrowhead
		arrow_head = hcat([[to  arrow_circ[:, i] fill(NaN, 3)] 
						  for i = 1:length(θ)]...)
		arrow_head = [arrow_head arrow_circ]
		plot!(plt, arrow_head[1, :], arrow_head[2, :], arrow_head[3, :];
			label = "", colorbar = :none, kwargs...)
		arrow_head_vert = [(to...,), [(arrow_circ[:, i]...,) for i = 1:length(θ)]...]
		mesh3d!(plt, arrow_head_vert; 
			connections = ([0 for i = 1:length(θ)], 
						   [i for i = 1:length(θ)], 
						   [length(θ), [i for i = 1:(length(θ) - 1)]...]), 
			label = "", colorbar = :none, kwargs...)
	end
	
end;

# ╔═╡ b3940e10-4427-11eb-16e4-51d63947f127
begin
	μ = mean(X, dims = 2)
	C = cov(X')
	Cₑᵢ = eigen(C)
	eig_order = sortperm(Cₑᵢ.values, rev = true)
	U = Cₑᵢ.vectors[:, eig_order]
	λs = Cₑᵢ.values[eig_order]
end;

# ╔═╡ 3dbd5e4e-487b-11eb-0fa6-2d920db0616b
begin
	μ, C
	x_mean_text = "\\left\\langle\\vec{x}\\right\\rangle"
	cov_text = "\\quad C = \\left\\langle\\left(\\vec{x} - \\left\\langle\\vec{x}\\right\\rangle\\right)\\left(\\vec{x} - \\left\\langle\\vec{x}\\right\\rangle\\right)^T\\right\\rangle"
	latexify(:($x_mean_text = $μ, $cov_text = $C))
end

# ╔═╡ 5bcb3570-4880-11eb-1166-f93c44ca898a
begin
	U, λs
	eig_vals_text(i) 	= "\\lambda_$i"
	eig_vecs_text(i) 	= "\\vec{u}^{\\left($i\\right)}"
	separator_text(i) 	= i == 1 ? "" : "\\ , \\quad"
	eig_vals_expr = [:($(separator_text(i) * eig_vals_text(i)) = $(λs[i])) for i = 1:N]
	eig_vecs_expr = [:($(separator_text(i) * eig_vecs_text(i)) = $(U[:, i])) for i = 1:N]
	md"""
	 $(latexify(eig_vals_expr...))
	
	 $(latexify(eig_vecs_expr...))
	"""
end

# ╔═╡ 05883900-4453-11eb-29a8-b749df34e2f4
M = Meta.parse(M_str);

# ╔═╡ 9d490b90-4429-11eb-1249-c32fb09df5d7
begin
	
	Uᵣₑ = U[:, 1:M]
	
	# Encoding
	Y = Uᵣₑ'*(X .- μ)
	
	# Decoding
	X̂ = Uᵣₑ*Y .+ μ
	
	# Error
	ε = sum(Cₑᵢ.values[eig_order[(M + 1):end]])
end;

# ╔═╡ 25ad9c40-455a-11eb-3eea-a182748bcff7
md"""
Average reconstruction error: ``\varepsilon = \left\langle\left\Vert\vec{x} - \hat{x}\right\Vert^2\right\rangle = `` $(latexify(ε))
"""

# ╔═╡ a046fa90-442f-11eb-1fe7-4702ffa40c40
begin
	redraw_sample
	
	sample = rand(1:P)
	εᵐ = norm(X[:, sample] - X̂[:, sample])^2
end;

# ╔═╡ b2418170-455b-11eb-0b52-1b69d503e7f5
begin
#= md"""
Sample reconstruction error: ``\varepsilon\left(\vec{x}^{\mu}\right) = \left\Vert\vec{x}^{\mu} - \hat{x}^{\mu}\right\Vert^2 = `` $(latexify(εᵐ))
""" =#
	sample, εᵐ
	sample_error_text = "\\varepsilon\\left(\\vec{x}^{$sample}\\right) = \\left\\Vert\\vec{x}^{$sample} - \\hat{x}^{$sample}\\right\\Vert^2"
	sample_error_expr = :($sample_error_text = $εᵐ)
	md"""
	Sample reconstruction error: $(latexify(sample_error_expr))
	"""
end

# ╔═╡ dfab8580-4429-11eb-1001-774f75463983
begin
	
	if N == 2
		
		gr(size = (670, 400))
		
		p_x = plot2d(X[:, sample], 
			xlabel = "\$x_1\$", ylabel = "\$x_2\$", 
			color = :lightgreen, 
			markersize = 10)
		plot2d!(p_x, X̂[:, sample], 
			color = :lightgreen, 
			markersize = 10)
		plot!(p_x, [X[1, sample], X̂[1, sample]], [X[2, sample], X̂[2, sample]], 
			color = :lightgreen, 
			label = "")
		plot2d!(p_x, X, 
			color = plot_cols[1])
		quiver!(μ[1]*[1, 1], μ[2]*[1, 1], quiver = (U[1, :], U[2, :]),
			color = plot_cols[3:4], linewidth = 3)
		annotate!(μ[1] + U[1, 1], μ[2] + U[2, 1], 
			text("\$\\vec{u}^{(1)}\$", :top, :right, plot_cols[3]))
		annotate!(μ[1] + U[1, 2], μ[2] + U[2, 2], 
			text("\$\\vec{u}^{(2)}\$", :top, :right, plot_cols[4]))
		set_2d_lims!(p_x)
		add_line!(p_x, μ, U[:, 1], 
			color = plot_cols[3], 
			linewidth = 1, 
			label = "")
		plot2d!(p_x, X̂,  
			color = plot_cols[2])
		scatter!([μ[1]], [μ[2]],
			color = plot_cols[5], 
			markersize = 7,
			markerstrokewidth = 0,
			series_annotations = [text("\$\\left\\langle\\vec{x}\\right\\rangle\$", :left, :top, plot_cols[5])],
			label = "")
		
		l = @layout [a [b{0.8h}; c{0.1h}]]
		
	elseif N == 3
		
		plotlyjs(size = (670, 500), extra_plot_kwargs = KW(:include_mathjax => "cdn"))
		
		p_x = plot3d(X[:, sample], 
			xlabel = "\$x_1\$", ylabel = "\$x_2\$", zlabel = "\$x_3\$", 
			color = :lightgreen, 
			markersize = 4)
		plot3d!(p_x, X̂[:, sample], 
			color = :lightgreen, 
			markersize = 4)
		plot!(p_x, [X[1, sample], X̂[1, sample]], [X[2, sample], X̂[2, sample]], [X[3, sample], X̂[3, sample]], 
			color = :lightgreen, 
			linewidth = 3, 
			label = "")
		plot3d!(p_x, X, 
			color = plot_cols[1])
		myquiver3d!(p_x, μ, U[:, 1], 
			color = plot_cols[3], 
			linewidth = 5, 
			arrowhead_scale = 0.1)
		myquiver3d!(p_x, μ, U[:, 2], 
			color = plot_cols[4], 
			linewidth = 5, 
			arrowhead_scale = 0.1)
		myquiver3d!(p_x, μ, U[:, 3], 
			color = plot_cols[7], 
			linewidth = 5, 
			arrowhead_scale = 0.1)
		# annotate!(μ[1] + U[1, 1], μ[2] + U[2, 1], 
		# 	text("\$\\vec{u}^{(1)}\$", :top, :right, plot_cols[3]))
		# annotate!(μ[1] + U[1, 2], μ[2] + U[2, 2], 
		# 	text("\$\\vec{u}^{(2)}\$", :top, :right, plot_cols[4]))
		set_3d_lims!(p_x)
		if M == 1
			add_line!(p_x, μ, U[:, 1], 
				color = plot_cols[3], 
				linewidth = 3, 
				label = "")
		elseif M == 2
			add_plane!(p_x, μ, U[:, 1], U[:, 2],
				color = plot_cols[4], colorbar = :none)
			add_line!(p_x, μ, U[:, 1], 
				color = plot_cols[3], 
				linewidth = 3, 
				label = "")
			add_line!(p_x, μ, U[:, 2], 
				color = plot_cols[4], 
				linewidth = 3, 
				label = "")
		end
		plot3d!(p_x, X̂,  
			color = plot_cols[2])
		scatter!([μ[1]], [μ[2]], [μ[3]],
			color = plot_cols[5], 
			markersize = 3,
			markerstrokewidth = 0,
			label = "")
		
		l = @layout [a{0.7w} [b{0.5h}; c{0h}]]
		
	end
	
	p_l = plot_legend(N)
	
	if M == 1
		p_y = plot1d([Y[sample]], 
			xlabel = "\$y\$", 
			color = :lightgreen,
			markersize = 10)
		plot1d!(p_y, Y, 
			color = plot_cols[6])
		plot!(p_y, [xlims(p_y)...], [0, 0], 
			color = plot_cols[3], 
			linewidth = 3,
			label = "")
		scatter!(p_y, [0], [0], 
			color = plot_cols[5], 
			markersize = 7,
			markerstrokewidth = 0,
			label = "")
	elseif M == 2
		p_y = plot2d(Y[:, sample], 
			xlabel = "\$y_1\$", ylabel = "\$y_2\$", 
			color = :lightgreen,
			markersize = 10)
		plot2d!(p_y, Y, 
			color = plot_cols[6])
		set_2d_lims!(p_y)
		plot!(p_y, [xlims(p_y)...], [0, 0], 
			color = plot_cols[3], 
			linewidth = 3,
			label = "")
		plot!(p_y, [0, 0], [ylims(p_y)...], 
			color = plot_cols[4],  
			linewidth = 3,
			label = "")
		scatter!(p_y, [0], [0], 
			color = plot_cols[5], 
			markersize = 7,
			markerstrokewidth = 0,
			label = "")
	end
	
	plot(p_x, p_y, p_l, layout = l)
end

# ╔═╡ Cell order:
# ╟─a05e9050-4426-11eb-0c69-79b297463ea2
# ╟─c461c17e-4452-11eb-29dc-1399eb26a0c0
# ╟─ee549620-4452-11eb-3204-0f0713dbf55e
# ╟─f3540d30-4426-11eb-0983-d74a034c19fa
# ╟─4a29cd80-4453-11eb-0786-3790d6ae767e
# ╟─3dbd5e4e-487b-11eb-0fa6-2d920db0616b
# ╟─5bcb3570-4880-11eb-1166-f93c44ca898a
# ╟─25ad9c40-455a-11eb-3eea-a182748bcff7
# ╟─b2418170-455b-11eb-0b52-1b69d503e7f5
# ╟─6bfef7ee-447b-11eb-1001-e3ba0a3683c9
# ╟─dfab8580-4429-11eb-1001-774f75463983
# ╟─74e7de00-442a-11eb-1083-593a38411b4a
# ╟─31ec2290-4430-11eb-1a19-1fd62361adef
# ╟─b3940e10-4427-11eb-16e4-51d63947f127
# ╟─9d490b90-4429-11eb-1249-c32fb09df5d7
# ╟─a046fa90-442f-11eb-1fe7-4702ffa40c40
# ╟─fbca5380-4452-11eb-1542-956601ca296c
# ╟─05883900-4453-11eb-29a8-b749df34e2f4
