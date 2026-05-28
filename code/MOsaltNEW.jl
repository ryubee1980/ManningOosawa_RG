using CairoMakie

function beta_inner(x)
    if x > 1.0
        return -2.0 * x * (x - 1.0)
    else
        return 0.0
    end
end

function integrate_to_tD(xi; ta=log(3.0/7.0), tD=7.6, nstep=5000)
    if tD <= ta
        return xi
    end

    h = (tD - ta) / nstep
    x = xi

    for _ in 1:nstep
        k1 = beta_inner(x)
        k2 = beta_inner(x + 0.5h*k1)
        k3 = beta_inner(x + 0.5h*k2)
        k4 = beta_inner(x + h*k3)

        x += h * (k1 + 2k2 + 2k3 + k4) / 6
    end

    return x
end

function plot_xeff_curves()
    ta = log(3.0 / 7.0)

    xis = range(0.5, 6.0, length=80)
    tDs = [7.6, 1.8, 0.7, 0.3]

    fig = Figure(size=(700, 500), fonts=(; regular="CMU Serif"))
    ax = Axis(fig[1, 1],
        xlabel=L"\xi",
        ylabel=L"\xi_{\mathrm{eff}}",
        title=L"RG estimate: $\xi_{\mathrm{eff}}=\xi_R(t_{\mathrm{Debye}})$"
    )
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18
    ax.titlesize=28

    for tD in tDs
        xeffs = [integrate_to_tD(xi; ta=ta, tD=tD) for xi in xis]
        scatterlines!(ax, xis, xeffs,
            label=L"t_{\mathrm{Debye}} = %$(tD)")
    end

    vlines!(ax, [1.0], linestyle=:dash, color=:black)
    axislegend(ax, position=(1.0,0.0),bboxcolor=:transparent, labelsize=22)

    save("xi_eff_RG_salt_cutoff.pdf", fig)
    return fig
end

fig = plot_xeff_curves()