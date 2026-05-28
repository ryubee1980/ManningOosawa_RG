using QuadGK
using CairoMakie

############################################################
# RG spherical
############################################################

function I_sphere(Q, ρ; ρL)
    val, err = quadgk(s -> s^2 * exp(Q / s), ρ, ρL)
    return val
end

function beta_sphere(Q, ρ; L=10000.0, lB=7.0)
    ρL = L / lB

    if ρ >= ρL
        return 0.0
    end

    I = I_sphere(Q, ρ; ρL=ρL)

    return -ρ^2 * Q * exp(Q / ρ) / I
end

function beta_vs_Q(Q_set, ρ; L=10000.0, lB=7.0)
    β = zeros(length(Q_set))
    for (i, Q) in enumerate(Q_set)
        β[i] = beta_sphere(Q, ρ; L=L, lB=lB)
    end
    return β
end

function Q_R(Q; ρ0=1.0, L=10000.0, lB=7.0, dρ=1e-3)
    ρL = L / lB

    ρs = ρ0:dρ:ρL
    QR = zeros(length(ρs))
    QR[1] = Q

    for i in 2:length(ρs)
        β = beta_sphere(QR[i-1], ρs[i-1]; L=L, lB=lB)
        QR[i] = QR[i-1] + β * dρ

        # 数値誤差で負にならないようにする
        if QR[i] < 0
            QR[i] = 0.0
        end
    end

    return ρs, QR
end

############################################################
# Plotting spherical RG
############################################################

function plot_RG_spherical(ρset, Q_set; L=10000.0, lB=7.0, ρ0=1.0,
                           length_Q=1000, dρ=1e-3)

    ρL = L / lB

    f = Figure(fonts = (; regular = "CMU Serif"), size = (1000, 500))
    Label(f[1,1:2], L"(b)RG analysis for spherical geometry ($\rho_a=a/\ell_B=1.5$, $\rho_L=$%$(round(ρL, digits=2)))", fontsize=28)

    ax = Axis(
        f[2,1],
        ylabel = L"\beta(Q_R;\rho_L,\rho)",
        xlabel = L"Q_R",
        title = L"(b-1) $\beta$-function"
    )

    ax.xlabelsize = 22
    ax.ylabelsize = 22
    ax.xticklabelsize = 18
    ax.yticklabelsize = 18
    ax.titlesize = 25

    ax2 = Axis(
        f[2,2],
        ylabel = L"Q_R(\rho)",
        xlabel = L"\rho",
        title = L"(b-2) RG flow of $Q_R$"
    )

    ax2.xlabelsize = 22
    ax2.ylabelsize = 22
    ax2.xticklabelsize = 18
    ax2.yticklabelsize = 18
    ax2.titlesize = 25

    ############################################################
    # inset (log scale)
    ############################################################

    axins = Axis(
        f[2,2],
        width = Relative(0.35),
        height = Relative(0.35),
        halign = 0.15,
        valign = 0.65,
        xscale = log10,
        backgroundcolor = :white
    )

    axins.xlabelsize = 14
    axins.ylabelsize = 14
    axins.xticklabelsize = 12
    axins.yticklabelsize = 12
    axins.titlesize = 14

    axins.title = "log-scale inset"

    for Q in Q_set
        ρs, QR = Q_R(Q; ρ0=ρ0, L=L, lB=lB, dρ=dρ)

        lines!(
            axins,
            ρs,
            QR,
            label = L"$Q=$%$(Q)"
        )
    end

    # inset の表示範囲
    xlims!(axins, 1, 30)

    Q_range = range(0.0, 25.0, length=length_Q)

    for ρ in ρset
        β = beta_vs_Q(Q_range, ρ; L=L, lB=lB)
        lines!(ax, Q_range, β, label=L"$\rho=$%$(ρ)")
    end

    for Q in Q_set
        ρs, QR = Q_R(Q; ρ0=ρ0, L=L, lB=lB, dρ=dρ)
        lines!(ax2, ρs, QR, label=L"$Q=$%$(Q)")
    end


    

    xlims!(ax, -0.1, 25.0)
    xlims!(ax2, ρ0, min(ρL, 1430.0))

    axislegend(ax, position=(0.0,0.0), bboxcolor=:transparent, labelsize=22)
    axislegend(ax2, position=(1.0,1.0), bboxcolor=:transparent, labelsize=22)

    #hlines!(ax, [0.0], linestyle=:dash, color=:black)

    save("RG_spherical.pdf", f)

    return f
end

ρset = [1.0, 3.0, 10.0, 100.0]
Q_set = [0.5, 1.0, 2.0, 5.0, 10.0]

plot_RG_spherical([1.5,1.8, 300.0,1000,1400], [30.0,10,5]; L=10000.0, lB=7.0, ρ0=1.5, dρ=1e-3)