using Roots
using CairoMakie

############################################################
# xi > 1 : trigonometric (Manning condensed regime)
############################################################

function solve_gamma(xi; a=3.0, L=1000.0)

    Δ = log(L / a)

    f(γ) = γ * Δ - atan(1 / γ) - atan((xi - 1) / γ)

    γmin = 1e-12
    γmax = π / Δ - 1e-12

    find_zero(f, (γmin, γmax), Bisection())
end


############################################################
# xi < 1 : hyperbolic regime
############################################################

function solve_alpha(xi; a=3.0, L=1000.0)

    

    # Correct boundary condition:
    # α Δ = atanh(1/α) - atanh((1 - xi)/α)
     Δ = log(L / a)

    f(α) = α * Δ -
           atanh(α / (1 - xi)) +
           atanh(α)

    αmin = 1e-12
    αmax = (1 - xi) - 1e-12

    find_zero(f, (αmin, αmax), Bisection())
end


############################################################
# Main solver
############################################################

function katchalsky_solution(xi; a=3.0, L=1000.0, lB=7.0)

    ta = log(a / lB)
    tL = log(L / lB)

    Δ = log(L / a)
    xi_c = Δ / (1 + Δ)

    @show xi_c

    if xi > xi_c

        γ = solve_gamma(xi; a=a, L=L)

        θL = atan(1 / γ)
        tM = tL - θL / γ

        θ = t -> γ * (t - tM)

        u = t -> 2 * (t - tL) +
                 2 * log(cos(θ(t)) / cos(θL))

        xi_eff = t -> 1 - γ * tan(θ(t))

        return u, xi_eff

    elseif xi < xi_c

        α = solve_alpha(xi; a=a, L=L)

         yL = atanh(α)
        t0 = tL + yL / α

        y = t -> α * (t0 - t)

        u = t -> 2 * (t - tL) +
                2 * log(sinh(y(t)) / sinh(yL))

        xi_eff = t -> 1 - α * coth(y(t))

        return u, xi_eff

    else

        t0 = tL - 1.0

        u = t -> 2 * (t - tL) -
                 2 * log((t - t0) / (tL - t0))

        xi_eff = t -> 1 - 1 / (t - t0)

        return u, xi_eff
    end
end

############################################################
# RG planar
############################################################
function beta_p(Sigma,l;L=10000.0,lB=7.0)
    lL= L / lB

    -Sigma^2/(1-exp(-Sigma*(lL-l)))
    
end

function beta_vs_Sigma(Sigma_set,l;L=10000.0,lB=7.0)
    lL= L / lB
    beta=zeros(length(Sigma_set))
    for (i, Sigma) in enumerate(Sigma_set)
        beta[i] = -Sigma^2/(1-exp(-Sigma*(lL-l)))
    end
    return beta
end

function Sigma_R(Sigma;L=10000.0,lB=7.0,dl=1e-4)
    lL= L / lB
    
    ls = 0:dl:lL
    Sigma_R=zeros(length(ls))
    Sigma_R[1]=Sigma
    for i in 2:length(ls)
        Sigma_R[i] = Sigma_R[i-1] + beta_p(Sigma_R[i-1],ls[i-1];L=L,lB=lB)*dl
    end

    return ls,Sigma_R
end

############################################################
# RG analysis
############################################################

function beta(xi,t;L=1000.0,lB=7.0)
    tL= log(L / lB)

    2*xi*(xi-1)/(exp(-2*(xi-1)*(tL-t))-1)
    
end

function beta_vs_xi(xi_set,t;L=1000.0,lB=7.0)
    tL= log(L / lB)
    beta=zeros(length(xi_set))
    for (i, xi) in enumerate(xi_set)
        beta[i] = 2*xi*(xi-1)/(exp(-2*(xi-1)*(tL-t))-1)
    end
    return beta
end

function xi_R(xi;a=3.0,L=1000.0,lB=7.0,dt=1e-4)
    tL= log(L / lB)
    ta= log(a / lB)
    ts = ta:dt:tL
    xi_R=zeros(length(ts))
    xi_R[1]=xi
    for i in 2:length(ts)
        xi_R[i] = xi_R[i-1] + beta(xi_R[i-1],ts[i-1];L=L,lB=lB)*dt
    end

    return ts,xi_R
end


############################################################
# Plotting graphs for planar
############################################################

function plot_RG_planar(lset,Sig_set;L=10000.0, R=3.0, lB=7.0,length_Sigma=1000)
    lL= L / lB

    f=Figure(fonts = (; regular = "CMU Serif"),size = (1000, 500))
    Label(f[1,1:2],L"(a) RG analysis for planar geometry ($\ell_L=$%$(round(lL, digits=2)))",fontsize=28)
    
    
    ax=Axis(f[2,1], ylabel=L"\beta(\Sigma_R;\ell_L,\ell)", xlabel=L"\Sigma_R",title=L"(a-1) $\beta$-function")
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18
    ax.titlesize=25
    
    ax2=Axis(f[2,2], ylabel=L"\Sigma_R(\ell)",xlabel=L"\ell",title=L"(a-2) RG flow of $\Sigma_R$")
    ax2.xlabelsize=22
    ax2.ylabelsize=22
    ax2.xticklabelsize=18
    ax2.yticklabelsize=18
    ax2.titlesize=25

    
    
    Sigma_range = range(0, 3.5, length=length_Sigma)
    for l in lset
        beta = beta_vs_Sigma(Sigma_range, l; L=L, lB=lB)

        

        lines!(ax,Sigma_range, beta,label=L"$\ell=$%$l")
   
    end

    
    
    for Sigma in Sig_set
        ls,SigmaR_t = Sigma_R(Sigma;L=L,lB=lB,dl=1e-4)

        

        lines!(ax2,ls, SigmaR_t,label=L"$\Sigma=$%$Sigma")
    
    end
    xlims!(ax,-0.1,3.5)
    xlims!(ax2,0,15)
    #ylims!(ax,-15,1)
    #ylims!(ax2,-15,1)

    axislegend(ax, position=(0.0,0.0),bboxcolor=:transparent, labelsize=22)
    axislegend(ax2, position=(1.0,1.0),bboxcolor=:transparent, labelsize=22)

    #vlines!(ax, [1.0],linestyle=:dash,color=:black)
    #vlines!(ax2, [1.0],linestyle=:dash,color=:black)

    save("RG_planar.pdf",f)

    
end


############################################################
# Plotting graphs cylinder
############################################################

function plot_katchalsky(xi_set; a=3.0, L=1e5, lB=7.0,length_t=1000)


    f=Figure(fonts = (; regular = "CMU Serif"),size = (1000, 500))
    Label(f[1,1:2],L"FKL solution ($L=$ %$(L/10000) $\mu$m)",fontsize=28)
    ax=Axis(f[2,1], ylabel=L"\xi_\mathrm{R}(t)", xlabel=L"t")
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18

    ax2=Axis(f[2,2], ylabel=L"u(t)", xlabel=L"t")
    ax2.xlabelsize=22
    ax2.ylabelsize=22
    ax2.xticklabelsize=18
    ax2.yticklabelsize=18

    for xi in xi_set
    
        u, xi_eff = katchalsky_solution(xi; a=a, L=L, lB=lB)

        ts = range(log(a / lB), log(L / lB), length=length_t)

        lines!(ax,ts, xi_eff.(ts),label=L"$\xi=$%$xi")
        lines!(ax2,ts,u.(ts),label=L"$\xi=$%$xi")
    end
    #xlims!(ax,0.024,5.1)
    #ylims!(ax,0.024,5.1)

    axislegend(ax, position=(1.0,1.0),bboxcolor=:transparent, labelsize=22)
    axislegend(ax2, position=(1.0,0.0),bboxcolor=:transparent, labelsize=22)

    hlines!(ax, [1.0],linestyle=:dash,color=:black)

    save("katchalsky$(L/10000).pdf",f)

    
end



function plot_compare(xi; a=3.0, L=1e5, lB=7.0,length_t=1000)


    f=Figure(fonts = (; regular = "CMU Serif"),size = (500, 500))
    ax=Axis(f[1,1], ylabel=L"\xi_\mathrm{R}(t)", xlabel=L"t",title="FKL vs RG")
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18

    
    
    u, xi_eff = katchalsky_solution(xi; a=a, L=L, lB=lB)

    ts = range(log(a / lB), log(L / lB), length=length_t)

    lines!(ax,ts, xi_eff.(ts),label="FKL")

    trange, xi_RG = xi_R(xi;a=a,L=L,lB=lB,dt=1e-4)
    lines!(ax, trange, xi_RG, label="RG")

    axislegend(ax, position=(1.0,1.0),bboxcolor=:transparent, labelsize=22)
    
    hlines!(ax, [1.0],linestyle=:dash,color=:black)

    save("xi_R.pdf",f)

    
end

function plot_beta(Lset; a=3.0, lB=7.0, length_xi=1000)


    f=Figure(fonts = (; regular = "CMU Serif"),size = (1000, 500))
    #Label(f[1,1:2],L"$\beta$ function)",fontsize=28)
    L=Lset[1]
    tL= log(L / lB)
    ax=Axis(f[1,1], ylabel=L"\beta(\xi_R;t_L,t)", xlabel=L"\xi_R",title=L"(a) $L=$%$(L/10000) $\mu$m ($t_L=$%$(round(tL, digits=2)))")
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18
    ax.titlesize=25
    L=Lset[2]
    tL= log(L / lB)
    ax2=Axis(f[1,2], xlabel=L"\xi_R",title=L"(b) $L=$%$(L/10000) $\mu$m ($t_L=$%$(round(tL, digits=2)))")
    ax2.xlabelsize=22
    ax2.ylabelsize=22
    ax2.xticklabelsize=18
    ax2.yticklabelsize=18
    ax2.titlesize=25

    L=Lset[1]
    tset=[3,7,9]
    xi_set = range(0, 3.5, length=length_xi)
    for t in tset
        beta = beta_vs_xi(xi_set, t; L=L, lB=lB)

        

        lines!(ax,xi_set, beta,label=L"$t=$%$t")
   
    end

    L=Lset[2]
    tset=[3,17,20]
    
    for t in tset
        beta = beta_vs_xi(xi_set, t; L=L, lB=lB)

        

        lines!(ax2,xi_set, beta,label=L"$t=$%$t")
   
    end
    xlims!(ax,-0.1,3.3)
    xlims!(ax2,-0.1,3.3)
    ylims!(ax,-15,1)
    ylims!(ax2,-15,1)

    axislegend(ax, position=(0.0,0.0),bboxcolor=:transparent, labelsize=22)
    axislegend(ax2, position=(0.0,0.0),bboxcolor=:transparent, labelsize=22)

    vlines!(ax, [1.0],linestyle=:dash,color=:black)
    vlines!(ax2, [1.0],linestyle=:dash,color=:black)

    save("beta.pdf",f)

    
end


function plot_RG(xi_set; a=3.0, Lset=[1e5,1e10], lB=7.0)


    
    f=Figure(fonts = (; regular = "CMU Serif"),size = (1000, 500))
    #Label(f[1,1:2],L"$\beta$ function)",fontsize=28)
    L=Lset[1]
    tL= log(L / lB)
    ax=Axis(f[1,1], ylabel=L"\xi_R(t)", xlabel=L"t",title=L"(a) $L=$%$(L/10000) $\mu$m ($t_L=$%$(round(tL, digits=2)))")
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18
    ax.titlesize=25

    L=Lset[2]
    tL= log(L / lB)
    ax2=Axis(f[1,2], xlabel=L"t",title=L"(b) $L=$%$(L/10000) $\mu$m ($t_L=$%$(round(tL, digits=2)))")
    ax2.xlabelsize=22
    ax2.ylabelsize=22
    ax2.xticklabelsize=18
    ax2.yticklabelsize=18
    ax2.titlesize=25

    
    L=Lset[1]

    
    for xi in xi_set
        trange, xi_RG = xi_R(xi;a=a,L=L,lB=lB,dt=1e-4)
        lines!(ax, trange, xi_RG, label=L"$\xi=$%$xi")
    end
    
    axislegend(ax, position=(1.0,1.0),bboxcolor=:transparent, labelsize=22)
    
    hlines!(ax, [1.0],linestyle=:dash,color=:black)


    L=Lset[2]
    for xi in xi_set
        trange, xi_RG = xi_R(xi;a=a,L=L,lB=lB,dt=1e-4)
        lines!(ax2, trange, xi_RG, label=L"$\xi=$%$xi")
    end

    axislegend(ax2, position=(1.0,1.0),bboxcolor=:transparent, labelsize=22)
    
    hlines!(ax2, [1.0],linestyle=:dash,color=:black)

    save("RG.pdf",f)

    
end


#xi_set = [0.5, 0.7, 0.9, 1.1, 2.0,3.0,4.0]
#plot_katchalsky(xi_set,L=1e5)
#plot_katchalsky(xi_set,L=1e10)
#plot_compare(2.0; a=3.0, L=1e12, lB=7.0,length_t=1000)
#plot_RG(xi_set; a=3.0, lB=7.0)
#plot_beta([1e5,1e10])

plot_RG_planar([2,1427,1428.3],[0.5,1.5,3.0])