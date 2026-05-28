using CairoMakie #####
using MySave:@savevar,@loadvar
using SpecialFunctions



Base.@kwdef struct Param{Tp1,Tp2,Tp3,Tp4,Tp5,Tp6,Tp7}
	ds::Tp1=1e-4 # fictitious time
	n0::Tp2=1e-10 # normalized by lB^3
	xi::Tp3=2.0 # Manning parameter
	dt::Tp4=0.02 # mesh size for t=log(r/lB)
	Nt::Tp5=520 # number of mesh points for t
	a::Tp6=3/7 # rod radius in unit of lB
	err_c::Tp7=10^-8 # error for convergence
end





function bc!(A,xi,dt)
    A[1] = A[2]-2*xi*dt
    A[end] = A[end-1]
    return
end



"""
calc_boost!(u,boost,kappa,dt)

"""
function calc_boost!(u,boost,kappa,t,dt)
	for i in 2:length(u)-1
		boost[i]=(u[i+1]+u[i-1]-2*u[i])*(1/dt)^2+0.5*kappa^2*(exp(-u[i]+2*t[i])-exp(u[i]+2*t[i]))
	end
	boost[1]=0
	boost[end]=0
	err=sum(abs.(boost))/length(boost)
	return err
end

"""
update!(u,u1,boost,ds)

"""

function update!(u,u1,boost,ds)
    
	@. u1 = u + ds*boost
	@. u = u1
	return
end







function main(pa::Param)
	
	(; ds,n0,xi,dt,Nt,a,err_c) = pa

	u=zeros(Float64,Nt)
	u1=zeros(Float64,Nt)
	boost=zeros(Float64,Nt)

	kappa=sqrt(8*pi*n0)
	t=range(log(a),step=dt,length=Nt)

	err=1.0
	bc!(u,xi,dt)
	while err>err_c
		
		err=calc_boost!(u,boost,kappa,t,dt)
		update!(u,u1,boost,ds)
		bc!(u,xi,dt)
	end
	
	xi_R=similar(u)
	xi_R[1]=xi
	xi_R[end]=0
	for i in 2:Nt-1
		xi_R[i]=0.25*(u[i+1]-u[i-1])/dt
	end
	
	tkappa=-log(kappa)

	
	idx_kappa=findmin(x -> abs(x - tkappa), t)[2]

	@savevar u xi_R t 
	
	println("Done.")
	@show xi,tkappa,xi_R[idx_kappa] t[idx_kappa]
	return xi, xi_R[idx_kappa], tkappa
end


function plot_xi_eff(xi_set,n0_set)

	f=Figure(fonts = (; regular = "CMU Serif"),size = (700, 500))
    #Label(f[1,1:2],L"(a) RG analysis for planar geometry ($\ell_L=$%$(round(lL, digits=2)))",fontsize=28)
    
    
    ax=Axis(f[1,1], ylabel=L"\xi_\mathrm{eff}", xlabel=L"\xi",title=L"Numerical result for $\xi_\mathrm{eff}=\xi_R(t_\mathrm{Debye})$")
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18
    ax.titlesize=25

	
	xieff=zeros(Float64,length(xi_set))
	
	tkappa=0
	for n0 in n0_set
		i=1
		for xi in xi_set
			pa=Param(xi=xi,n0=n0)
			x,xieff[i],tkappa=main(pa)
			i+=1
		end
		scatter!(ax,xi_set,xieff,label=L"$t_\mathrm{Debye}=$%$(round(tkappa, digits=2))",markersize=25)
	end
	vlines!(ax, [1.0],linestyle=:dash,color=:black)

	axislegend(ax, position=(1.0,0.0),bboxcolor=:transparent, labelsize=22)

	save("xi_eff.pdf",f)
end


function plot_profile(xi_set,n0_1,n0_2)

	f=Figure(fonts = (; regular = "CMU Serif"),size = (1000, 500))
    #Label(f[1,1:2],L"(a) RG analysis for planar geometry ($\ell_L=$%$(round(lL, digits=2)))",fontsize=28)
    
	kappa1=sqrt(8*pi*n0_1)
	tD1=-log(kappa1)

	kappa2=sqrt(8*pi*n0_2)
	tD2=-log(kappa2)
    
    ax=Axis(f[1,1], ylabel=L"\xi_R(t)", xlabel=L"t",title=L"(a) $t_\mathrm{Debye}=$%$(round(tD1, digits=2))")
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18
    ax.titlesize=25

	ax2=Axis(f[1,2], ylabel=L"\xi_R(t)", xlabel=L"t",title=L"(b) $t_\mathrm{Debye}=$%$(round(tD2, digits=2))")
    ax2.xlabelsize=22
    ax2.ylabelsize=22
    ax2.xticklabelsize=18
    ax2.yticklabelsize=18
    ax2.titlesize=25

	
	
	
	tkappa=0
	xieff=0
	for xi in xi_set
		pa=Param(xi=xi,n0=n0_1)
		x,xieff,tkappa=main(pa)
		t,xi_R=@loadvar t xi_R
		lines!(ax,t,xi_R,label=L"$\xi=$%$xi")
	end
	vlines!(ax, [tD1],linestyle=:dash,color=:black)

	for xi in xi_set
		pa=Param(xi=xi,n0=n0_2)
		x,xieff,tkappa=main(pa)
		t,xi_R=@loadvar t xi_R
		lines!(ax2,t,xi_R,label=L"$\xi=$%$xi")
	end
	vlines!(ax2, [tD2],linestyle=:dash,color=:black)

	axislegend(ax, position=(0.3,0.9),bboxcolor=:transparent, labelsize=22)
	axislegend(ax2, position=(1.0,1.0),bboxcolor=:transparent, labelsize=22)

	save("profile_salt.pdf",f)
end

function plot_profile_LPB(xi,n0_1,n0_2)

	f=Figure(fonts = (; regular = "CMU Serif"),size = (700, 500))
    #Label(f[1,1:2],L"(a) RG analysis for planar geometry ($\ell_L=$%$(round(lL, digits=2)))",fontsize=28)
    
	kappa1=sqrt(8*pi*n0_1)
	tD1=-log(kappa1)

	kappa2=sqrt(8*pi*n0_2)
	tD2=-log(kappa2)
    
    ax=Axis(f[1,1], ylabel=L"\xi_R(t)", xlabel=L"t",title=L"Profiles for large and small salt conc. ($\xi=$%$xi)")
    ax.xlabelsize=22
    ax.ylabelsize=22
    ax.xticklabelsize=18
    ax.yticklabelsize=18
    ax.titlesize=25

	
	tkappa=0
	xieff=0
	
	pa=Param(xi=xi,n0=n0_1)
	x,xieff,tkappa=main(pa)
	t,xi_R=@loadvar t xi_R
	lines!(ax,t,xi_R,label=L"$t_\mathrm{Debye}=$%$(round(tD1, digits=2))")
	scatter!(ax,tD1,xieff)
	tl=tD1:0.02:t[end]
	var=exp.(tl .- tD1)
	xiR_lin=@. (xieff/besselk(1,1.0))*var*besselk(1,var)
	lines!(ax,tl,xiR_lin,linestyle=:dash,color=:black,linewidth=2)

	pa=Param(xi=xi,n0=n0_2)
	x,xieff,tkappa=main(pa)
	t,xi_R=@loadvar t xi_R
	lines!(ax,t,xi_R,label=L"$t_\mathrm{Debye}=$%$(round(tD2, digits=2))")
	scatter!(ax,tD2,xieff)
	tl=tD2:0.02:t[end]
	var=exp.(tl .- tD2)
	xiR_lin=@. (xieff/besselk(1,1.0))*var*besselk.(1,var)
	lines!(ax,tl,xiR_lin,linestyle=:dash,color=:black,label=L"LPB $\xi_R^\mathrm{lin}(t)$",linewidth=2)
	

	

	axislegend(ax, position=(0.9,0.9),bboxcolor=:transparent, labelsize=22)
	

	save("profile_salt_LPB.pdf",f)
end


###### xi vs xi_eff ########
#plot_xi_eff([0.5,1,2,3,4,5,6],[1e-8,1e-6,1e-4,1e-3])

##### xi_R(t) ########
#plot_profile([0.5,1,2,3,5],1e-8,1e-2)

##### xi_R(t) and the solution of linearized PB ########
plot_profile_LPB(2.0, 1e-8,1e-3)