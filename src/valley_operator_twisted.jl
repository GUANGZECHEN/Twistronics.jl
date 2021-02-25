#using LinearAlgebra
#using DelimitedFiles
#using Statistics
#using Plots
#using DelimitedFiles
#using Random
#using PyCall
#using SparseArrays

#plt = pyimport("matplotlib.pyplot")

#include("Geometry_twisted.jl")


### for triangular lattice, use pi-flux hopping!!!

function get_valley_operator_twisted_bilayer(R,r,mode,theta,r0,layer_factor=1)  # d: inter-layer distance, lambda: parameter, mode="uniform_hop" or "pi-flux"
    t0=now()
    N=size(R,1)

    I,J,value=ones(Int,N*13),ones(Int,N*13),zeros(Complex,N*13)
    n=1
    for i=1:N
        for j=1:N
            dR=R[i]-(R[j]+r)
            z=dR[3]
            L=norm(dR)
                               
            if abs(z)<0.01 && L<2
                t=get_valley_operator_hopping(R[i],R[j]+r,mode,theta,r0,layer_factor) 
                if abs(t)>0
                    I[n]=i
                    J[n]=j
                    value[n]=t
                    n+=1
                end                           
            end
        end
    end
    P=sparse(I,J,value,N,N)
    P=dropzeros(P)
    t1=now()
    println("time to construct P: ",t1-t0)
    return P
end

function get_valley_operator_hopping(R,R2,mode,theta,r0,layer_factor=1)
    factor=1
    U=rotation_along_z_3D(-(theta))
    if R[3]==1
        R=U*(R-r0)
        R2=U*(R2-r0)
        factor=layer_factor
    end
    x=R[1]-R2[1]
    y=R[2]-R2[2]
    if mode=="uniform_hop"    #valley operator for honeycomb
        hop=0
        if round(x/sqrt(3),digits=3)==1 && round(y,digits=3)==0
            hop=im
        elseif round(x/sqrt(3),digits=3)==-1 && round(y,digits=3)==0
            hop=-im
        elseif round(x*2/sqrt(3),digits=3)==-1 && round(y,digits=3)==3/2
            hop=im
        elseif round(x*2/sqrt(3),digits=3)==1 && round(y,digits=3)==-3/2
            hop=-im
        elseif round(x*2/sqrt(3),digits=3)==-1 && round(y,digits=3)==-3/2
            hop=im
        elseif round(x*2/sqrt(3),digits=3)==1 && round(y,digits=3)==3/2
            hop=-im
        end
        return hop*layer_factor/(3*sqrt(3))
    elseif mode=="pi_flux"
        hop=0
        if round(x,digits=3)==1 && round(y,digits=3)==0
            hop=-im
        elseif round(x,digits=3)==-1 && round(y,digits=3)==0
            hop=im
        elseif round(x,digits=3)==0 && round(y/sqrt(3),digits=3)==-1
            hop=-im
        elseif round(x,digits=3)==0 && round(y/sqrt(3),digits=3)==1
            hop=im
        end
        return hop*factor/4                              
    else
        println("invalid mode for valley operator")
    end
end

function get_valley_operator_inter_twisted_bilayer(g::geometry_twisted,mode,layer_factor=1)
    R=g.sites
    theta=g.twist_angle
    r0=g.inter_alignment

    P_intra=get_valley_operator_twisted_bilayer(R,[0,0,0],mode,theta,r0,layer_factor)
    P_tx=get_valley_operator_twisted_bilayer(R,g.inter_vector.x,mode,theta,r0,layer_factor)
    P_ty=get_valley_operator_twisted_bilayer(R,g.inter_vector.y,mode,theta,r0,layer_factor)
    P_txy=get_valley_operator_twisted_bilayer(R,g.inter_vector.xy,mode,theta,r0,layer_factor)
    P=get_hamiltonian(P_intra,P_tx,P_ty,P_txy,g)
    return P    
end

function get_Pk(P::get_hamiltonian,k::Array{Float64,1})
    P0=exp(im*dot(k,P.geometry.inter_vector.x))*P.tx+exp(im*dot(k,P.geometry.inter_vector.y))*P.ty+exp(im*dot(k,P.geometry.inter_vector.xy))*P.txy
    Pk=P.intra+P0+adjoint(P0)  
    return Pk 
end

function get_Pk_ribbon(P::get_hamiltonian,k::Array{Float64,1})   # here k represents the phase difference
    @assert k[1]*k[2]==0 "please input 1D momenta"
    k[1]==0 ? P0=exp(im*k[2])*P.ty : P0=exp(im*k[1])*P.tx
        
    Pk=P.intra+P0+adjoint(P0)  
    return Pk
end

function plot_P(P,R)
    plt.figure(figsize=(6,6),dpi=80)
    N=size(R,1)
    for i=1:N
        for j=i:N
            if imag(P[i,j])>0.1
                plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="red",linewidth=norm(P[i,j])*1)
            elseif imag(P[i,j])<-0.1
                plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="blue",linewidth=norm(P[i,j])*1)
            end
        end
    end
    plt.axis([-5,5,-5,5])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end


function get_layer_operator_twisted_bilayer(R)  # d: inter-layer distance, lambda: parameter, mode="uniform_hop" or "pi-flux"
    N=size(R,1)

    L=spzeros(Complex,N,N)   
    for i=1:N
        if R[i][3]==0
            L[i,i]=1
        else
            L[i,i]=-1
        end
    end
    return L
end
