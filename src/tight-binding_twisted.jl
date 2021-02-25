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

struct get_hamiltonian
    intra
    tx
    ty
    txy
    geometry::geometry_twisted
end


#for twisted lattice, number of sites very big, need sparse matrix, to avoid reallocation of data, create it from [row,col,val], and preassign large enough length to row/col/val.
function get_H_twisted_bilayer(R,r,t_xy,t_z,d,lambda,mode,theta,r0,interlayer_bias=0)  # d: inter-layer distance, lambda: parameter, mode="uniform_hop" or "pi-flux"
    t0=now()
    N=size(R,1)

    #H=spzeros(Complex,N,N)                        # H_QSL is spin-degen
    I,J,value=ones(Int,N*100),ones(Int,N*100),zeros(Complex,N*100)
    n=1
    for i=1:N
        if R[i][3]!=0
           I[n]=i
           J[n]=i
           value[n]=interlayer_bias/2
           n+=1
        else
           I[n]=i
           J[n]=i
           value[n]=-interlayer_bias/2
           n+=1
        end

        for j=1:N
            dR=R[i]-(R[j]+r)
            z=dR[3]
            L=norm(dR)
                               
            if abs(z)>0.01
                if abs(t_z*(z^2/L^2)*exp(-lambda*(L-d)))>0.0001
                    I[n]=i
                    J[n]=j
                    value[n]=t_z*(z^2/L^2)*exp(-lambda*(L-d))
                    n+=1
                    #H[i,j]=H[i,j]+t_z*(z^2/L^2)*exp(-lambda*(L-d))         # from j to i
                end
            elseif 0.9<L<1.1
                t=get_hopping(R[i],R[j]+r,t_xy,mode,theta,r0) 
                I[n]=i
                J[n]=j
                value[n]=t
                n+=1                           
            end
        end
    end
    H=sparse(I,J,value,N,N)
    H=dropzeros(H)
    t1=now()
    println("time to construct H: ",t1-t0)
    return H
end

# for H_inter, H.tx+ty is always 0,only H.tx H.ty and H.tx-ty is nonzero, because theta(a1,a2)=60degree
function get_H_inter_twisted_bilayer(g::geometry_twisted, t_xy,t_z,lambda,mode,interlayer_bias=0,theta2=0)    #theta2: additional rotation of second layer w.r.t. first layer, can be n*2pi/3
    R=g.sites
    d=g.inter_distance
    theta=g.twist_angle
    theta+=theta2
    r0=g.inter_alignment
    H_intra=get_H_twisted_bilayer(R,[0,0,0],t_xy,t_z,d,lambda,mode,theta,r0,interlayer_bias)
    H_tx=get_H_twisted_bilayer(R,g.inter_vector.x,t_xy,t_z,d,lambda,mode,theta,r0)
    H_ty=get_H_twisted_bilayer(R,g.inter_vector.y,t_xy,t_z,d,lambda,mode,theta,r0)
    H_txy=get_H_twisted_bilayer(R,g.inter_vector.xy,t_xy,t_z,d,lambda,mode,theta,r0)
    H=get_hamiltonian(H_intra,H_tx,H_ty,H_txy,g)
    return H    
end

function get_Hk(H::get_hamiltonian,k::Array{Float64,1})
    H0=exp(im*dot(k,H.geometry.inter_vector.x))*H.tx+exp(im*dot(k,H.geometry.inter_vector.y))*H.ty+exp(im*dot(k,H.geometry.inter_vector.xy))*H.txy
    Hk=H.intra+H0+adjoint(H0)  
    return Hk 
end

function get_Hk_ribbon(H::get_hamiltonian,k::Array{Float64,1})   # here k represents the phase difference
    @assert k[1]*k[2]==0 "please input 1D momenta"
    k[1]==0 ? H0=exp(im*k[2])*H.ty : H0=exp(im*k[1])*H.tx
        
    Hk=H.intra+H0+adjoint(H0)  
    return Hk
end

function get_hopping(R,R2,t_xy,mode,theta,r0)
    factor=1
    if mode=="uniform_hop"
        if R[3]==1
            factor=1
        end
        return t_xy*factor
    elseif mode=="pi_flux"
        U=rotation_along_z_3D(-(theta))
        if R[3]==1
            R=U*(R-r0)
            R2=U*(R2-r0)
            factor=1
        end
        x=round(R[1]-R2[1],digits=2)
        y=R[2]-R2[2]
        if x*y<-0.01                             # with this method, it is better that the original lattice passes (0,0)
            hop=t_xy
        elseif x>0.1&&y>0.1
            if mod(round(-R[2]*2/sqrt(3)),2)==0  # from B to A, R==A
                hop=-t_xy
            else
                hop=t_xy
            end
        elseif x<-0.1&&y<-0.1
            if mod(round(-R[2]*2/sqrt(3)),2)==0
                hop=t_xy
            else
                hop=-t_xy
            end
        else
            if mod(round(-R[2]*2/sqrt(3)),2)==0
                hop=t_xy
            else
                hop=-t_xy
            end
        end
        return hop*factor                              
    else
        println("invalid mode, default to uniform_hop")
        return t_xy
    end
end

function plot_H(H,R)
    plt.figure(figsize=(6,6),dpi=80)
    N=size(H,1)
    for i=1:N
        for j=1:N
            if real(H[i,j])>0.1
                plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="red",linewidth=norm(H[i,j])*1)
            elseif real(H[i,j])<-0.1
                plt.plot((R[i][1],R[j][1]), (R[i][2],R[j][2]),color="blue",linewidth=norm(H[i,j])*1)
            end
        end
    end
    plt.axis([-5,5,-5,5])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end

function plot_H_2(H_inter,R,inter_vector)
    plt.figure(figsize=(6,6),dpi=80)
    n=size(H_inter,1)
    
    for i=1:n
        H=H_inter[i]
        r=inter_vector[i]
        R2=[x+r for x in R]

        N=size(H,1)
        for i=1:N
            for j=1:N
                if real(H[i,j])>0.1
                    plt.plot((R[i][1],R2[j][1]), (R[i][2],R2[j][2]),color="red",linewidth=norm(H[i,j])*1)
                elseif real(H[i,j])<-0.1
                    plt.plot((R[i][1],R2[j][1]), (R[i][2],R2[j][2]),color="blue",linewidth=norm(H[i,j])*1)
                end
            end
        end
    end
    plt.axis([-5,5,-5,5])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end



function test_tight_binding_twisted()
    n_a1=50
    n_a2=50
    m=1
    r=1
    r0=[0,0,0]
    d=1
    lattice,mode="triangular","pi_flux"
    BC="PBC"

    g=get_twisted_lattice_2(lattice,n_a1,n_a2,m,r,d,r0)
    println(g.twist_angle*180/pi)
    println(size(g.sites,1))
    
    t_xy=1
    t_z=0
    lambda=10
    H=get_H_inter_twisted_bilayer(g,t_xy,t_z,lambda,mode,0,pi)
    Hk=get_Hk(H,[0.0,0.0,0.0])
    plot_H(Hk,g.sites)
end

#test_tight_binding_twisted()
