#using LinearAlgebra
#using DelimitedFiles
#using Statistics
#using Plots
#using DelimitedFiles
#using Random
#using SparseArrays
#using Arpack

#include("Geometry_twisted.jl")
#include("tight-binding_twisted.jl")
#include("valley_operator_twisted.jl")
#using PyCall
#plt = pyimport("matplotlib.pyplot")


#using Dates

function get_k_points(b1,b2,n_k)
    k_points=[b2*i/n_k+(b1-b2/2)/2 for i=0:n_k]
    return k_points
end

function get_band_twisted(H::get_hamiltonian,P::get_hamiltonian,k_points::Array{Array{Float64,1},1},n_bands,dim=2) # dim=1 for ribbon geometry
    n_k=size(k_points,1)
    momenta=Array{Float64}(undef,Int(n_k*n_bands))
    energies=Array{Float64}(undef,Int(n_k*n_bands))
    valley_polarization=Array{Float64}(undef,Int(n_k*n_bands))
    index=1
    for i=1:n_k
        k=k_points[i]
        if dim==2
            Hk=get_Hk(H,k)
            Pk=get_Pk(P,k)#get_layer_operator_twisted_bilayer(R_unitcell)#
        else
            Hk=get_Hk_ribbon(H,k)
            Pk=get_Pk_ribbon(P,k)#get_layer_operator_twisted_bilayer(R_unitcell)#           
        end
        
        eigvals,eigvecs=eigs(Hk,nev=n_bands,sigma=0.0001,which=:LM,maxiter=10000)
        #F=eigen(Matrix(H))
        #eigvals,eigvecs=F.values,F.vectors
        #eigvals=real(eigvals)
        #n_bands=size(eigvals,1)
        for i_band=1:n_bands
            momenta[index]=i
            energies[index]=real(eigvals[i_band])
            psi=eigvecs[:,i_band]
            valley_polarization[index]=real(adjoint(psi)*Pk*psi)
            index+=1
        end
    end
    return momenta, energies, valley_polarization
end

function get_band_twisted_2(H::get_hamiltonian,P::get_hamiltonian,k_points::Array{Array{Float64,1},1},n_bands,dim=2) # dim=1 for ribbon geometry
    n_k=size(k_points,1)
    momenta=Array{Float64}(undef,Int(n_k*n_bands))
    energies=Array{Float64}(undef,Int(n_k*n_bands))
    valley_polarization=Array{Float64}(undef,Int(n_k*n_bands))
    index=1
    for i=1:n_k
        k=k_points[i]
        if dim==2
            Hk=get_Hk(H,k)
            Pk=get_Pk(P,k)#get_layer_operator_twisted_bilayer(R_unitcell)#
        else
            Hk=get_Hk_ribbon(H,k)
            Pk=get_Pk_ribbon(P,k)#get_layer_operator_twisted_bilayer(R_unitcell)#           
        end
        Hk=Hk+0.05*Pk
        eigvals,eigvecs=eigs(Hk,nev=n_bands,sigma=0.0001,which=:LM,maxiter=10000)
        #F=eigen(Matrix(H))
        #eigvals,eigvecs=F.values,F.vectors
        #eigvals=real(eigvals)
        #n_bands=size(eigvals,1)
        for i_band=1:n_bands
            momenta[index]=i
            energies[index]=real(eigvals[i_band])
            psi=eigvecs[:,i_band]
            valley_polarization[index]=real(adjoint(psi)*Pk*psi)
            index+=1
        end
    end
    return momenta, energies, valley_polarization
end

function get_eigen_twisted(H::get_hamiltonian,k_points::Array{Array{Float64,1},1},n_bands)
    n_k=size(k_points,1)
    energies=Array{Complex}(undef,n_k*n_bands)
    wavefuncs=Array{Any}(undef,n_k*n_bands)
    t0=now()
    ii=1
    for i=1:n_k
        k=k_points[i]            
        
        Hk=get_Hk(H,k)
  
        eigvals,eigvecs=eigs(Hk,nev=n_bands,sigma=0.0001,which=:LM,maxiter=10000)
            
        for i_band=1:n_bands
            energies[ii]=eigvals[i_band]
            wavefuncs[ii]=eigvecs[:,i_band]
            ii+=1
        end
    end

    t1=now()
    open("log.txt","a") do io                     
        println(io,"summed over", n_k," k points, time spent: ",t1-t0)
    end

    return energies, wavefuncs
end

function get_DOS(energies,omega,eta)
    DOS=0
    for E in energies
        DOS=DOS+1/(E-omega+im*eta)
    end
    DOS=-1/pi*imag(DOS)
    return real(DOS)
end

function get_LDOS(energies,wavefuncs,omega,eta)
    n_sites=size(wavefuncs[1],1)
    n_energies=size(energies,1)

    LDOS=zeros(Complex{Float64},n_sites)
    for i=1:n_energies
        psi=wavefuncs[i]
        E=energies[i]
        for j=1:n_sites
            LDOS[j]=LDOS[j]+abs2(psi[j])/(omega-E+im*eta)
        end
    end
    LDOS=real(-1/pi*imag(LDOS))
    return LDOS
end

function plot_LDOS(R_unitcell,LDOS)
    plt.figure(figsize=(6,6),dpi=80)
    N=size(R_unitcell,1)
    x=Array{Float64}(undef,N)
    y=Array{Float64}(undef,N)
    for i=1:N
        x[i]=R_unitcell[i][1]
        y[i]=R_unitcell[i][2]
    end
    #println(size(x,1),size(y,1))
    colors=[(1,0,0),(0,1,0),(0,0,1)]
    cm=colormap.LinearSegmentedColormap.from_list("my_list", colors)
    sc=plt.scatter(x, y,c=LDOS,cmap=cm)
    plt.colorbar(sc)
    plt.axis([-40,40,-40,40])
    plt.xlabel("")
    plt.ylabel("")
    #plt.savefig()
    plt.show()
end

