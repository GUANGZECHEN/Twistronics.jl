#using Twistronics
include("../../src/Twistronics_local.jl")
np=pyimport("numpy")

n_a1=90
n_a2=90
m=17          # in version 0.1.0, when m is even, it does not work: the way the I construct the lattice destroys bipartite nature of honeycomb (will fix in later versions)
r=1
r0=[0,0,0]   # inter-layer alignment
d=1          # inter-layer distance
lattice,mode="honeycomb","uniform_hop"
# m and r determine twist angle

# build the geometry, it contains all sites in the supercell 
g=get_twisted_lattice_2(lattice,n_a1,n_a2,m,r,d,r0,mode)
println(g.inter_vector)
#plot_R(g.sites) # plot the sites 

t_xy=-1 # intra-layer hopping
t_z=0.3 # inter-layer hopping
bias=0.1 # interlayer bias
lambda=10 # this parameter determines the decay of inter-layer hopping w.r.t distance
    
H=get_H_inter_twisted_bilayer(g,t_xy,t_z,lambda,mode,bias,0)  # get the Hamiltonian of twisted bilayer
P=get_valley_operator_inter_twisted_bilayer(g,mode)
#Twistronics.plot_H(H.intra,g.sites)              # Plot the intra-unitcell Hamiltonian to see if hoppings are correct

b1,b2=get_reciprocal_vector(g.inter_vector.x,g.inter_vector.y) # get reciprocal vector for the supercell
k_points=[i/300*(b1+2*b2) for i=1:300]  # plot bandstrcture along Gamma-K-K'-Gamma
#k_points=[i/100*(b1-b2) for i=1:100]
n_bands=26                            # compute 40 eigenvalues for each k point, using Arpack
momenta, energies, valley_polarization=get_band_twisted(H,P,k_points,n_bands) # get bandstructure with valley-polarization
np.savetxt("bands.OUT",np.transpose([momenta, energies, valley_polarization]))

# plot bandstructure
plt.figure(figsize=(6,6),dpi=80)
colors=[(1,0,0),(0,1,0),(0,0,1)]
cm=colormap.LinearSegmentedColormap.from_list("my_list", colors) 
sc=plt.scatter(momenta, energies,c=valley_polarization,cmap=cm,vmin=-1,vmax=1)  # color indicates valley polarization
plt.colorbar(sc)
plt.axis([-1,300+1,-0.1,0.1])
plt.xlabel("k")
plt.ylabel("E")
plt.show()

# Note: for TBG, all low energy states are valley-degenerate
