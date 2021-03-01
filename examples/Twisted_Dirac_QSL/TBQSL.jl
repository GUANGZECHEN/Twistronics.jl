using Twistronics

n_a1=200
n_a2=200
m=9       
r=1
r0=[0,0,0]   # inter-layer alignment
d=1          # inter-layer distance
lattice,mode="triangular","pi_flux"
# m and r determine twist angle

# build the geometry, it contains all sites in the supercell 
g=Twistronics.get_twisted_lattice_2(lattice,n_a1,n_a2,m,r,d,r0,mode)

#Twistronics.plot_R(g.sites) # plot the sites 

t_xy=-1 # intra-layer hopping
t_z=0.36 # inter-layer hopping
lambda=10 # this parameter determines the decay of inter-layer hopping w.r.t distance
    
H=Twistronics.get_H_inter_twisted_bilayer(g,t_xy,t_z,lambda,mode)  # get the Hamiltonian of twisted bilayer
P=Twistronics.get_valley_operator_inter_twisted_bilayer(g,mode)    # get the valley operator of the twisted bilayer

#Twistronics.plot_H(H.intra,g.sites)              # Plot the intra-unitcell Hamiltonian to see if hoppings are correct

b1,b2=Twistronics.get_reciprocal_vector(g.inter_vector.x,g.inter_vector.y) # get reciprocal vector for the supercell
n_k=100
k_points_1=[i/n_k*b1/2 for i=0:n_k-1]
k_points_2=[b1/2+i/n_k*b2/2 for i=0:n_k-1]
k_points_3=[(b2+b1)/2-i/n_k*(b2+b1)/2 for i=0:n_k]
k_points=[k_points_1;k_points_2;k_points_3]          # plot bands along Gamma-M-M'-Gamma

n_bands=40                              # compute 40 eigenvalues for each k point, using Arpack
momenta, energies, valley_polarization=Twistronics.get_band_twisted(H,P,k_points,n_bands) # get bandstructure with valley-polarization

# plot bandstructure
Twistronics.plt.figure(figsize=(6,6),dpi=80)
colors=[(1,0,0),(0,1,0),(0,0,1)]
cm=Twistronics.colormap.LinearSegmentedColormap.from_list("my_list", colors) 
sc=Twistronics.plt.scatter(momenta, energies,c=valley_polarization,cmap=cm,vmin=-1,vmax=1)  # color indicates valley polarization
Twistronics.plt.colorbar(sc)
Twistronics.plt.axis([-1,300+1,-0.3,0.3])
Twistronics.plt.xlabel("k")
Twistronics.plt.ylabel("E")
Twistronics.plt.show()


