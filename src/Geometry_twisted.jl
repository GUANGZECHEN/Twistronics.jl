struct unit_vector
    x
    y
    xy
end

struct geometry_twisted
    lattice::String
    sites
    inter_vector::unit_vector
    dimension::Int
    twist_angle::Float64
    inter_distance
    inter_alignment
end

function get_theta(m,r,a1,a2)
    theta=acos((3*m^2+3*m*r+r^2/2)/(3*m^2+3*m*r+r^2))
    A1=m*a2+(m+r)*a1
    A2=-(m+r)*a2+(2*m+r)*a1                                     # the a1,a2 in literature is the opposite definition from here
    n=3*m^2+3*m*r+r^2
    return theta, A1, A2, n
end

function get_parameter_for_theta(theta)
    m=1
    r=1
    theta2=acos((3*m^2+3*m*r+r^2/2)/(3*m^2+3*m*r+r^2))

    while theta2>theta
        m+=1
        theta2=acos((3*m^2+3*m*r+r^2/2)/(3*m^2+3*m*r+r^2))
    end
    m-=1
    theta2=acos((3*m^2+3*m*r+r^2/2)/(3*m^2+3*m*r+r^2))
    factor=theta2/theta

    return m, factor
end

function rotation_along_z_3D(theta)   
    U=[[cos(theta),sin(theta),0] [-sin(theta),cos(theta),0] [0,0,1]]     #in this way U=anti-clock rotation of theta
    return U
end

function geometry_simple_3D_rotated(n,m,a1,a2,z,theta,r0)
    U=rotation_along_z_3D(theta)
    a1=U*a1
    a2=U*a2
    R0=[0,0,z]-(n-1)/2*a1-(m-1)/2*a2+r0                                  # s.t. the second layer is rotated w.r.t. [0,0]
    R=Array{Any}(undef,n*m)
    ii=1
    for i=0:n-1
        for j=0:m-1
            R[ii]=i*a1+j*a2+R0
            ii=ii+1
        end
    end
    return R,a1,a2
end

function geometry_bipartite_3D_rotated(n,m,a1,a2,a3,z,theta,r0)
    U=rotation_along_z_3D(theta)
    a1=U*a1
    a2=U*a2
    a3=U*a3
    R0=[0,0,z]-(n-1)/2*a1-(m-1)/2*a2+r0
    R=Array{Any}(undef,2*n*m)
    ii=1
    for i=0:n-1
        for j=0:m-1
            for k=0:1
                R[ii]=i*a1+j*a2+k*a3+R0
                ii=ii+1
            end
        end
    end
    return R,a1,a2
end

function nearest(site_a,site_b)
    return 0.1<norm(site_a-site_b)<1.1
end

function next_nearest(site_a,site_b)
    return 1.1<norm(site_a-site_b)<1.9
end

function next_next_nearest(site_a,site_b)
    return 1.9<norm(site_a-site_b)<2.1
end

function plot_R(R)
    plt.figure(figsize=(6,4),dpi=100)
    N=size(R,1)
    x=Array{Float64}(undef,Int(N))
    y=Array{Float64}(undef,Int(N))
    x2=Array{Float64}(undef,Int(N))
    y2=Array{Float64}(undef,Int(N))
    j1=1
    j2=1
    for i=1:N
        if R[i][3]==0
            x[j1]=R[i][1]
            y[j1]=R[i][2]
            j1+=1
        else
            x2[j2]=R[i][1]
            y2[j2]=R[i][2]
            j2+=1
        end
    end
    if j1>1
        x=x[1:j1-1]
        y=y[1:j1-1]
        plt.scatter(x,y,color="red",s=10)
    end
    if j2>1
        x2=x2[1:j2-1]
        y2=y2[1:j2-1]
        plt.scatter(x2,y2,color="blue",s=10)
    end
    plt.axis([-50,50,-30,30])
    plt.xlabel("")
    plt.ylabel("")
    plt.show()
end

function add_vacancy(R,vacant_site)
    deleteat!(R,vacant_site)
    return R
end

function merge_R(R,R2)
    M=size(R,1)
    N=size(R2,1)
    RR=Array{Any}(undef,M+N)
    for i=1:M
        RR[i]=R[i]
    end
    for j=1:N
        RR[N+j]=R2[j]
    end
    return R
end

function check_site_in_unit_cell(site,A1,A2)      #check whether site is in unitcell [-1/2,1/2]*A1+[-1/2,1/2]*A2, site can be 3D
    num=A1[1]*A2[2]-A1[2]*A2[1]
    num2=site[1]*A1[2]-site[2]*A1[1]
    num3=site[1]*A2[2]-site[2]*A2[1]              
    alpha=num3/num
    beta=-num2/num                                # site=alpha*a1+beta*a2
    if -0.50<round(alpha,digits=10)<=0.50 && -0.50<round(beta,digits=10)<=0.50
        return true
    else 
        return false
    end
end


function check_R(R_unitcell,inter_vector,R)
    n_inter=size(inter_vector,1)
    R=[round.(x,digits=3) for x in R]
    #println(R)
    N=size(R_unitcell,1)
    for ii=2:n_inter
        r=inter_vector[ii]
        for i=1:N
            R0=R_unitcell[i]+r
            R0=round.(R0,digits=3)
            #println(R0)
            if !(R0 in R)
                println(R0)
                #println(R)
                return false
            end
        end
    end

    return true
end

function plot_R_unitcell(R_unitcell,inter_vector)
    n_inter=size(inter_vector,1)
    N=size(R_unitcell,1)
    for ii=2:n_inter
        r=inter_vector[ii]
        for i=1:N
            push!(R_unitcell,R_unitcell[i]+r)
        end
    end
    plot_R(R_unitcell)
end

function get_reciprocal_vector(A1,A2) 
    A3=[0,0,1]
    B1=2*pi*cross(A2,A3)/dot(A1,cross(A2,A3))
    B2=2*pi*cross(A1,A3)/dot(A2,cross(A1,A3))
    #b1=[B1[1],B1[2]]
    #b2=[B2[1],B2[2]]
    return B1,B2
end

function get_inter_vector_3D(n,m,a1,a2)
    inter_vector=[[0,0,0],n*a1,m*a2,m*a2-n*a1]
    return inter_vector
end


function get_lattice_3D_rotated(lattice,n_a1,n_a2,z,theta,r0,a)   #theta in units of 1, z: layer-index  #add confinement in unitcell, r0: parallel alignment of layer, a=lattice constant   
    if lattice=="triangular"
        a1=[1,0,0]*a
        a2=[1/2,-sqrt(3)/2,0]*a
        #println(a1,a2)
        R,a1,a2=geometry_simple_3D_rotated(n_a1,n_a2,a1,a2,z,theta,r0)
    elseif lattice=="square"
        a1=[1,0,0]*a
        a2=[0,1,0]*a
        R,a1,a2=geometry_simple_3D_rotated(n_a1,n_a2,a1,a2,z,theta,r0)
    elseif lattice=="honeycomb"
        a1=[sqrt(3),0,0]*a
        a2=[sqrt(3)/2,-3/2,0]*a
        a3=[sqrt(3)/2,-1/2,0]*a
        R,a1,a2=geometry_bipartite_3D_rotated(n_a1,n_a2,a1,a2,a3,z,theta,r0)
    else
        println("invalid lattice")
    end
    
    return R,a1,a2
end

function get_twisted_lattice_2(lattice,n_a1,n_a2,m,r=1,d=1,r0=[0,0,0],mode="pi_flux")
    mode=="pi_flux" ? a=2 : a=1
    R,a1,a2=get_lattice_3D_rotated(lattice,n_a1,n_a2,0,0,[0,0,0],a)
    theta, A1, A2, n_sites=get_theta(m,r,a1,a2)
    R2,a1_2,a2_2=get_lattice_3D_rotated(lattice,n_a1,n_a2,d,theta,r0,a)
    R_unitcell,n1,n2=(mode=="pi_flux" ? get_unit_cell_2(R,R2,A1,A2,a1,a2,a1_2,a2_2,n_sites) : get_unit_cell(R,R2,A1,A2,n_sites))
    N=size(R_unitcell,1)

    @assert n1==n2 "number of sites on different layers mismatch",n1,n2
    @assert n1+n2==N "primitive cell too small"
    inter_vector=unit_vector(A1,A2,A1-A2)
    g=geometry_twisted(lattice, R_unitcell, inter_vector, 2, theta, d, r0)

    return g
end

function get_twisted_lattice_ribbon(lattice,n_a1,n_a2,m,r=1,d=1,r0=[0,0,0],mode="pi_flux",n_A1=1,n_A2=1)
    mode=="pi_flux" ? a=2 : a=1
    R,a1,a2=get_lattice_3D_rotated(lattice,n_a1,n_a2,0,0,[0,0,0],a)
    theta, A1, A2, n_sites=get_theta(m,r,a1,a2)
    R2,a1_2,a2_2=get_lattice_3D_rotated(lattice,n_a1,n_a2,d,theta,r0,a)
    R_unitcell,n1,n2=(mode=="pi_flux" ? get_unit_cell_2(R,R2,A1,A2,a1,a2,a1_2,a2_2,n_sites,n_A1,n_A2) : get_unit_cell(R,R2,A1,A2,n_sites,n_A1,n_A2))
    N=size(R_unitcell,1)

    @assert n1==n2 "number of sites on different layers mismatch",n1,n2
    @assert n1+n2==N "primitive cell too small"
    inter_vector=unit_vector(n_A1*A1,n_A2*A2,n_A1*A1-n_A2*A2)
    g=geometry_twisted(lattice, R_unitcell, inter_vector, 2, theta, d, r0)

    return g
end

function get_unit_cell_2(R,R2,A1,A2,a1,a2,a1_2,a2_2,n_sites,n_A1=1,n_A2=1)
    A1=n_A1*A1
    A2=n_A2*A2
    n_sites=n_sites*8*n_A1*n_A2   # for triangular, the factor is 2*4=8
    R_unitcell=Array{Any}(undef,n_sites)
    N=size(R,1)
    n1=0
    n2=0
    n=1
    for i=1:N
        if check_site_in_unit_cell(R[i],A1,A2)
            n1=n1+4
            R_unitcell[n]=R[i]
            R_unitcell[n+1]=R[i]+a1/2
            R_unitcell[n+2]=R[i]+a2/2
            R_unitcell[n+3]=R[i]+a1/2+a2/2
            n=n+4
        end
        if check_site_in_unit_cell(R2[i],A1,A2)
            n2=n2+4
            R_unitcell[n]=R2[i]
            R_unitcell[n+1]=R2[i]+a1_2/2
            R_unitcell[n+2]=R2[i]+a2_2/2
            R_unitcell[n+3]=R2[i]+a1_2/2+a2_2/2
            n=n+4
        end
    end

    return R_unitcell,n1,n2
end

function get_unit_cell(R,R2,A1,A2,n_sites,n_A1=1,n_A2=1)
    A1=n_A1*A1
    A2=n_A2*A2
    n_sites=n_sites*4*n_A1*n_A2 # 4 for honeycomb (bipartite)
    R_unitcell=Array{Any}(undef,n_sites)
    N=size(R,1)
    n1=0
    n2=0
    n=1
    for i=1:N
        if check_site_in_unit_cell(R[i],A1,A2)
            n1=n1+1
            R_unitcell[n]=R[i]
            n=n+1
        end
        if check_site_in_unit_cell(R2[i],A1,A2)
            n2=n2+1
            R_unitcell[n]=R2[i]
            n=n+1
        end
    end

    return R_unitcell,n1,n2
end

function remove_layer_R(R,layer_index)
    N=size(R,1)
    for i=1:N
        j=N+1-i
        if R[j][3]==layer_index
            deleteat!(R,j)
        end
    end
    return R
end

function remove_layer_g(g,layer_index)
    R=remove_layer_R(g.sites,layer_index)
    return geometry_twisted(g.lattice,R,g.inter_vector,g.dimension,g.twist_angle,g.inter_distance,g.inter_alignment)
end

function test_geometry_twisted()
    n_a1=50
    n_a2=50
    m=1
    r=1
    r0=[0,0,0]
    d=1
    lattice,mode="triangular","pi_flux"
    #lattice,mode="honeycomb","uniform_hop"

    g=get_twisted_lattice_2(lattice,n_a1,n_a2,m,r,d,r0,mode)
    println(g.twist_angle*180/pi)
    println(size(g.sites,1))
    #println(typeof(g))
    remove_layer_g(g,1)
    println(size(g.sites,1))
    plot_R(g.sites)
end



#end  #end module

