#####################           TEST            ################################
#This program is a test code for the Liesegang package, it is a simple model based
#on the Stochastic Rotation Dynamics model.
#Here the idea is to have two species of particles, A and B
#loading the package Liesegang.jl
include("Liesegang.jl")
include("types.jl")
include("func_sim.jl")
using StatsBase
using Plots
pyplot() #plotting packaty
#plot(sin(x))
################################################################################
#defining the parameters
Lx = 100
Ly = 10 #size of the space
dim = [Lx,Ly]
a = 1.0 #size of the boxes, default = 1
m = [1.0, 2.0] #masses
np = [500,500] #number ojuf particles
number_per_cell=Int64(ceil((np[1]+np[2])/(Lx*Ly)))
println("number_per_cell = ",number_per_cell,", Lx = ", Lx)
println("Ly = ", Ly, ", np1 = ", np[1], ", np2 = ", np[2])
ntp = 4 #number of species.
Tr = 1/3 #reference temperature
τ = 1.73 #1.73
kr = 0.9 # probability of reaction or reaction rate.
tmax = 10;  II = 0; NC_th = 500; ND_th = 40;
k_nuc = 10; ## Condn for first nucleation
k_nuc1 = 1; ### condn for nucleation in near by cells
k_agg = 1; ### condn for aggregation in same cell
H = 1.2  #### Enthalpy
vloss_wall = 0.8;
angles = [60.0,90.0,120.0,180.0]
################################################################################
###########                       INITIALIZING                       ###########

parts = vcat([particle(1,dim, m[1],1) for _ in 1:np[1]], [particle(Lx,dim,m[2],2) for _ in 1:np[2]]) #initializing the particles.
#normalizing the momentum
norm_momentum!(parts)
#now the temperature to the reference Tr
norm_temperature!(parts, Tr)
#initializing the boxes
boxes = [box(ntp,i) for i in 1:(Lx * Ly)]
# print(length(boxes), boxes)
partsd1=Array{particle,1}()
# parts_nuc=Array{particle,1}()
# parts_agg=Array{particle,1}()
################################################################################
#########################    now the simulation...   ###########################
# tic()
anim = @animate for t in 1:tmax
  #############    streaming step   ###########################################
  pa=filter(x->x.tp==1, parts);pb=filter(x->x.tp==2, parts);pc=filter(x->x.tp==3, parts)
  pd=filter(x->x.tp==4, parts);T_mass1 = total_mass(parts); T_mass2 = total_mass(partsd1); T_mass=T_mass1 + T_mass2
  ff=open("data.d", "a");println(ff, "Time = $t")
  NC=length(pc);  ND=length(pd);
  if t>0;  println("\n\n#############################      Time = $t         ##########################################")
    println("NA= ",length(pa),'\t',"NB= ",length(pb),'\t',"NC= ",length(pc),'\t',"ND= ",length(pd),"\t","NDD= ", length(partsd1),"\t","M= ",T_mass,);end
    mass_pa=total_mass(pa);mass_pb=total_mass(pb);mass_pc=total_mass(pc);mass_pd=total_mass(pd);
    mass_pc_pd = mass_pc + mass_pd;
    println("Mass: pa=$mass_pa, pb=$mass_pb, pc=$mass_pc, pd=$mass_pd, pc+pd=$mass_pc_pd")
    ############################  for continuous rate of particles   #############################
    close(ff)
    if II == 1
      NA= filter(x->x.pos[1] < 1, parts);NB= filter(x->x.pos[1] > (Lx-1), parts);
      na=length(NA);nb=length(NB);println("length of parts ", length(parts), ", NA= $na and NB = $nb")
      if na < np[1];  diff1=np[1]-na;  parts_a = [particle(1,dim, m[1],1) for _ in 1:diff1];push!(parts,parts_a...);end
      if nb < np[2];  diff2=np[2]-nb;  parts_b = [particle(Lx,dim, m[2],2) for _ in 1:diff2];push!(parts,parts_b...);end
      NA1= filter(x->x.pos[1] < 1, parts);NB1= filter(x->x.pos[1] > (Lx-1), parts);
    end
    # na1=length(NA1);nb1=length(NB1);println("NA1= $na1 and NB1 = $nb1")
    # println("length of parts ", length(parts));
    ###################### streaming #####################################
    pd_walls=getpos_slip!(parts, τ, dim, vloss_wall, k_nuc, ND)
    if !isempty(pd_walls);push!(partsd1, pd_walls...); end
    #first the grid is shifted
    shift_grid!(parts, a, dim)
    # @time shift_grid!
    #now label the particles in the boxes
    get_box(parts, boxes, Lx)
    # pause()
    ############################################################
    # pd_walls = zero_velo(parts, boxes)
    # if !isempty(pd_walls);push!(partsd1, pd_walls...); end
    ################################################################
    print("reaction time "); tic();
    for (i, box1) in enumerate(boxes) #cycling over the boxes
      parbox = filter(x-> x.indbox == i, parts) #selecting the particles that are in the box.
      # filter!(x -> x.mass != 0.0, parbox)
      if isempty(parbox); continue; end  #ignore next steps if the box is empty (come out form the for loop )
      collide_mc(parbox, angles)
      if countnz(box1.np) > 1 #this is if there is more than one type of particle in the box
        collide_sc(parbox, ntp, angles)
        nr = prob_box(box1.np, kr) #compute the probabilites of the transitions
        if nr != 0
          pc = reac_box!(parbox, nr, H);
          push!(parts, pc...)
        end
      end
    end
    toc();
    ###############################################################################3
    print("nucleation time ");tic();
    if NC > NC_th
      for (i, box2) in enumerate(boxes)
        filter!(x -> x.mass != 0.0, parts)
        parbox = filter(x -> x.indbox == i, parts)
        if box2.np[3] >= k_nuc && box2.np[4] >= 0
          # println("\n########################### main nucleation ")
          println("Nucleation ", box2.np, " index = $i")
          pd = nucleate_c(parbox, k_nuc)
          if isempty(pd);continue;end
          push!(parts, pd...)
        end
        #########################################################################
        if box2.np[3] >= k_agg && box2.np[4] > 0 && ND > ND_th
          # println("\n########################### main AGGREGATION ")
          println("AGGREGATION ", box2.np, " index = $i")
          pd1 = aggre_sbox(parbox, k_agg)
          push!(parts, pd1...)
        end
        if box2.np[4] > 0 &&  ND > ND_th
          for j=[-Lx,-1,+1,Lx]; i1 = i+j;
            if 0 < i1 < Lx*Ly
              parbox1 = filter(x-> x.indbox == i1, parts)
              pc = filter(x-> x.tp == 3, parbox1)
              pd = filter(x-> x.tp == 4, parbox1)
              # println("flag1 auto catalyst , i= $i, i1= $i1")
              nc=length(pc);nd=length(pd)
              if length(pc) >= k_nuc1 && length(pd) == 0
                println("Catalyst ", box2.np, " index = $i", " index1 = $i1")
                pd2 = nucleate_c(pc, k_nuc1)
                push!(parts,pd2...)
              end
            end
          end
        end
      end
    end
    toc();
    ################################################################################
    pd11 = filter(x-> x.tp == 4, parts); Totalmass_D=total_mass(pd11);# println("pdmass = ",Totalmass_D)
    if Totalmass_D > 700 && ND > ND_th;
      for p in pd11
        p.pgrid[1] = p.pos[1]
        p.pgrid[2] = p.pos[2]
      end
      get_box(pd11, boxes, Lx)
      # println("length nd = $(length(pd11))")
      for (i, box2) in enumerate(boxes)
        parabox = filter(x-> x.indbox == i, pd11); Mass = total_mass(parabox)
        if Mass > 30;println("totalmass_D = $Totalmass_D, mass= $Mass, i=$i, length = $(length(parabox))")
          for  p in parabox; p.vel[1] = 0.0; p.vel[2]=0.0;end
        end
      end
    end
    ##################################################
    pd_walls1 = filter(x -> x.vel == [0.0,0.0], pd11)
    if !isempty(pd_walls1);push!(partsd1, pd_walls1...) ; end;
    filter!(x -> x.vel != [0.0,0.0], parts);
    ##############################################################################
    filter!(x -> x.mass != 0.0, parts)
    # get_box(parts, boxes, Lx)
    ########### shifting back the particles to their original places ############
    shiftback_grid!(parts)
    x = grap_pos(parts,1)
    y = grap_pos(parts,2)
    z = grap_pos(parts,3)
    z1= grap_pos(parts,4)
    z2= grap_pos(partsd1,4)
    #vx = [parts[i].vel[1]/3 for i in 1:np] #dividing the vectors by a factor of 3 just for the visualization.
    #vy = [parts[i].vel[2]/3 for i in 1:np]
    scatter(x[:,1],x[:,2], xlims = (0,Lx), ylims = (0,Ly), size = (1200,200),legend=false)
    # scatter(x[:,1],x[:,2], xlims = (0,Lx), ylims = (0,Ly), size = (Lx*(1.25),Ly*30), legend=false)
    scatter!(y[:,1],y[:,2])# xlims = (0,Lx), ylims = (0,Ly))
    scatter!(z[:,1],z[:,2], title= "t=$(tmax), knu=$(k_nuc), kagg=$(k_agg), k_nuc1=$(k_nuc1), H=$(H), vloss=$(vloss_wall), ppc=$(number_per_cell), ND_th= $ND_th, NC_th=$NC_th, np=$(np[1]), t=$t")
    scatter!(z1[:,1],z1[:,2], c=:yellow, xlabel="N=fixed")
    scatter!(z2[:,1],z2[:,2], c=:black)

    # println(length(pa),'\t',length(pb),'\t',length(pc),'\t',length(pd),"\t",T_mass, "\t\t\t\tt = ",I)
  end
  gif(anim, "z_t_$(tmax)_knu_$(k_nuc)_kagg_$(k_agg)_k_nuc1_$(k_nuc1)_H_$(H)_vloss_$(vloss_wall)_ppc_$(number_per_cell).gif", fps = 9)
    # toc()
