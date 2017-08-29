################################################################################
##################### MODULES FOR MULTI PARTICLE COLLISIONS ####################
################################################################################
#This is to get the number of box where the particle is.
using StatsBase
function get_box(parts::Array{particle,1}, boxes::Array{box,1}, Lx::Int64)
  for box in boxes
    box.np[:] = 0
  end
  for p in parts
    p.indbox = ceil(p.pgrid[1]) + Lx * (ceil(p.pgrid[2])-1)
    if p.indbox < 0
      println(p)
    end
    boxes[p.indbox].np[p.tp] += 1
    # println(boxes[p.indbox].np[p.tp], p.tp, p.indbox)
  end
end
################################################################################
#rotation function
function rotate_vec(v::Array{Float64,1}, α::Float64)
  R = [cos(α) sin(α) ; -sin(α) cos(α)]
  return R*v
end
################################################################################
#computing the momentum of the boxes
function box_vel(parts::Array{particle,1},boxes::Array{box,1})
  for (i, box) in enumerate(boxes) #this is for enumerating the boxes
    box.vel = zeros(2)
    tmass = 0 #initializating the total mass of the box
    for p in filter(x-> x.indbox == i, parts) #the loop is made in the particles inside the box i
      #println("inside the filter ", i)
      box.vel += p.mass * p.vel #sums the momentums of all particles
      tmass += p.mass #sums the mass
    end
    if tmass != 0
      box.vel /= tmass #normalize the momentum of the box with the mass of the particles.
    end
  end
end
################################################################################
#This is for the collisions of the particles of equal mass or same particles
function box_velmc(parts::Array{particle,1},boxes::Array{box,1}, m::Array{Float64,1})
  for (i, box) in enumerate(boxes) #this is for enumerating the boxes
    if countnz(box.np) <= 1; continue; end
    for j in m #cycling over the different masses.
      box.vel = zeros(2)
      tmass = 0 #initializating the total mass of the box
      for p in filter(x-> x.indbox == i && x.mass == j , parts) #the loop is made in the particles inside the box i
        #println("inside the filter ", i)
        box.vel += p.mass * p.vel #sums the momentums of all particles
        tmass += p.mass #sums the mass
      end
      if tmass != 0
        box.vel /= tmass #normalize the momentum of the box with the mass of the particles.
      end

    end
  end
end
################################################################################
#collision in the new way.
function collide_mc(parts::Array{particle,1}, angles::Array{Float64})
  α = rand(angles)rand([-1,1])
  vel = zeros(2)
  tmass = 0 #initializating the total mass of the box
  for p in parts #the loop is made in the particles inside the box
    vel += p.mass * p.vel #sums the momentums of all particles
    tmass += p.mass #sums the mass
  end
  vel /= tmass #normalize the momentum of the box with the mass of the particles.
  for p in parts #loop over all particles
    v = p.vel - vel #extraction of the velocity of the box
    #α = 90.0*rand([-1,1]) #the rotation angle
    vn = rotate_vec(v, α) #rotation of the vector
    p.vel = vn + vel #adding the vector and the velocity of the box
  end
end
function collide_sc(parts::Array{particle,1}, tp::Int64, angles::Array{Float64})
  for j =1:tp
    α = rand(angles)*rand([-1,1])
    vel = zeros(2)
    tmass = 0
    sp = filter(x-> x.tp == j, parts)
    if isempty(sp); continue; end
    for p in sp
      vel += p.mass * p.vel
      tmass += p.mass
    end
    vel /= tmass
    for p in sp #loop over all particles
      v = p.vel - vel #extraction of the velocity of the box
      #α = 90.0*rand([-1,1]) #the rotation angle
      vn = rotate_vec(v, α) #rotation of the vector
      p.vel = vn + vel #adding the vector and the velocity of the box
    end
  end
end
#function of the shifting, actually you shift the positions of the particles.
function shift_grid!(parts::Array{particle,1},a::Float64, dim::Array{Int64,1})
  δx = rand()*rand(-a/2:a/2)
  δy = rand()*rand(-a/2:a/2)
  for p in parts
    if p.pos[1] + δx < 0 || p.pos[1] + δx > dim[1]
      p.pgrid[1] = p.pos[1]
    else
      p.pgrid[1] = p.pos[1] + δx
    end
    p.pgrid[2] = mod(p.pos[2] + δy, dim[2])
  end
end
#now we need to shift back the particles.
function shiftback_grid!(parts::Array{particle,1})
  for p in parts
    p.pgrid[1] = p.pos[1]
    p.pgrid[2] = p.pos[2]
  end
end
################################################################################
#computing the new velocities of the particles, collisions.
function parts_vels!(parts::Array{particle,1}, boxes::Array{box,1}, angles::Array{Float64,1})
  for p in parts #loop over all particles
    v = p.vel - boxes[p.indbox].vel #extraction of the velocity of the box
    α = rand(angles)*rand([-1,1]) #the rotation angle
    vn = rotate_vec(v, α) #rotation of the vector
    p.vel = vn + boxes[p.indbox].vel #adding the vector and the velocity of the box
  end
end
################################################################################
#getting new positions of the particles with periodic boundary conditions (pbc)
function getpos_pbc!(parts::Array{particle,1}, τ::Float64, dim::Array{Int64,1})
  for p in parts
    p.pos[1] = mod(p.pos[1] + p.vel[1] * τ , dim[1]) #to get the periodic boundary conditions.
    p.pos[2] = mod(p.pos[2] + p.vel[2] * τ , dim[2])
  end
end
################################################################################
#this is the ininitialization part, first the normalization of the total momentum.
function norm_momentum!(parts::Array{particle,1})
  vt = zeros(2) #sum of momentum
  mt = 0 #total mass
  for p in parts #summing all velocities and masses
    vt += p.mass * p.vel
    mt += p.mass
  end
  vt /= mt #normalizing the sum
  for p in parts #applying the normalization of momentum
    p.vel = p.vel - vt
  end
end
################################################################################
#now the normalization of temperature
function norm_temperature!(parts::Array{particle,1}, Tr::Float64)
  T = 0 #temperature
  for p in parts #this is the calculation of the Temperature
    T += p.mass * norm(p.vel)^2 / (2 * length(parts))
  end
  for p in parts
    p.vel = p.vel * sqrt(Tr / T) #multiplying each velocity for the factor of sqrt(Ttarget/Tactual)
  end
end
################################################################################
#this function compute probabilities of P(n(t+1)|n(t)) and returns how many particles are going to react.
function prob_box(np::Array{Int64}, kr::Float64)
  # println("func prob box np", np)
  nr = min(np[1],np[2]) #estimating minimum or maximum of reactions
  if nr != 0
    p = zeros(nr+1) #array of probabilities
    p[1] = 1 - kr #no-reaction case
    for i=1:nr
      p[i+1] = kr * ( factorial(big(np[1])) / factorial(big(np[1])-i) * factorial(big(np[2])) / factorial(big(np[2])-i) * factorial(big(np[3])) / factorial(big(np[3])+i) )
    end
    s = sum(p)
    if s != 1
      p = p / s
    end
    cs = cumsum(p); a = rand()
    out = findfirst(sort([a;cs]),a)-1
  else
    out = 0
  end
  return out
end
################################################################################
function col_box(box::box,parts::Array{particle,1}, m::Array{Float64,1})
  parbox = filter(x-> x.indbox == box.ind, parts) #selecting the particles that are in the box.
  #println(i,'\t',length(parbox))
  if isempty(parbox) == false #ignore next steps if the box is empty
    collide_mc(parbox)
    if countnz(box.np) > 1
      collide_sc(parbox, m)
    end
  end
end
################################################################################
#next function is for the birth and death process of the reaction.
function reac_box!(parts::Array{particle,1}, nr::Int64, H::Float64)
  pa = filter(x -> x.tp == 1, parts)
  pb = filter(x -> x.tp == 2, parts)
  # println(length(pa), " pa=", pa,"\n", length(pb), " pb =", pb, "\nnr= ", nr)
  pc = Array{particle,1}()
  for i =1:nr #nr = number of reactions
    a = rand(pa); b = rand(pb) # choosing particles a and b at random
    cm = a.mass + b.mass #c mass
    x = a.mass/cm * a.pos[1] + b.mass/cm * b.pos[1] #position of new particle
    y = a.mass/cm * a.pos[2] + b.mass/cm * b.pos[2]
    vx = a.mass/cm * a.vel[1] + b.mass/cm * b.vel[1] #velocity of the new particle, conserving the momentum
    vy = a.mass/cm * a.vel[2] + b.mass/cm * b.vel[2]
    c = particle([1,1], cm, 3) #creating the particle c
    c.pos = [x,y]; c.vel = [vx, vy] #positions and velocities.
    push!(pc, c) #adding the particle c to the
    a.mass = 0.0; b.mass = 0.0
    filter!(x -> x.mass != 0.0, pa); filter!(x -> x.mass != 0.0, pb)  #removing the particles that reacted.
  end
  #  println(pa,"\t", pb, "\n pc =",pc)
  ################################## Enthalpy ############3333333
  H_tot= length(pc)*H
  va_enthalpy=sqrt(H_tot)*1.73;
  vb_enthalpy=sqrt(H_tot)*1.22;
  vc_enthalpy=sqrt(H_tot)*1.00;
  ################## Velocity increment for pa ############################
  for p in pa
    # println("flag Enthalpy gain")
    # println("velo1 $(p.vel)")
    p.vel = va_enthalpy*p.vel
    # println("velo2 $(p.vel)")
  end
  ################## Velocity increment for pb ############################33333
  for p in pb
    # println("flag Enthalpy gain")
    p.vel = vb_enthalpy*p.vel
  end
  ################# Velocity increment for pc ############################33333
  for p in pc
    # println("flag Enthalpy gain")
    p.vel = vc_enthalpy*p.vel
  end
  filter!(x -> x.mass != 0.0, pc)
  # println("function reation end # pa= ",length(pa)," pb=  ",length(pb)," pc= ",length(pc))
  return pc
end
################################################################################
#next function is to graph the positions
function grap_pos(parts::Array{particle,1},tp::Int64)
  x = filter(x -> x.tp == tp, parts)
  xy = zeros(length(x),2)
  for (i,p) in enumerate(x)
    xy[i,:] = collect(p.pos)
  end
  return xy
end
################################################################################
#This is the nucleation function
# parts=per_c_box
function nucleate_c(parts::Array{particle,1}, ks::Int64)
  pc = filter(x -> x.tp == 3, parts)
  pd = filter(x -> x.tp == 4, parts)
  nc = length(pc);
  no_seed = div(length(pc),ks)
  pd=Array{particle,1}()
  for i=1:no_seed # number of seed particles
    nc=length(pc)
    # println("i= $i and no_seed =$no_seed, nc =$nc, k =$ks")
    # for p in pc; println(p.pos, p.mass);end
    nr_seed  = sample(pc,ks, replace = false) #selection of ks random particles
    mc_seed=nr_seed[1].mass * ks
    d = particle([1,1], mc_seed, 4)
    x1=0.0;y1=0.0;v1=0.0;v2=0.0
    for p in nr_seed
      x1+=p.pos[1];v1+=p.vel[1];
      y1+=p.pos[2];v2+=p.vel[2];
    end
    x2=x1/length(nr_seed);    v11=v1/length(nr_seed);
    y2=y1/length(nr_seed);    v22=v2/length(nr_seed);
    d.pos=[x2, y2]; d.vel=[v11, v22];
    push!(pd,d)
    for p in nr_seed;  p.mass=0.0; end
    filter!(x -> x.mass != 0.0, pc)
  end
  return pd
end
############### AGGREGATION  #######################################################
function aggre_sbox(parts::Array{particle,1}, k_agg::Int64)
  pc = filter(x -> x.tp == 3, parts)
  pd = filter(x -> x.tp == 4, parts)
  total_massD=total_mass(pd)
  pd1 = Array{particle,1}()
  pcmass=total_mass(pc);
  nd_seed = rand(pd)
  nd_mass = nd_seed.mass;
  md_aggr = pcmass + nd_seed.mass;
  d_aggr = particle([1,1], md_aggr, 4)
  v1=0.0;v2=0.0;
  for p in pc; v1+=p.vel[1]; v2+=p.vel[2];  end
  v1=v1+nd_seed.vel[1];
  v2=v1+nd_seed.vel[2];
  v11=v1/(length(pc)+1);
  v22=v2/(length(pc)+1);
  d_aggr.vel=[v11, v22]; d_aggr.pos = nd_seed.pos    #   [x1, x2];
  for p in pc;  p.mass = 0.0; end
  nd_seed.mass = 0.0;
  push!(pd1,d_aggr)
  filter!(x -> x.mass != 0.0, pc)
  filter!(x -> x.mass != 0.0, pd)
  return pd1
end
###############################################################################
function wall_nd(parts::Array{particle,1})
  pd = filter(x -> x.tp == 4, parts)
  pd_wall = filter(x -> x.vel == [0.0,0.0], pd)
  filter!(x -> x.vel == [0.0,0.0], pd)
  return pd_wall
end
function getpos_slip!(parts::Array{particle,1}, τ::Float64, dim::Array{Int64,1}, vloss_wall::Float64, ks::Int64, ND::Int64)
  # pe=Array(particle,1)()
  d_part = filter(x -> x.tp == 4, parts)
  c_part = filter(x -> x.tp == 3, parts)
  mass_d=total_mass(d_part)
  for p in  parts
      #doing the bouncing on the x axis
    if p.pos[1] + p.vel[1] * τ > dim[1]
      dif = p.pos[1] + p.vel[1] * τ - dim[1]
        p.pos[1] = dim[1] - dif #it moves back
        p.vel[1] = -p.vel[1] #changing direction of velocity in x
      elseif p.pos[1] + p.vel[1] * τ < 0
        dif = abs(p.pos[1] + p.vel[1] * τ)
          p.pos[1] = 0 + dif
          p.vel[1] = -p.vel[1] #changin direction of velocity in x
        else
          p.pos[1] = p.pos[1] + p.vel[1] * τ
        end
        #######################now the bouncing on the y axis
        if p.pos[2] + p.vel[2] * τ > dim[2]
          dif = p.pos[2] + p.vel[2] * τ - dim[2]
            p.pos[2] = dim[2] - dif
            p.vel[2] = -p.vel[2]*vloss_wall
          elseif p.pos[2] + p.vel[2] * τ < 0
            dif = abs(p.pos[2] + p.vel[2] * τ)
              p.pos[2] = 0 + dif
              p.vel[2] = -p.vel[2]*vloss_wall
            else
              p.pos[2] = p.pos[2] + p.vel[2] * τ
            end
          end
          # Totalmass_D=total_mass(d_part); println(Totalmass_D)
          # if Totalmass_D > 1500; println("nd = $(length(d_part))")
          #   for (i, box2) in enumerate(boxes)
          #     parabox = filter(x-> x.indbox == i, d_part); Mass = total_mass(parabox)
          #     if Mass > 50;println("totalmass_D = $Totalmass_D, mass= $Mass, i=$i, length = $(length(parabox))")
          #       for  p in parabox;p.vel[1] = 0.0; p.vel[2]=0.0;end
          #     end
          #   end
          # end
          ##################################################
          # pd_walls1 = filter(x -> x.vel == [0.0,0.0], d_part)
          # if !isempty(pd_walls1);push!(partsd1, pd_walls1...) ; end;
          # filter!(x -> x.vel != [0.0,0.0], parts);
            # if total_massD > 400  # && p.mass >= c_part[1].mass*(ks-6)
            #   for p in d_part
            #     Mass = Mass
            #   # p.pos[2] = dim[2] #it stick there
            #   p.vel[1] = 0.0;p.vel[2] = 0.0;
            #   # println("pdmass ", p.mass,  "\t",p.indbox,"\t",p.pos)
            # end
          pd_walls = filter(x -> x.vel == [0.0,0.0], parts)
          filter!(x -> x.vel != [0.0,0.0], parts)
          return pd_walls
        end
  ##########################33
        function total_mass(parts::Array{particle,1})
          mass=0.0
          for p in parts
            mass = mass + p.mass
          end
          return mass
        end
        ################################################################################
        function zero_velo(parts::Array{particle,1}, boxes::Array{box,1} )
          pd11 = filter(x-> x.tp == 4, parts); Totalmass_D=total_mass(pd11); println(Totalmass_D)
          if Totalmass_D > 1500; println("nd = $(length(pd11))")
            for (i, box2) in enumerate(boxes)
              parabox = filter(x-> x.indbox == i, pd11); Mass = total_mass(parabox)
              if Mass > 50;println("totalmass_D = $Totalmass_D, mass= $Mass, i=$i, length = $(length(parabox))")
                for  p in parabox; p.vel[1] = 0.0; p.vel[2]=0.0;end
              end
            end
          end
          pd_walls = filter(x -> x.vel == [0.0,0.0], pd11)
          filter!(x -> x.vel != [0.0,0.0], parts)
          return pd_walls
        end
