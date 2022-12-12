############################### Description ############################################

# This script accompanies Smith et al., (in revision), "Tracing energy inputs into the 
# seafloor using carbonate sediments". The annotated code below reproduces Figure 3 
# in the main text using values from Table 1. The code for the physical model presented
# below is only slightly modified from the supplementary files in Trower et al., (2019) 
# "The Origin of Carbonate Mud." 

################### Load databases and import spreadsheet data #########################

# Clear any variables from the global environment.

  remove(list= ls())

# Load libraries. All of these must be downloaded to run the code.
 
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(forcats)
  library(readxl)
  library(ggpubr)

# Optional: Set the working directory. 

   setwd("C:/Users/bsmit/OneDrive/Desktop/Phanerozoic bioturbation")
  #setwd("...")
 
# Load data from Excel file and replace unused entries with NA

  data <- read_excel("Smith_etal_Dataset.xlsx", 
                   sheet = "data", 
                   na=c('','NA','#NAME?')) %>%
                   filter(grav>-5)

  data <- filter(data,!is.na(grav)|!is.na(crack_prop))

# Organize data for plotting

  DATAL <- data %>% select(OID,valid_binomial:general_group,grav,crack_prop) %>%
                     gather(model,value,grav,crack_prop)

############################## Define Variables ########################################

  D1 <- seq(from = 10, to = 1500, by = 5)           # range of grain sizes [um]
  
  D1 <- D1*10^-6;                                   # convert grain size to [m]

  ustar1 <-seq(from = 0.001, to = 0.15, by = 0.001) # range of shear velocities [m/s]

  grid <- expand.grid(ustar1,D1)                    # X-Y grid of different combinations
                                                    # of grain size and shear velocity

  f_int = 0.25                                      # Intermittency factor between 0 and
                                                    # 1 that represents the proportion 
                                                    # of time that grains are in motion 
                                                    # versus at rest. Default values of 
                                                    # 0.1 < f < 0.25 are taken from 
                                                    # modern studies of ooid shoals in 
                                                    # Trower et al., (2018)

  rho_s = 2850                                      # particle density in [kg/m^3]
  
  rho_f = 1035                                      # density of seawater in [kg/m^3]
  
  R = (rho_s - rho_f)/rho_f                         # submerged specific density 
                                                    # [unitless]
                                                    
  kv = 25*10^6                                      # coefficient that accounts for 
                                                    # differences between particles and 
                                                    # the bed surface (Scheingross et 
                                                    # al., 2014) [unitless]
                                                    # 
  g = 9.81                                          # gravitational acceleration[ m/s^2]
  
  nu = 1.3*10^-6                                    # kinematic viscosity of water 
                                                    # [m^2/s]
                                                    
  young = 144*10^9                                  # young's modulus [kg/m/s^2]
  
  strength = 1*10^6                                 # tensile strength of the bed 
                                                    # [kg/m/s^2]

  tauc = 0.03                                       # Critical Shields number.  
                                                    # 0.03 is good for sand.
                                                    
  gaurds2 = 1                                       # this sets limit to Ub if  = 1
  
  Stc = 10                                          # Critical stokes number below which
                                                    # collisions are viscously damped. 
                                                    # Default value of 10 is based on 
                                                    # experiments by Trower et al., 
                                                    # (2017) in EPSL.

  CE   <- 2.032577e+05                              # Value for chemical energy. See 
                                                    # main text for details.

  km = 0.92                                         # fraction of abrasion partitioned 
                                                    # to mud 
                                                   
  H = 50                                            # Water depth [m]. Matters more for 
                                                    # shallow flows than deep ones, 
                                                    # since grains in shallow flows may
                                                    # not reach terminal settling 
                                                    # velocity

################# Summarize bioturbation energy from spreadsheet ######################

  # Subset of complete data for crack propagation model and assign column names
  
  names <-c('group','rate_of_work', 'type')
  
  subset_crack <- data.frame(data$general_group,data$crack_prop)
  subset_crack$type <- 'crack propogation model'
  names(subset_crack) <- names
  
  subset_crack <-subset_crack[complete.cases(subset_crack$rate_of_work), ]
  subset_crack <- subset_crack%>%mutate(group = fct_reorder(group, 
                                                            desc(rate_of_work), 
                                                            .fun = 'median'))
  
  # Subset of complete data for gravitational model and assign column names
  
  subset_grav <- data.frame(data$general_group,data$grav)
  subset_grav$type <- 'gravitational model'
  names(subset_grav) <- names
  
  subset_grav <-subset_grav[complete.cases(subset_grav$rate_of_work), ]
  subset_grav <- subset_grav%>%mutate(group = fct_reorder(group, 
                                                          desc(rate_of_work), 
                                                          .fun = 'median'))
  
  # Merge subsetted data for plotting
  
  to_plot <- rbind(subset_grav,subset_crack)
  to_plot <- to_plot[complete.cases(to_plot$rate_of_work), ]
  
  to_plot <-to_plot%>%mutate(group = fct_reorder(group, 
                                                 desc(rate_of_work), 
                                                 .fun = 'median'))
  
  crack_sub <- subset(to_plot, to_plot$type == 'crack propogation model')
  grav_sub <- subset(to_plot, to_plot$type == 'gravitational model')
  
  crack_summary <-aggregate(rate_of_work~group, data = crack_sub, median)
  grav_summary <- aggregate(rate_of_work~group, data = grav_sub, median)
  
  Q25 <- quantile(grav_sub$rate_of_work, 0.25, na.rm = TRUE) # Lower quartile
  Q50 <- quantile(grav_sub$rate_of_work, 0.5, na.rm = TRUE)  # Median
  Q75 <- quantile(grav_sub$rate_of_work, 0.75, na.rm = TRUE)  # Upper quartile
  
####################### Calculate Settling Velocities ##################################

  CSF = 0.8                                         # Corey shape factor: 1 is for 
                                                    # spheres, 0.8 is for natural
  
  PS = 3.5                                          # 6 is for spheres, 3.5 is for 
                                                    # natural

  Dstar = (R*g*D1^3)/(nu^2)

  X = log10(Dstar)

  R1 = -3.76715+1.92944*X - 0.09815*(X^2) - 0.00575*(X^3) + 0.00056*(X^4)

  R2 = log10(1-((1-CSF)/0.85))-(((1-CSF)^2.3)*tanh(X-4.6)) + 0.3*
  (0.5-CSF)*((1-CSF)^2)*(X-4.6)

  R3 = (0.65-((CSF/2.83)*tanh(X-4.6)))^(1+((3.5-PS)/2.5))

  Wstar = R3*10^(R2+R1)

  ws1 = (R*g*nu*Wstar)^(1/3)

# Variables for impact rate eqn

  Rep  = (R*g*D1)^(1/2)*D1/nu                       # [unitless]
  A_GP = 1.3*10^-7                                  # constant from Garcia and Parker
  A1   = 0.36                                       # [unitless]
  Vp   = pi/6*D1^3                                  # [m^3]

  eps_v = kv*strength^2/(2*young)                   # kinetic energy per unit volume 
                                                    # eroded [kg/m/s^2]


# Pre-allocate empty vectors for computational efficiency in loops
 
  Ewi = rep(0,length(D1))
  c_b = matrix(0,length(D1),length(ustar1))
  Ir = matrix(0,length(D1),length(ustar1))
  Erate = matrix(0,length(D1),length(ustar1))
  Ewi = matrix(0,length(D1),length(ustar1))

##################### Calculate energy fluxes to the bed ############################### 

# Note: The nested for loops are from the original suspension/abrasion calculation. 
# For simplicity, the calculation is the same, and the units are adjusted in the 
# section after the for loop.

  for (position in 1:length(ustar1)){
  
    nn = position 
  
    ustar = ustar1[nn]
  
    Z = ustar/ws1*Rep^0.6 # [unitless]
    c_b[1:length(D1),nn] <- A_GP*Z^5/(1+A_GP/0.3*Z^5) # [unitless]

    Ir[1:length(D1),nn] <- A1*c_b[1:length(D1),nn]/Vp # impact rate (without w_i) 
                                                      # [1/m^3]

    for (point in 1:length(D1)){
  
       mm = point

       D = D1[mm]
       ws = ws1[mm]
       tau = ustar^2/(R*g*D)
       tstage = tau/tauc
       c_b1 = c_b[mm,nn]



    ##SUSP ABR CALC

      cdrag = (4/3)*(R*g*D)/(ws^2) # %[unitless]

      # compute flow velocity
      z0 = 3*D/30 # This is a roughness coefficient that needs to be set according to
                  # grainsize
                  
      dz = (H-z0)/1000 # [m]
      z=0
      z = seq(from = z0, to = H, by = dz) # [m]
      flow_vel = (ustar/0.41) * log(z/z0) # [m/s]
      flow_z = z # %[m]
      Uf = sum((ustar/0.41) * log(z/z0))*dz/H # [m/s]

      # compute bed load height and velocity
      hb = 0
      Us = 0
      
      if (tstage >1){
          hb = D*1.44*(tstage-1)^0.5                # height of the bed load layer [m]
          Us = (R*g*D)^0.5*1.56*(tstage-1)^0.56     # bed load velocity [m/s]
                    }
          Us1=Us
 
      # don't let particle velocity exceed fluid velocity

      if (gaurds2 == 1){
          Us[Us>Uf] <-Uf}

      # compute suspended load profile

      if (hb < H) {hb[hb<D]=D                      # sets the minimum height of bed 
                                                   # load layer [m]

          b = hb                                   # bottom of the suspended load layer 
                                                   # - same as height of bedload [m]
                                                   # 
          betta = 2;                               # Based on Scheingross et al. (2014) 
                                                   # best fit
                                                    
          P = ws/(0.41*ustar*betta)                # Rouse number [unitless]

       # This Log scale cleans up the integration
         
          di5 = 0
          i5=0
          res = 1000
          di5 = (log(H)-log(b))/res
          i5 = seq(from = log(b), to = log(H), by = di5) #[log(b):di5:log(H)]
          z=0
          dz=0
          z <- exp(i5)
          z[length(z)]<-H  # [m]
          dz = diff(z) #[m]
          dz = c(dz[1],dz)
          a1=sum((((1-(z[z>z0]/H))/(1-(b/H)))*(b/z[z>z0]))^P*log(z[z>z0]/z0)*dz[z>z0])/ 
            (Uf*H) * (ustar/0.41)
          cb = c_b1

        # find concentration profile
        
          c=0
          c = cb*(((1-(z/H))/(1-(b/H)))*(b/z))^P 
          c[1] = cb
          c[z==H]=0;

        # calculate the fall distance

          gradc <- -diff(c);
          gradc <- c(0,gradc)
          Hfall <- (1/cb)*sum(z*gradc) # [m]
          
          } else{
          hb = H
          cb = 1/(Us.*hb)
          Hfall = hb   
          a1=0
      }


      if (cb == 0)
      {Hfall = 0}


      if (P < 0.8){
      Hfall = H/5}

    # Probability for fluctuations
    
    sig = ustar                                    # [m/s]
    dx = sig/100                                   # the number of bins to subdivide 
                                                   # [m/s]
                                                   
    X = seq(from =-6*sig, to = 6*sig, by = dx)     # Spread distribution for six sigma 
                                                   # = w'/ws [m/s]
                                                   
    f = dnorm(X,0,sig)                             # centered at zero normal Gaussian 
                                                   # distribution [m/s]
                                                
    X = X/ws                                       # Normalize as w'/ws same as psi 
                                                   # [unitless]

    # Calculate impact velocity due to gravity and turbulence

    Scos = 1  #cosine of the angle of the bed for impacts.  Assume flat for ocean

    wfall = Scos*((2*(2/3)*D*g/cdrag*R)*
                   (1-exp(-cdrag*rho_f/rho_s*(Hfall/Scos)/(2/3*D))))^0.5 # [m/s]
   
    wfall[Hfall<=(0.5*D)]=0

    psifall = wfall/ws                             # [unitless]
    settlematrix = 0
    settlematrix = psifall + X                     # [unitless]
    settlematrix1=settlematrix
    settlematrix[settlematrix<0] = 0               # no negative impacts
    psifall_turb = sum((settlematrix)*f)*dx 
    psi_fall3 = sum((settlematrix^3)*f)*dx
    E1 = psi_fall3;                                # erosion with turbulence

    # Stokes number correction

    wi_st = settlematrix
    wi_st[(D*wi_st*ws*rho_s/(18*nu*rho_f))<Stc] = 0; 

    psi_fall3_st = sum((wi_st^3)*f)*dx
    E1_st = psi_fall3_st                           # erosion with turbulence and 
                                                   # stokes correction

    if (tstage <=1) {E1_st = 0}                    # Set erosion to very small if 
                                                   # particles stop moving


    if (Hfall == 0){
    E1_st = 0}

  
    Ewi[mm,nn] = E1_st*(g*D)^(3/2)                 # [m^3/s^3]

  }

# Gather results and finish calculation in correct units

  V_i = 1/2*Vp*rho_s/eps_v                         # volume eroded per impact 
                                                   # (without w_i) [m*s^2]

  Erate[1:length(D1),nn] <- V_i*Ir[1:length(D1),nn]*Ewi[1:length(D1),nn]
  
                                                   # Original calculation, gives unit 
                                                   # of [m^3/m^2/s]

 }

  Erate = (eps_v*Erate*60*60*24*365)               # Energy conversion to J/m^2/yr
  
  Erate = f_int*Erate                              # adjust for intermittent motions



  Erate[Erate == 0] <- NA
  reshape <- as.vector(t(Erate))
  grid$Var3 <- reshape

  grid <-na.omit(grid)
  
  colnames(grid) <- (c('Shear_velocity','Grain_size', 'KE_sed'))
  
  
########################### Chemical energy change #####################################
  
  # The constants listed below are from Table 1 in the text. We consider a case where
  # aragonite cements are growing in seawater at 25C.
  
  k = 12.9*10^(-6)*24*365                    # rate constant (mol m-2 yr-1) from 
                                             # Zhong and Mucci (1989)
  
  n = 2.26                                   # power law constant (-) from Zhong 
                                             # and Mucci (1989)
  
  M = .1                                     # molar mass of aragonite (kg mol-1)
  
  rho_c = 2850                               # density of aragonite (kg m-3)
  
  sigma = 283e-3                             # aragonite surface energy (J m-2) in 
                                             # water as calculated by Sun et al., 
                                             # (2015). Note that this is likely 
                                             # a conservative upper limit as Morse 
                                             # and Berner (1974) assume 0.08 J m-2
                                             # for calcite, and aragonite is more
                                             # soluble.
  R_f = 5
  grid$roughness = (6)*(1/rho_c)*
                   (1/grid$Grain_size)*R_f   # A roughness factor from Walter and Morse
                                             # (1984) based on echinoid fragments.
  
  omega = 5                                  # Use Omega_arag values of 3, 5, and 10.    
  
  rwk_z = 0.05                               # reworking depths (m), typically 2-5 cm.
                                             # Average for all Thayer data is 5 cm. 
  
  phi = .40                                  # porosity(m3 m-3)
  
  
  grid$S_0 = grid$roughness*(1-phi)*rho_c    # Initial reactive surface area (m2 m-3)
  
  grid$K_C <- sigma*2/3*grid$S_0^2*k*M*      # Run the chemical model for all grided
    rwk_z*(omega-1)^n /                      # values of Omega, roughness, and rwk_z.
    (rho_c) 
  
  grid$ratio <- grid$KE_sed/grid$K_C         # Add a column for an energy ratio 
                                             # between physical and chemical 
                                             # energy
  
  
############################ Panel 1: No Bioturbation ##################################
  
  
  ustar_draw_1 <- ws1/(2.5*.41*betta)             # based on the Rouse number, ~2.5 is 
  # the threshold  for suspension
  
  ustar_draw_2 <- ws1/(7.5*.41*betta)             # based on the Rouse number, ~7.5 is 
  # the threshold for washload
  
  ustar_draw <-data.frame(D1*10^6,ustar_draw_1,ustar_draw_2)
  
  p1 <-   ggplot() +  
    
    geom_tile(data = grid, 
              aes(x = Shear_velocity,
                  y = 1000000*Grain_size,
                  fill = ratio)) + 
    
    scale_fill_gradient2(
      breaks = c(1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3), 
      trans = "log",
      low = "blue",
      mid = 'yellow',
      high = 'red',
      limits = c(10^(-3), 10^(4))) +
    
    geom_contour(data = grid, aes(x = Shear_velocity,
                                  y = 1000000*Grain_size,
                                  z = ratio ), color = "black", breaks = c(1), size = 1.5)+
    
    ylim(min(1000000*D1),max(1000000*D1)) + xlim(0, max(ustar1))+ theme_classic()+
    labs (x = "Shear velocity (m/s)", y = "Grain size (um)", fill = 'Energy Ratio') +
    ggtitle('Sediment transport only') + theme(plot.title = element_text(hjust = 0.5))
  
############################ Panel 2: Bioturbation Histogram ###########################
  
  p2 <- ggplot(data = grav_sub, aes (x = rate_of_work)) + 

        geom_vline(aes(xintercept=Q25), linetype="dashed") +
        geom_vline(aes(xintercept=Q50), linetype="solid", size = 0.75) + 
        geom_vline(aes(xintercept=Q75), linetype="dashed") +
        theme_classic() + 
        geom_histogram(color = 'black', fill = 'white', alpha = 1) + 
        geom_histogram(color = 'black', fill = 'blue', alpha = 0.5) +
        ylim(0,20) + xlab('log rate of work') + 
        ggtitle('Bioturbation work rates \n (gravity model)') + 
        xlab('Log10 rate of work [J/m2yr]') + ylab('Count')  +
        theme(plot.title = element_text(hjust = 0.5))
    
    
 
############################ Panel 3 Q25 Bioturbation ##################################
  p3 <-   ggplot() +  
    
    geom_tile(data = grid, 
              aes(x = Shear_velocity,
                  y = 1000000*Grain_size,
                  fill = (KE_sed + 10^(Q25))/K_C)) + 
    
    scale_fill_gradient2(
      breaks = c(1e-3,1e-2,1e-1,1e0,1e1,1e2, 1e3), 
      trans = "log",
      low = "blue",
      mid = 'yellow',
      high = 'red',
      limits = c(10^(-3), 10^(4))) +
    
    geom_contour(data = grid, aes(x = Shear_velocity,
                                  y = 1000000*Grain_size,
                                  z = (KE_sed + 10^(Q25))/K_C), 
                                  color = "black", breaks = c(1), size = 1.5)+
    
    ylim(min(1000000*D1),max(1000000*D1)) + xlim(0, max(ustar1))+ theme_classic()+
    labs (x = "Shear velocity (m/s)", y = '', fill = 'Energy Ratio') +
    ggtitle('Transport + bioturbation \n (first quartile)') + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position = 'none')
  
  
############################ Panel 4 Q50 Bioturbation ##################################
  p4 <-   ggplot() +  
    
    geom_tile(data = grid, 
              aes(x = Shear_velocity,
                  y = 1000000*Grain_size,
                  fill = (KE_sed + 10^(Q50))/K_C)) + 
    
    scale_fill_gradient2(
      breaks = c(1e-3,1e-2,1e-1,1e0,1e1,1e2, 1e3), 
      trans = "log",
      low = "blue",
      mid = 'yellow',
      high = 'red',
      limits = c(10^(-3), 10^(4))) +
    
    geom_contour(data = grid, aes(x = Shear_velocity,
                                  y = 1000000*Grain_size,
                                  z = (KE_sed + 10^(Q50))/K_C), 
                 color = "black", breaks = c(1), size = 1.5)+
    
    ylim(min(1000000*D1),max(1000000*D1)) + xlim(0, max(ustar1))+ theme_classic()+
    labs (x = "Shear velocity (m/s)", y = '', fill = 'Energy Ratio') +
    ggtitle('Transport + bioturbation \n (median)') + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = 'none')
  
  
############################ Panel 5 Q75 Bioturbation ##################################
  p5 <-   ggplot() +  
    
    geom_tile(data = grid, 
              aes(x = Shear_velocity,
                  y = 1000000*Grain_size,
                  fill = (KE_sed + 10^(Q75))/K_C)) + 
    
    scale_fill_gradient2(
      breaks = c(1e-3,1e-2,1e-1,1e0,1e1,1e2, 1e3), 
      trans = "log",
      low = "blue",
      mid = 'yellow',
      high = 'red',
      limits = c(10^(-3), 10^(4))) +
    
    geom_contour(data = grid, aes(x = Shear_velocity,
                                  y = 1000000*Grain_size,
                                  z = (KE_sed + 10^(Q75))/K_C), 
                 color = "black", breaks = c(1), size = 1.5)+
    
    ylim(min(1000000*D1),max(1000000*D1)) + xlim(0, max(ustar1))+ theme_classic()+
    labs (x = "Shear velocity (m/s)", y = '', fill = 'Energy Ratio') +
    ggtitle('Transport + bioturbation \n (upper quartile)') + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = 'none')
  
  
########################## Composite Figures ##################################  
p_out <- ggarrange(ggarrange(p1, p2, nrow = 1, ncol = 2),
          ggarrange(p3, p4, p5, ncol = 3, nrow =1 ),
                    nrow = 2, align = 'v')
  
p_out
  
  ggsave('Figure_3.eps', device = 'pdf', plot = p_out, width = 17, 
         height = 12, units = c('cm'))
