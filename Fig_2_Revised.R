############################### Description ############################################

  # This script accompanies Smith et al., (submitted), "Tracing energy inputs into the 
  # seafloor using carbonate sediments". The annotated code below reproduces Figure 2 
  # in the main text using data from the supplementary Excel file. 

################### Load databases and import spreadsheet data #########################

  # Clear any variables from the global environment.

    remove(list= ls())

  # Load libraries. All of these must be downloaded to run the code. 

    library(tidyr)
    library(dplyr)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(forcats)
    library(readxl)

  # Optional: Set the working directory. The working directory must include 
  # the input file(s) (with headers) saved as a .csv file. 

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
  
########################### Bioturbation energy calculation ############################

  # This section uses the Thayer data to calculate energy associated with crack 
  # propagation and excavation mechanisms. For convenience, these numbers are already
  # included in the supplemental spreadsheet, but the section below shows how they were 
  # calculated using values in other columns. 

  # Unit conversions into mks units:
    
    data$log10_volume_m3    <- as.numeric(data$log10_volume) - 9     
                                               # units conversion: mm^3-> m^3
    
    data$longest_axis_m     <- as.numeric(data$longest_axis_mm)*10^-3    
                                               # unit conversion: mm -> m
    
    data$narrowest_axis_m     <- as.numeric(data$narrowest_axis_mm)*10^-3    
                                               # unit conversion: mm -> m
    
    data$log10_max_length_m <- as.numeric(data$log10_max_length) - 3  
                                               # unit conversion: mm -> m
    
    data$thayer_rework      <- as.numeric(data$pop.rwk.ann.literxm2_yr)*(10^-3)
                                               # unit conversions: 
                                               # liters/m2/year-> m^3/m^2/year
    
    dim <- dim(data)
    n = dim[1]
    G_c <- 385                                 # Fracture toughness of sediment 
                                               # in Pa*m^0.5, taken from Dorgan et al.,
                                               # (2011) after Johnson et al. (2002)
    
    phi = 0.4                                  # Assume 40% porosity
    
    
    data$radius <- ((3*(10^data$log10_volume_m3))/
                      (4*pi*10^data$log10_max_length_m))^0.5 
                                               # If a body volume and maximum length are
                                               # reported rather than a radius, then we 
                                               # assume an ellipsoid with a circular 
                                               # cross section. Under this assumption, 
                                               # we can solve for the radius of the 
                                               # burrow geometrically

    
    data$flag_1 <- is.na(data$radius)          # set flags for missing entries for this
                                               # method
    
                                               # Check to see if there are other ways 
                                               # to calculate volume
    
    for (i in 1:n){
      
      row = i
      
      if (data[row, "flag_1"] == "TRUE")
        
      {data[row, "radius"] <- .5*data[row, "narrowest_axis_m"]}
    }
    
    
    # Make new columns for the crack propagation and gravitational work models     
    
    data$crack_prop <- log10(2*G_c*data$thayer_rework/
                               (pi*data$radius*data$ind.rwk.z.cm*10^(-2)))
    
    data$grav<- log10(.5*(1- phi)*(2650-1035)*9.81*data$ind.rwk.z.cm*10^(-2)*
                        data$thayer_rework)
    
   # write.csv(data,
           #   "C:/Users/bsmit/OneDrive/Desktop/Phanerozoic bioturbation\\Thayer_updated_out.csv", row.names = FALSE)
    
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
    
######################### Bioturbation energy calculation ##############################

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

    roughness = c(11.6, 41, 450)               # surface area per mass (m2 kg-1) 
                                               # for 513 um coral, 81 um echinoids, and 
                                               # 5 um calcite rhombs after Walter and 
                                               # Morse (1984) 

    omega = c(3,5,10,16)                       # Use Omega_arag values of 3, 5, and 10.    

    rwk_z = c(0.02, 0.05)                      # reworking depths (m), typically 2-5 cm.
                                               # Average for all Thayer data is 5 cm. 

    phi = .40                                  # porosity(m3 m-3)


    grid <- expand.grid(rwk_z = rwk_z,         # Make a grid with different combinations
                        roughness = roughness, # Omega, roughness, and rwk_z. 
                        omega = omega)

    grid$Vc = 1*1*grid$rwk_z*(1-phi)           # Carbonate volume (m3)

    grid$S_0 = grid$roughness*(1-phi)*rho_c    # Initial reactive surface area per unit
                                               # volume (m2 m-3)

    grid$CE <- sigma*2/3*grid$S_0^2*k*M*       # Run the chemical model for all grided
              grid$rwk_z*(grid$omega-1)^n /    # values of Omega, roughness, and rwk_z.
              (rho_c) 

######################## Plot results and save output file #############################
   
   # Assign plot elements and handles
   
   P <- DATAL %>%
     group_by(general_name,model) %>%
     mutate(grpmed=median(value,na.rm=T)) %>%
     ungroup() %>%
     mutate(xax=fct_reorder(general_name,grpmed)) %>%
     ggplot(aes(xax,value,fill=model)) +
     
     annotate(geom = "rect", 
              ymin = log10(grid$CE[2]), xmax = Inf,       # Band for 500 um sand
              ymax = log10(grid$CE[14]), xmin = -Inf,
              fill = "green", alpha = 0.2) +
     

     
     annotate(geom = "rect",                              # Band for 81 um vf sand
              ymin = log10(grid$CE[4]), xmax = Inf,
              ymax = log10(grid$CE[16]), xmin = -Inf,
              fill = "yellow", alpha = 0.2) +
     
     
  
     annotate(geom = "rect", 
             ymin = log10(grid$CE[6]), xmax = Inf,
              ymax = log10(grid$CE[18]), xmin = -Inf,
              fill = "orange", alpha = 0.2) +             # Band for 5 um mud
     
  
     
     geom_boxplot() + theme_classic() + 
     facet_grid(~general_group,switch='x',scales='free_x',space='free_x',
                labeller=label_wrap_gen()) +
     scale_y_continuous(breaks=seq(-10,10,2)) +
     labs(y = 'Log10 Energy Flux\n(J m-2 yr-1',x='') +
     theme(axis.text.y=element_text(size=8),
           panel.grid.major.x=element_blank(),
           panel.grid.minor.x=element_blank(),
           panel.grid.minor.y=element_blank(),
           axis.text.x=element_text(angle=45,
                                    hjust=1,vjust=1),
           legend.position=c(0.2,0.1),
           axis.line.y.left=element_line(),
           strip.text=element_text(size=8,hjust=0.5),
           panel.border=element_blank(),
           panel.spacing=unit(0.25, "lines"),
           strip.background=element_blank(),
           strip.placement="inside")
   
    q <- ggplotGrob(P)
    
    lg <- linesGrob(
          x=unit(c(0,1),"npc"), y=unit(c(1,1),"npc"), 
          gp=gpar(col="black", lwd=4))
    
    for (k in grep("strip-b",q$layout$name)) {
          q$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
    }
    
grid.arrange(q) 

ggsave('Figure_2.pdf', plot = last_plot(), width = 30, 
       height = 14, units = c('cm'), dpi = 300)