# Guillaume Lobet - University of Liege


# ROOT-FIT
# Julkowska et al. (2014). 
# Capturing Arabidopsis root architecture dynamics with ROOT-FIT reveals 
# diversity in responses to salinity. 
# Plant Physiology. doi:10.1104/pp.114.248963


shinyServer(
  
  function(input, output) {
    
    #------------------------------------------------------
    # LOAD THE USER DATA
    
    Data <- reactive({
      
      inFile <- input$data_file   
      if (is.null(inFile)) return(NULL)
      data <- read.csv(inFile$datapath)  
      return(data)
      
    })   
    
    #------------------------------------------------------
    # PROCESS THE DATA
    
    Results <- reactive({
      
      if(input$runROOTFIT == 0){return()}
    
      
      rs <- Data()      
#       rs <- read.csv("~/Desktop/test_data.csv")
      
      if(input$software == "RSML_reader"){
        # Get the data in correct shape
        rs$length <- as.numeric(rs$length)
        rs$image <- as.character(rs$image)
        
        # Parse the name
        # For this step, we need to have a standart naming procedure
        rs$genotype <- "-"
        rs$treatment <- "-"
        rs$dag <- 0
        for(i in 1:nrow(rs)){
          name <- strsplit(strsplit(rs$image[i], split = "\\.")[[1]][1], split=input$name_sep)[[1]]
          for(j in 1:length(name)){
            if(grepl(input$gen_tag, name[j])) rs$genotype[i] <- gsub(input$gen_tag, "", name[j])
            if(grepl(input$tr_tag, name[j])) rs$treatment[i] <- gsub(input$tr_tag, "", name[j])
            if(grepl(input$dag_tag, name[j])) rs$dag[i] <- as.numeric(gsub(input$dag_tag, "", name[j]))
          }
        }
        
        rs$dag <- rs$dag - min(rs$dag)
        
        # Create lists for the different variables
        l.images <- unique(rs$image) # Get the list of images
        l.roots <- unique(rs$root) # Get the list of roots
        l.prim.roots <- unique(rs$root[rs$root_order == 0]) # Get the list of primary roots
        l.dag <- unique(rs$dag) # Get the list of days
        l.genotype <- unique(rs$genotype) # Get the list of genotypes
        l.treatment <- unique(rs$treatment) # Get the list of treatments
        
        # Create a new table for each primary root and each dag
        root <- ddply(rs[rs$root_order == 0,], .(root, dag), summarize, 
                           image = image,
                           genotype = unique(genotype)[1],
                           treatment = unique(treatment)[1])
        
        
        # For each of the primary, get the lateral root length and the number of laterals.
        for(d in l.dag){
          for(p in l.prim.roots){
            root$prim_length[root$dag == d & root$root == p] <- sum(rs$length[rs$root == p & rs$dag == d])
            root$tot_lat_length[root$dag == d & root$root == p] <- sum(rs$length[rs$parent == p & rs$dag == d])
            root$mean_lat_length[root$dag == d & root$root == p] <- mean(rs$length[rs$parent == p & rs$dag == d])
            root$lat_number[root$dag == d & root$root == p] <- length(rs$length[rs$parent == p & rs$dag == d]) 
          }
        }
        # For each of the primary, ge tthe total root length
        root$tot_length <- root$prim_length + root$tot_lat_length
        
        for(i in 6:ncol(root)){
          root[,i][is.nan(root[,i])] <- 0
        }
        
        # For each primary, get the growth for the different variables
        for(p in l.prim.roots){
          temp <- root[root$root == p,]
          root$prim_growth[root$root == p] <- temp$prim_length - temp$prim_length[1]
          root$lat_growth[root$root == p] <- temp$mean_lat_length - temp$mean_lat_length[1]
          root$lat_number_growth[root$root == p] <- temp$lat_number - temp$lat_number[1]
          root$tot_growth[root$root == p] <- temp$tot_length - temp$tot_length[1]
        }

    # ----------------------------------------------------------------------

    }else{
                
        root <- data.frame(root = paste(rs$Plate.name, "-", rs$Media, "-", rs$Genotype,"-", rs$Plant.id, sep=""))
        root$genotype <- factor(rs$Genotype)
        root$treatment <- factor(rs$Media)
        root$dag <- rs$Age - min(rs$Age)
        root$image <- rs$Picture.name
        # Create lists for the different variables
        l.images <- unique(root$image) # Get the list of images
        l.roots <- unique(root$root) # Get the list of roots
        l.prim.roots <- unique(root$root) # Get the list of primary roots
        l.dag <- unique(root$dag) # Get the list of days
        l.genotype <- unique(root$genotype) # Get the list of genotypes
        l.treatment <- unique(root$treatment) # Get the list of treatments
        

        # Create a new table for each primary root and each dag
        root$prim_length <- rs$MR.path.length
        root$lat_number <- rs$Number.LR.MR
        root$tot_lat_length <- rs$Total.root.size - root$prim_length        
        root$mean_lat_length <- 0
        for(i in 6:ncol(root)){
          root[,i][is.nan(root[,i])] <- 0
          root[,i][is.na(root[,i])] <- 0
        }
        root$mean_lat_length[root$lat_number > 0] <- root$tot_lat_length[root$lat_number > 0] / root$lat_number[root$lat_number > 0]        
        root$tot_length <-  rs$Total.root.size 
       
        
        # For each primary, get the growth for the different variables
        for(p in l.prim.roots){
          temp <- root[root$root == p,]
          root$prim_growth[root$root == p] <- temp$prim_length - temp$prim_length[1]
          root$lat_growth[root$root == p] <- temp$mean_lat_length - temp$mean_lat_length[1]
          root$lat_number_growth[root$root == p] <- temp$lat_number - temp$lat_number[1]
          root$tot_growth[root$root == p] <- temp$tot_length - temp$tot_length[1]
        }
        
    }


        # Create a new table for each primary in order to retrieve the different factors
        growth.data <- ddply(root, .(root), summarize, 
                             genotype = unique(genotype)[1],
                             treatment = unique(treatment)[1])
        
        
# ------------------------------------------------------------------------------------------        
# ------------- Find the best fit for the data
# ------------------------------------------------------------------------------------------        

        best_fit <- data.frame(type=rep(c("Quadratic", 'Linear', 'Exponential'), 3), organ=rep(c("prim", "lat", "count"), each=3), r2 = rep(0,9))
        
        if(input$fitting == "Find best"){
          n_sample <-sample(l.prim.roots, min(input$sample, length(l.prim.roots)))
          for(i in n_sample){
            temp <- root[root$root ==i,]
            if(nrow(temp) > 2){
              # QUADRATIC
                # Primary factor
                fit <- lm(sqrt(temp$prim_growth) ~ temp$dag)
                best_fit$r2[1] <- sum(best_fit$r2[1], summary(fit)$r.squared, na.rm=T)
                  
                # Lateral factor
                fit <- lm(sqrt(temp$lat_growth) ~ temp$dag)
                best_fit$r2[4] <- sum(best_fit$r2[2], summary(fit)$r.squared, na.rm=T)
                  
                # Lateral count factor
                fit <- lm(sqrt(temp$lat_number_growth) ~ temp$dag)
                best_fit$r2[7] <- sum(best_fit$r2[3], summary(fit)$r.squared, na.rm=T)
                  
                  
              ## LINEAR
                # Primary factor
                fit <- lm((temp$prim_growth) ~ temp$dag)
                best_fit$r2[2] <- sum(best_fit$r2[4], summary(fit)$r.squared, na.rm=T)
                
                # Lateral factor
                fit <- lm((temp$lat_growth) ~ temp$dag)
                best_fit$r2[5] <-sum(best_fit$r2[5], summary(fit)$r.squared, na.rm=T)
                
                # Lateral count factor
                fit <- lm((temp$lat_number_growth) ~ temp$dag)
                best_fit$r2[8] <- sum(best_fit$r2[6], summary(fit)$r.squared, na.rm=T)
                
                
              ## EXPONENTIAL
                # Primary factor                
                tryCatch({mod <- nls(prim_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
                  RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
                  TSS <- sum((temp$prim_growth - mean(temp$prim_growth))^2)  # Total sum of squares
                  best_fit$r2[3] <- sum(best_fit$r2[7], (1 - (RSS.p/TSS)), na.rm=T)  # R-squared measure
                  }, warning = function(w) {                    
                  }, error = function(e) {
                  }, finally = {
                  })      
                # Lateral factor
                tryCatch({mod <- nls(lat_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
                  RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
                  TSS <- sum((temp$lat_growth - mean(temp$lat_growth))^2)  # Total sum of squares
                  best_fit$r2[6] <- sum(best_fit$r2[8], (1 - (RSS.p/TSS)), na.rm=T)
                }, warning = function(w) {                    
                }, error = function(e) {
                }, finally = {
                })      
                # Lateral count factor
                tryCatch({mod <- nls(lat_number_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
                  RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
                  TSS <- sum((temp$lat_number_growth - mean(temp$lat_number_growth))^2)  # Total sum of squares
                  best_fit$r2[9] <- sum(best_fit$r2[9], (1 - (RSS.p/TSS)), na.rm=T)
                }, warning = function(w) {                    
                }, error = function(e) {
                }, finally = {
                })      
            }
          }
          
          fitting_results <- best_fit
          fitting_results$r2 <- fitting_results$r2 / (min(as.numeric(input$sample), length(l.prim.roots)))
          best_fit <- ddply(best_fit, .(organ), summarise, r_sqrt = max(r2), type=type[r2 == max(r2)])
          
        }else{
          best_fit <- best_fit[c(1:3),]
          best_fit$type = rep(input$fitting, nrow(best_fit))
          fitting_results <- best_fit
        }
        best_fit <<- best_fit
        
        message(paste("Best fit for primary = ", best_fit$type[best_fit$organ == "prim"]))
        message(paste("Best fit for laterals = ", best_fit$type[best_fit$organ == "lat"]))
        message(paste("Best fit for count = ", best_fit$type[best_fit$organ == "count"]))
        
        #input <- data.frame(fitting = "Find best")
        
# ------------------------------------------------------------------------------------------        
# ------------- Get the estimators
# ------------------------------------------------------------------------------------------

        for(p in l.prim.roots){   
          temp <- root[root$root == p,]
          
          # Primary factor          
          if(input$fitting == "Quadratic" ||(input$fitting == "Find best" && best_fit$type[best_fit$organ == "prim"] == "Quadratic")){
            fit <- lm(sqrt(temp$prim_growth) ~ temp$dag)
            growth.data$prim_factor[growth.data$root == p] <- fit$coefficients[2]^2
            growth.data$prim_factor2[growth.data$root == p] <- 0
            growth.data$prim_r2[growth.data$root == p] <- summary(fit)$r.squared
            
          }else if(input$fitting == "Linear" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Linear")){
            fit <- lm((temp$prim_growth) ~ temp$dag)
            growth.data$prim_factor[growth.data$root == p] <- fit$coefficients[2]
            growth.data$prim_factor2[growth.data$root == p] <- 0
            growth.data$prim_r2[growth.data$root == p] <- summary(fit)$r.squared
            
          }else if(input$fitting == "Exponential" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Exponential")){
            tryCatch({
              mod <- nls(prim_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
              growth.data$prim_factor[growth.data$root == p] <- coef(mod)[1]
              growth.data$prim_factor2[growth.data$root == p] <- coef(mod)[2]
              RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
              TSS <- sum((temp$prim_growth - mean(temp$prim_growth))^2)  # Total sum of squares
              growth.data$prim_r2[growth.data$root == p] <- (1 - (RSS.p/TSS))
            }, warning = function(w) {  
              growth.data$prim_r2[growth.data$root == p] <- 0
              growth.data$prim_factor[growth.data$root == p] <-  0
              growth.data$prim_factor2[growth.data$root == p] <-  0
            }, error = function(e) {
              growth.data$prim_r2[growth.data$root == p] <- 0
              growth.data$prim_factor[growth.data$root == p] <-  0
              growth.data$prim_factor2[growth.data$root == p] <-  0
              }, finally = {
            }) 
          }
          
          
          # Lateral factor        
          if(input$fitting == "Quadratic" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "lat"] == "Quadratic")){
            fit <- lm(sqrt(temp$lat_growth) ~ temp$dag)
            growth.data$lat_factor[growth.data$root == p] <- fit$coefficients[2]^2
            growth.data$lat_factor2[growth.data$root == p] <- 0
            growth.data$lat_r2[growth.data$root == p] <- summary(fit)$r.squared
            
          }else if(input$fitting == "Linear" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "lat"] == "Linear")){
            fit <- lm((temp$lat_growth) ~ temp$dag)
            growth.data$lat_factor[growth.data$root == p] <- fit$coefficients[2]
            growth.data$lat_factor2[growth.data$root == p] <- 0
            growth.data$lat_r2[growth.data$root == p] <- summary(fit)$r.squared
            
          }else if(input$fitting == "Exponential" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "lat"] == "Exponential")){
            tryCatch({
              mod <- nls(lat_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
              growth.data$lat_factor[growth.data$root == p] <- coef(mod)[1]
              growth.data$lat_factor2[growth.data$root == p] <- coef(mod)[2]
              RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
              TSS <- sum((temp$lat_growth - mean(temp$lat_growth))^2)  # Total sum of squares
              growth.data$lat_r2[growth.data$root == p] <- (1 - (RSS.p/TSS))            
            }, warning = function(w) {  
              growth.data$lat_factor[growth.data$root == p] <- 0
              growth.data$lat_factor2[growth.data$root == p] <- 0
              growth.data$lat_r2[growth.data$root == p] <- 0
            }, error = function(e) {
              growth.data$lat_factor[growth.data$root == p] <- 0
              growth.data$lat_factor2[growth.data$root == p] <- 0
              growth.data$lat_r2[growth.data$root == p] <- 0
            }, finally = {
            })
          }
          
           
          
          # Count factor          
          if(input$fitting == "Quadratic" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "count"] == "Quadratic")){
            fit <- lm(sqrt(temp$lat_number_growth) ~ temp$dag)
            growth.data$lat_num_factor[growth.data$root == p] <- fit$coefficients[2]^2
            growth.data$lat_num_factor2[growth.data$root == p] <- 0
            growth.data$lat_num_r2[growth.data$root == p] <- summary(fit)$r.squared
            
          }else if(input$fitting == "Linear" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "count"] == "Linear")){
            fit <- lm((temp$lat_number_growth) ~ temp$dag)
            growth.data$lat_num_factor[growth.data$root == p] <- fit$coefficients[2]
            growth.data$lat_num_factor2[growth.data$root == p] <- 0
            growth.data$lat_num_r2[growth.data$root == p] <- summary(fit)$r.squared 
            
          }else if(input$fitting == "Exponential" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "count"] == "Exponential")){
            tryCatch({
              mod <- nls(lat_number_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
              growth.data$lat_num_factor[growth.data$root == p] <- coef(mod)[1]
              growth.data$lat_num_factor2[growth.data$root == p] <- coef(mod)[2]              
              RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
              TSS <- sum((temp$lat_number_growth - mean(temp$lat_number_growth))^2)  # Total sum of squares
              growth.data$lat_num_r2[growth.data$root == p] <- (1 - (RSS.p/TSS))   
            }, warning = function(w) {  
              message(paste("Warning: ",w))
              growth.data$lat_num_r2[growth.data$root == p] <- 0
              growth.data$lat_num_factor[growth.data$root == p] <- 0
              growth.data$lat_num_factor2[growth.data$root == p] <- 0
            }, error = function(e) {
              message(paste("Error: ",e))
              growth.data$lat_num_r2[growth.data$root == p] <- 0
              growth.data$lat_num_factor[growth.data$root == p] <- 0
              growth.data$lat_num_factor2[growth.data$root == p] <- 0
            })
          }
          
        }
                
        growth.data <- growth.data[!is.na(growth.data$prim_factor),]
        growth.data <- growth.data[!is.nan(growth.data$lat_r2),]

        # Merge the factors by genotype / treatment
        factor <- ddply(growth.data, .(genotype, treatment), summarize,
                        prim = mean(prim_factor, na.rm = T),
                        prim2 = mean(prim_factor2, na.rm = T),
                        prim_r2 = mean(prim_r2, na.rm = T),
                        prim_sd = sd(prim_factor, na.rm = T),
                        lat = mean(lat_factor, na.rm = T),
                        lat2 = mean(lat_factor2, na.rm = T),
                        lat_r2 = mean(lat_r2, na.rm = T),
                        lat_sd = sd(lat_factor, na.rm = T),
                        lat_number = mean(lat_num_factor, na.rm = T),
                        lat_number2 = mean(lat_num_factor2, na.rm = T),
                        lat_number_r2 = mean(lat_num_r2, na.rm = T),
                        lat_number_sd = sd(lat_num_factor, na.rm = T))
        
        factor$prim_fit <- rep(best_fit$type[best_fit$organ == "prim"][1], nrow(factor))
        factor$lat_fit <- rep(best_fit$type[best_fit$organ == "lat"][1], nrow(factor))
        factor$count_fit <- rep(best_fit$type[best_fit$organ == "count"][1], nrow(factor))
        
# ------------------------------------------------------------------------------------------        
# ------------- Merge the data by genotype / treamtent / dag
# ------------------------------------------------------------------------------------------        

data <- ddply(root, .(genotype, treatment, dag), summarize,
                      mean_prim = mean(prim_length, na.rm = T),
                      sd_prim = sd(prim_length, na.rm = T),
                      mean_lat = mean(mean_lat_length, na.rm = T),
                      sd_lat = sd(mean_lat_length, na.rm = T),
                      mean_lat_number = mean(lat_number, na.rm = T),
                      sd_lat_number = sd(lat_number, na.rm = T),
                      mean_tot = mean(tot_length, na.rm = T),
                      sd_tot = sd(tot_length, na.rm = T)
        )
        
# ------------------------------------------------------------------------------------------        
# ------------- Compute the modelled data
# ------------------------------------------------------------------------------------------

        for(g in l.genotype){
          for(t in l.treatment){
            temp <- data[data$genotype == g & data$treatment == t & data$dag == min(data$dag),]
            
            # ------------------------------------------------------------------------------------------
            message(paste("Gen:", g, "/ TR:", t))
            # Model for the primary            
            message("Modelling primary growth")
            if(input$fitting == "Quadratic" ||(input$fitting == "Find best" && best_fit$type[best_fit$organ == "prim"] == "Quadratic")){
              data$prim_model[data$genotype == g & data$treatment == t] <- 
                temp$mean_prim + factor$prim[factor$genotype == g & factor$treatment == t] * (l.dag^2)
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Linear")){
              data$prim_model[data$genotype == g & data$treatment == t] <- 
                temp$mean_prim + factor$prim[factor$genotype == g & factor$treatment == t] * (l.dag)
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Exponential")){
              a <- factor$prim[factor$genotype == g & factor$treatment == t]
              b <- factor$prim2[factor$genotype == g & factor$treatment == t]
              data$prim_model[data$genotype == g & data$treatment == t] <- a + exp(b*(l.dag))
            }
            
            # Model for the laterals            
            if(input$fitting == "Quadratic" ||(input$fitting == "Find best" && best_fit$type[best_fit$organ == "lat"] == "Quadratic")){
              message("Modelling quadratic lateral growth")
              data$lat_model[data$genotype == g & data$treatment == t] <- 
                temp$mean_lat + factor$lat[factor$genotype == g & factor$treatment == t] * ((l.dag-min(l.dag))^2)
              
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "lat"] == "Linear")){
              message("Modelling linear lateral growth")
              data$lat_model[data$genotype == g & data$treatment == t] <- 
                temp$mean_lat + factor$lat[factor$genotype == g & factor$treatment == t] * (l.dag-min(l.dag))
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Exponential")){
              a <- factor$lat[factor$genotype == g & factor$treatment == t]
              b <- factor$lat2[factor$genotype == g & factor$treatment == t]
              message("Modelling exponential lateral growth")
              data$lat_model[data$genotype == g & data$treatment == t] <- a + exp(b*(l.dag))
            }
            
            message("Modelling lateral number")
            # Model for the lateral count           
            if(input$fitting == "Quadratic" ||(input$fitting == "Find best" && best_fit$type[best_fit$organ == "count"] == "Quadratic")){
              data$lat_number_model[data$genotype == g & data$treatment == t] <- 
                temp$mean_lat_number + factor$lat_number[factor$genotype == g & factor$treatment == t] * ((l.dag-min(l.dag))^2)  
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "count"] == "Linear")){
              data$lat_number_model[data$genotype == g & data$treatment == t] <- 
                temp$mean_lat_number + factor$lat_number[factor$genotype == g & factor$treatment == t] * (l.dag-min(l.dag))
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Exponential")){
              a <- factor$lat_number[factor$genotype == g & factor$treatment == t]
              b <- factor$lat_number2[factor$genotype == g & factor$treatment == t]
              data$lat_number_model[data$genotype == g & data$treatment == t] <- a + exp(b*(l.dag))
            }
            
            # Model for the total root length
            data$tot_model[data$genotype == g & data$treatment == t] <- 
              data$prim_model[data$genotype == g & data$treatment == t] + 
              data$lat_model[data$genotype == g & data$treatment == t] * 
              data$lat_number_model[data$genotype == g & data$treatment == t]
          }
        }
        
        
      
# Return the results

        rootfit <<- list()
        rootfit$data <<- data
        rootfit$factor <<- factor
        rootfit$growth <<- growth.data[,-1]
        rootfit$fitting_results <<- fitting_results

        message("RESULTS EXPORTED")

        return(rootfit)            
      
    })
  
    
    # Plot the different growth factors
    factorPlot <- function(){
      rs  <- Results()$growth
     # rs <- rootfit$growth
     n <- min(input$n_plot, length(unique(rs$genotype)))
     gens <- unique(rs$genotype)[1:n]
     message(gens)
     rs <- rs[rs$genotype %in% gens,]
     
     rs <<- rs
     
     if(input$plotting == "Treatment"){
        plot1 <- ggplot(rs, aes(factor(treatment), prim_factor, fill=treatment)) +
          geom_boxplot() + 
          facet_grid(genotype ~ .)+
          labs(x = "Treatment", y = "Primary growth factor [-]") +
          ggtitle(paste("Primary root growth [", best_fit$type[best_fit$organ == "prim"] ,"]")) + 
          theme_bw()   
        
        plot2 <- ggplot(rs, aes(factor(treatment), lat_factor, fill=treatment)) +
          geom_boxplot() + 
          facet_grid(genotype ~ .)+
          labs(x = "Treatment", y = "Lateral growth factor [-]") +
          ggtitle(paste("Lateral root growth [", best_fit$type[best_fit$organ == "lat"] ,"]")) + 
          theme_bw()   
        
        plot3 <- ggplot(rs, aes(factor(treatment), lat_num_factor, fill=treatment)) +
          geom_boxplot() + 
          facet_grid(genotype ~ .)+
          labs(x = "Treatment", y = "Lateral number factor [-]") +
          ggtitle(paste("Lateral root number [", best_fit$type[best_fit$organ == "count"] ,"]")) + 
          theme_bw()   
     }else{
       plot1 <- ggplot(rs, aes(factor(genotype), prim_factor, fill=genotype)) +
         geom_boxplot() + 
         facet_grid(treatment ~ .)+
         labs(x = "Genotype", y = "Primary growth factor [-]") +
         ggtitle(paste("Primary root growth [", best_fit$type[best_fit$organ == "prim"] ,"]")) + 
         theme_bw()   
       
       plot2 <- ggplot(rs, aes(factor(genotype), lat_factor, fill=genotype)) +
         geom_boxplot() + 
         facet_grid(treatment ~ .)+
         labs(x = "Genotype", y = "Lateral growth factor [-]") +
         ggtitle(paste("Lateral root growth [", best_fit$type[best_fit$organ == "lat"] ,"]")) + 
         theme_bw()   
       
       plot3 <- ggplot(rs, aes(factor(genotype), lat_num_factor, fill=genotype)) +
         geom_boxplot() + 
         facet_grid(treatment ~ .)+
         labs(x = "Genotype", y = "Lateral number factor [-]") +
         ggtitle(paste("Lateral root number [", best_fit$type[best_fit$organ == "count"] ,"]")) + 
         theme_bw()  
     }
      grid.arrange(plot1, plot2, plot3, nrow=2, ncol=2)
      
      
    }
    
    output$factorPlot <- renderPlot({
      print(factorPlot())
    }, height=600)

    output$downloadPlot2 <- downloadHandler(
      filename = "rootfit_factor.png",
      content = function(file) {
        png(file, width = 1000, height=1200)
        factorPlot()
        dev.off()
      })  
    #------------------------------------------------------
    # Plot the model results
    
    modelPlot <- function(){      
      
      if(input$runROOTFIT == 0){return()}
      
      data  <- Results()$data
      n <- min(input$n_plot, length(unique(data$genotype)))
      gens <- unique(data$genotype)[1:n]
      data <- data[data$genotype %in% gens,]
      
      pd <- position_dodge(0.5) # move them .05 to the left and right
      
      if(input$plotting == "Treatment"){
        plot1 <- ggplot(data, aes(dag, mean_prim, colour=treatment, ymax=max(mean_prim)), title="Primary root length") +
          geom_line(size=1, position=pd, ylim=range(data$mean_prim)) + 
          geom_line(aes(dag, prim_model, colour=treatment), size=1, position=pd, lty=3) + 
          facet_grid(genotype ~ .)+
          geom_errorbar(aes(ymin=(mean_prim-sd_prim), ymax=(mean_prim+sd_prim)), width=0, position=pd, size=1) +
          xlab("time [days after treatment]") +
          ylab("length [cm]") + 
          ggtitle(paste("Primary root growth [", best_fit$type[best_fit$organ == "prim"] ,"]")) + 
          theme_bw()
        
        plot2 <- ggplot(data, aes(dag, mean_lat, colour=treatment, ymax=max(mean_lat)), title="Primary root length") +
          geom_line(size=1, position=pd) + 
          geom_line(aes(dag, lat_model, colour=treatment), size=1, position=pd, lty=3) + 
          facet_grid(genotype ~ .)+
          geom_errorbar(aes(ymin=mean_lat-sd_lat, ymax=mean_lat+sd_lat), width=0, position=pd, size=1) +
          xlab("time [days after treatment]") +
          ylab("length [cm]") + 
          ggtitle(paste("Lateral root growth [", best_fit$type[best_fit$organ == "lat"] ,"]")) + 
          theme_bw()
        
        plot3 <- ggplot(data, aes(dag, mean_lat_number, colour=treatment, ymax=max(mean_lat_number)), title="lat_numberary root length") +
          geom_line(size=1, position=pd) + 
          geom_line(aes(dag, lat_number_model, colour=treatment), size=1, position=pd, lty=3) + 
          facet_grid(genotype ~ .)+
          geom_errorbar(aes(ymin=mean_lat_number-sd_lat_number, ymax=mean_lat_number+sd_lat_number), width=0, position=pd, size=1) +
          xlab("time [days after treatment]") +
          ylab("number of roots [-]") +
          ggtitle(paste("Lateral root count [", best_fit$type[best_fit$organ == "count"] ,"]")) + 
          theme_bw()
        
        plot4 <- ggplot(data, aes(dag, mean_tot, colour=treatment, ymax=max(mean_tot)), title="totary root length") +
          geom_line(size=1, position=pd) + 
          geom_line(aes(dag, tot_model, colour=treatment), size=1, position=pd, lty=3) + 
          facet_grid(genotype ~ .)+
          geom_errorbar(aes(ymin=mean_tot-sd_tot, ymax=mean_tot+sd_tot), width=0, position=pd, size=1) +
          xlab("time [days after treatment]") +
          ylab("length [cm]") + 
          ggtitle("Total root growth") + 
          theme_bw()
      }else{
        plot1 <- ggplot(data, aes(dag, mean_prim, colour=genotype, ymax=max(mean_prim))) +
          geom_line(size=1, position=pd, ylim=range(data$mean_prim)) + 
          geom_line(aes(dag, prim_model, colour=genotype), size=1, position=pd, lty=3) + 
          facet_grid(treatment ~ .)+
          geom_errorbar(aes(ymin=(mean_prim-sd_prim), ymax=(mean_prim+sd_prim)), width=0, position=pd, size=1) +
          xlab("time [days after treatment]") +
          ylab("length [cm]") + 
          ggtitle(paste("Primary root growth [", best_fit$type[best_fit$organ == "prim"] ,"]")) + 
          theme_bw()
        
        plot2 <- ggplot(data, aes(dag, mean_lat, colour=genotype, ymax=max(mean_lat))) +
          geom_line(size=1, position=pd) + 
          geom_line(aes(dag, lat_model, colour=genotype), size=1, position=pd, lty=3) + 
          facet_grid(treatment ~ .)+
          geom_errorbar(aes(ymin=mean_lat-sd_lat, ymax=mean_lat+sd_lat), width=0, position=pd, size=1) +
          xlab("time [days after treatment]") +
          ylab("length [cm]") + 
          ggtitle(paste("Lateral root growth [", best_fit$type[best_fit$organ == "lat"] ,"]")) + 
          theme_bw()
        
        plot3 <- ggplot(data, aes(dag, mean_lat_number, colour=genotype, ymax=max(mean_lat_number))) +
          geom_line(size=1, position=pd) + 
          geom_line(aes(dag, lat_number_model, colour=genotype), size=1, position=pd, lty=3) + 
          facet_grid(treatment ~ .)+
          geom_errorbar(aes(ymin=mean_lat_number-sd_lat_number, ymax=mean_lat_number+sd_lat_number), width=0, position=pd, size=1) +
          xlab("time [days after treatment]") +
          ylab("number of roots [-]") +
          ggtitle(paste("Lateral root count [", best_fit$type[best_fit$organ == "count"] ,"]")) + 
          theme_bw()
        
        plot4 <- ggplot(data, aes(dag, mean_tot, colour=genotype, ymax=max(mean_tot))) +
          geom_line(size=1, position=pd) + 
          geom_line(aes(dag, tot_model, colour=genotype), size=1, position=pd, lty=3) + 
          facet_grid(treatment ~ .)+
          geom_errorbar(aes(ymin=mean_tot-sd_tot, ymax=mean_tot+sd_tot), width=0, position=pd, size=1) +
          xlab("time [days after treatment]") +
          ylab("length [cm]") + 
          ggtitle("Total root growth") + 
          theme_bw()
      }
      grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
      
        
    }
    
    output$modelPlot <- renderPlot({
      print(modelPlot())
    }, height=600)
    output$downloadPlot <- downloadHandler(
      filename = "rootfit_fitting.png",
      content = function(file) {
        png(file, width = 1000, height=1200)
        modelPlot()
        dev.off()
    })       
    
    #------------------------------------------------------
    # Visualize the result of the matching as a table
    
    
    output$results <- renderTable({
      if (is.null(input$data_file) || input$runROOTFIT == 0) { return()}
      rs <- Results()$growth
      rs

    })


    #------------------------------------------------------
    # Visualize the result of the matching as a table
    
    
    output$fitting_results <- renderTable({
      if (is.null(input$data_file) || input$runROOTFIT == 0) { return()}
      rs <- Results()$fitting_results
      rs
    })

    #------------------------------------------------------
    # Visualize the result of the matching as a table
    
    
    output$factors <- renderTable({
      if (is.null(input$data_file) || input$runROOTFIT == 0) { return()}
      rs <- Results()$factor
      rs
      
    })
    
    #------------------------------------------------------
    # Download th result from the matching
    
    output$downloadData <- downloadHandler(
      filename = function() {"ROOT-FIT_results.csv"},
      content = function(file) {
        write.csv(Results()$growth, file)
      }
    )

  #------------------------------------------------------
  # Download th result from the matching
  
  output$downloadFactors <- downloadHandler(
    filename = function() {"ROOT-FIT_factors.csv"},
    content = function(file) {
      write.csv(Results()$factor, file)
    }
  )
    
  }
)