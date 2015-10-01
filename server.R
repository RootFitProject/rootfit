# Guillaume Lobet - University of Liege


# ROOT-FIT
# Julkowska et al. (2014). 
# Capturing Arabidopsis root architecture dynamics with ROOT-FIT reveals 
# diversity in responses to salinity. 
# Plant Physiology. doi:10.1104/pp.114.248963


shinyServer(
  
  function(input, output, clientData, session) {
   
#------------------------------------------------------
#------------------------------------------------------
# UPDATE THE UI DATA
#------------------------------------------------------
#------------------------------------------------------

     observe({
       
       rs <- Data()  
       
       #rs <- read.csv("~/Desktop/data_magda/RILS_small.csv")       
       
       gens <- na.omit(unique(factor(rs$Genotype)))
       trs <- na.omit(factor(unique(rs$Media)))
       vars <- colnames(rs)
       
       # Genotype list
       s_options <- list()
       for(g in gens) s_options[[g]] <- g
       updateSelectInput(session, "ref_gen", choices = s_options)  
       
       
       # Treatment list
       t_options <- list()
       for(t in trs) t_options[[t]] <- t
       updateSelectInput(session, "ref_cond", choices = t_options)  
       
       # Starting time
       updateNumericInput(session, "time_of_treatment", value = min(rs$Age, na.rm = T))
       
       # Maximal number of observation for the fitting
       nObs <- nrow(rs)
       updateNumericInput(session, "sample", value = round(nObs/2), min = min(10, nObs), max = nObs, step = 10)
       
       # Genotype check box
#       updateCheckboxGroupInput(session, "to_plot", choices = s_options, selected = s_options)
       updateSelectizeInput(session, 'to_plot', choices = as.character(s_options), server = TRUE)

      # updateSelectizeInput(session, 'to_analyse', choices = as.character(vars), server = TRUE)
       
    })
      
#------------------------------------------------------
#------------------------------------------------------
# LOAD THE USER DATA
#------------------------------------------------------
#------------------------------------------------------
    
    Data <- reactive({
      
      inFile <- input$data_file   
      if (is.null(inFile)) return(NULL)
      data <- read.csv(inFile$datapath) 
      
#      data <- read.csv("~/Desktop/data_magda/RILS_small.csv")
      
      cols <- c("Age", "Plate.name", "Media", "Genotype", "Plant.id", "MR.path.length", "Number.LR.MR", "Total.root.size")
      for(co in cols){
        if(length(data[[co]]) == 0){
          message(paste0(co, " column is missing in the dataset"))
        }
      }

      return(data)
      
    })   
    
#------------------------------------------------------
#------------------------------------------------------
# PROCESS THE DATA
#------------------------------------------------------
#------------------------------------------------------
    
    Results <- reactive({
      
      if(input$runROOTFIT == 0){return()}
    
      time_of_treatment <- as.numeric(input$time_of_treatment)
      
      rs <- Data()      
        #rs <- read.csv("~/Desktop/RILS_small.csv")
        
        # Remove NA's
        ind <- c()
        for(i in ncol(rs)) ind <- c(ind, which(is.na(rs[,i])))
        if(length(ind) > 0) rs <- rs[-ind,]
        
        root <- data.frame(root = paste(rs$Plate.name, "-", rs$Media, "-", rs$Genotype,"-", rs$Plant.id, sep=""))
        root$genotype <- factor(rs$Genotype)
        root$treatment <- factor(rs$Media)

        # Create lists for the different variables
        l.roots <- unique(root$root) # Get the list of roots
        l.prim.roots <- unique(root$root) # Get the list of primary roots
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
        
        root$dag <- rs$Age - time_of_treatment
        root <- root[root$dag >= 0,]
        l.dag <- unique(root$dag) # Get the list of days
        
        root <- root[!is.na(root$prim_length),]
        
        # Create a new table for each primary in order to retrieve the different factors
        growth.data <- ddply(root, .(root), summarize, 
                             genotype = unique(genotype)[1],
                             treatment = unique(treatment)[1])
        
        
# ------------------------------------------------------------------------------------------        
# ------------- Find the best fit for the data
# ------------------------------------------------------------------------------------------        

        best_fit <- data.frame(type=rep(c("Quadratic", 'Linear', 'Exponential'), 3), organ=rep(c("prim", "lat", "count"), each=3), r2 = rep(0,9))
        
        #input <- data.frame(fitting = "Linear", sample=10)
        
        if(input$fitting == "Find best"){
          n_sample <-sample(l.prim.roots, min(input$sample, length(l.prim.roots)))
          #n_sample <-sample(l.prim.roots,10)
          
          for(p in n_sample){
            temp <- root[root$root ==p,]
            if(nrow(temp) > 2){
              # QUADRATIC
                temp$dag2 <- temp$dag^2
                
                # Primary factor
                fit <- lm(temp$prim_growth ~ temp$dag2)
                best_fit$r2[1] <- sum(best_fit$r2[1], summary(fit)$r.squared, na.rm=T)
                  
                # Lateral factor
                fit <- lm(temp$lat_growth ~ temp$dag2)
                best_fit$r2[4] <- sum(best_fit$r2[4], summary(fit)$r.squared, na.rm=T)
                  
                # Lateral count factor
                fit <- lm(temp$lat_number_growth ~ temp$dag2)
                best_fit$r2[7] <- sum(best_fit$r2[7], summary(fit)$r.squared, na.rm=T)
                  
                  
              ## LINEAR
                # Primary factor
                fit <- lm((temp$prim_growth) ~ temp$dag)
                best_fit$r2[2] <- sum(best_fit$r2[2], summary(fit)$r.squared, na.rm=T)
                
                # Lateral factor
                fit <- lm((temp$lat_growth) ~ temp$dag)
                best_fit$r2[5] <-sum(best_fit$r2[5], summary(fit)$r.squared, na.rm=T)
                
                # Lateral count factor
                fit <- lm((temp$lat_number_growth) ~ temp$dag)
                best_fit$r2[8] <- sum(best_fit$r2[8], summary(fit)$r.squared, na.rm=T)
                
                
              ## EXPONENTIAL
                # Primary factor                
                tryCatch({
                  mod <- nls(prim_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
                  RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
                  TSS <- sum((temp$prim_growth - mean(temp$prim_growth))^2)  # Total sum of squares
                  best_fit$r2[3] <- sum(best_fit$r2[3], (1 - (RSS.p/TSS)), na.rm=T)  # R-squared measure
                  }, warning = function(w) {                    
                  }, error = function(e) {
                  }, finally = {
                  })      
                # Lateral factor
                tryCatch({
                  mod <- nls(lat_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
                  RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
                  TSS <- sum((temp$lat_growth - mean(temp$lat_growth))^2)  # Total sum of squares
                  best_fit$r2[6] <- sum(best_fit$r2[6], (1 - (RSS.p/TSS)), na.rm=T)
                }, warning = function(w) {                    
                }, error = function(e) {
                }, finally = {
                })      
                # Lateral count factor
                tryCatch({
                  mod <- nls(lat_number_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
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
          best_fit <- best_fit[c(1,4,7),]
          best_fit$type = rep(input$fitting, nrow(best_fit))
          fitting_results <- best_fit
        }
        best_fit <<- best_fit
        
        message(paste("Best fit for primary = ", best_fit$type[best_fit$organ == "prim"]))
        message(paste("Best fit for laterals = ", best_fit$type[best_fit$organ == "lat"]))
        message(paste("Best fit for count = ", best_fit$type[best_fit$organ == "count"]))
        

# ------------------------------------------------------------------------------------------        
# ------------- Get the estimators
# ------------------------------------------------------------------------------------------


        for(p in l.prim.roots){   
          temp <- root[root$root == p,]
          if(nrow(temp) > 2){
            # Primary factor          
            if(input$fitting == "Quadratic" ||(input$fitting == "Find best" && best_fit$type[best_fit$organ == "prim"] == "Quadratic")){
              fit <-lm(sqrt(temp$prim_growth) ~  temp$dag)
              growth.data$prim_a[growth.data$root == p] <- fit$coefficients[1]
              growth.data$prim_b[growth.data$root == p] <- fit$coefficients[2]
              growth.data$prim_r2[growth.data$root == p] <- summary(fit)$r.squared
              
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Linear")){
              fit <- lm((temp$prim_growth) ~ temp$dag)
              growth.data$prim_a[growth.data$root == p] <- fit$coefficients[1]
              growth.data$prim_b[growth.data$root == p] <- fit$coefficients[2]
              growth.data$prim_r2[growth.data$root == p] <- summary(fit)$r.squared
              
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Exponential")){
              tryCatch({
                mod <- nls(prim_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
                growth.data$prim_a[growth.data$root == p] <- coef(mod)[1]
                growth.data$prim_b[growth.data$root == p] <- coef(mod)[2]
                RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
                TSS <- sum((temp$prim_growth - mean(temp$prim_growth))^2)  # Total sum of squares
                growth.data$prim_r2[growth.data$root == p] <- (1 - (RSS.p/TSS))
              }, warning = function(w) {  
                growth.data$prim_r2[growth.data$root == p] <- -1
                growth.data$prim_a[growth.data$root == p] <-  -1
                growth.data$prim_b[growth.data$root == p] <-  -1
              }, error = function(e) {
                growth.data$prim_r2[growth.data$root == p] <- -1
                growth.data$prim_a[growth.data$root == p] <- -1
                growth.data$prim_b[growth.data$root == p] <- -1
                }, finally = {
              }) 
            }
            
            
            # Lateral factor        
            if(input$fitting == "Quadratic" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "lat"] == "Quadratic")){
              fit <- lm(sqrt(temp$lat_growth)  ~ temp$dag)           
              growth.data$lat_a[growth.data$root == p] <- fit$coefficients[1]
              growth.data$lat_b[growth.data$root == p] <- fit$coefficients[2]
              growth.data$lat_r2[growth.data$root == p] <- summary(fit)$r.squared
              
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "lat"] == "Linear")){
              fit <- lm(temp$lat_growth ~ temp$dag^2)
              growth.data$lat_a[growth.data$root == p] <- fit$coefficients[1]
              growth.data$lat_b[growth.data$root == p] <- fit$coefficients[2]
              growth.data$lat_r2[growth.data$root == p] <- summary(fit)$r.squared
              
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "lat"] == "Exponential")){
              tryCatch({
                mod <- nls(lat_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))                
                growth.data$lat_a[growth.data$root == p] <- coef(mod)[1]
                growth.data$lat_b[growth.data$root == p] <- coef(mod)[2]
                RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
                TSS <- sum((temp$lat_growth - mean(temp$lat_growth))^2)  # Total sum of squares
                growth.data$lat_r2[growth.data$root == p] <- (1 - (RSS.p/TSS))            
              }, warning = function(w) {  
                growth.data$lat_a[growth.data$root == p] <- -1
                growth.data$lat_b[growth.data$root == p] <- -1
                growth.data$lat_r2[growth.data$root == p] <- -1
              }, error = function(e) {
                growth.data$lat_a[growth.data$root == p] <- -1
                growth.data$lat_b[growth.data$root == p] <- -1
                growth.data$lat_r2[growth.data$root == p] <- -1
              }, finally = {
              })
            }
            
             
            
            # Count factor          
            if(input$fitting == "Quadratic" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "count"] == "Quadratic")){
              fit <- lm(sqrt(temp$lat_number_growth) ~ temp$dag)
              growth.data$lat_num_a[growth.data$root == p] <- fit$coefficients[1]
              growth.data$lat_num_b[growth.data$root == p] <- fit$coefficients[2]
              growth.data$lat_num_r2[growth.data$root == p] <- summary(fit)$r.squared
              
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "count"] == "Linear")){
              fit <- lm((temp$lat_number_growth) ~ temp$dag)
              growth.data$lat_num_a[growth.data$root == p] <- fit$coefficients[1]
              growth.data$lat_num_b[growth.data$root == p] <- fit$coefficients[2]
              growth.data$lat_num_r2[growth.data$root == p] <- summary(fit)$r.squared 
              
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" && best_fit$type[best_fit$organ == "count"] == "Exponential")){
              tryCatch({
                mod <- nls(lat_number_growth ~ a + exp( b * dag), data = temp, start = list(a = 0, b = 0))
                growth.data$lat_num_a[growth.data$root == p] <- coef(mod)[1]
                growth.data$lat_num_b[growth.data$root == p] <- coef(mod)[2]              
                RSS.p <- sum(residuals(mod)^2)  # Residual sum of squares
                TSS <- sum((temp$lat_number_growth - mean(temp$lat_number_growth))^2)  # Total sum of squares
                growth.data$lat_num_r2[growth.data$root == p] <- (1 - (RSS.p/TSS))   
              }, warning = function(w) {  
                message(paste("Warning: ",w))
                growth.data$lat_num_r2[growth.data$root == p] <- 0
                growth.data$lat_num_a[growth.data$root == p] <- 0
                growth.data$lat_num_b[growth.data$root == p] <- 0
              }, error = function(e) {
                message(paste("Error: ",e))
                growth.data$lat_num_r2[growth.data$root == p] <- 0
                growth.data$lat_num_a[growth.data$root == p] <- 0
                growth.data$lat_num_b[growth.data$root == p] <- 0
              })
            }
          }
          
        }
                
        growth.data <- growth.data[!is.na(growth.data$prim_a),]
        growth.data <- growth.data[!is.nan(growth.data$lat_r2),]
        growth.data <<- growth.data
        
        ## Comput the relative growth data
        relative.growth.data <- growth.data
        genotype.growth.data <- growth.data
        condition.growth.data <- growth.data
        for(i in 4:ncol(growth.data)){
          mean <- mean(growth.data[,i][growth.data$genotype == input$ref_gen & growth.data$treatment == input$ref_cond])
          relative.growth.data[,i] <- growth.data[,i] / mean 
          
          for(g in l.treatment){
            temp <- growth.data[growth.data$treatment == g,]
            mean <- mean(temp[,i][temp$genotype == input$ref_gen])
            genotype.growth.data[,i][genotype.growth.data$treatment == g] <- temp[,i] / mean 
          }
          
          for(g in l.genotype){
            temp <- growth.data[growth.data$genotype == g,]
            mean <- mean(temp[,i][temp$treatment == input$ref_cond])
            condition.growth.data[,i][condition.growth.data$genotype == g] <- temp[,i] / mean 
          }
        }
        
        
        

        # Merge the factors by genotype / treatment
        factor <- ddply(growth.data, .(genotype, treatment), summarize,
                        prim_a = mean(prim_a, na.rm = T),
                        prim_b = mean(prim_b, na.rm = T),
                        prim_r2 = mean(prim_r2, na.rm = T),
                        lat_a = mean(lat_a, na.rm = T),
                        lat_b = mean(lat_b, na.rm = T),
                        lat_r2 = mean(lat_r2, na.rm = T),
                        lat_number_a = mean(lat_num_a, na.rm = T),
                        lat_number_b = mean(lat_num_b, na.rm = T),
                        lat_number_r2 = mean(lat_num_r2, na.rm = T)
                        )
        
        factor$prim_fit <- rep(best_fit$type[best_fit$organ == "prim"][1], nrow(factor))
        factor$lat_fit <- rep(best_fit$type[best_fit$organ == "lat"][1], nrow(factor))
        factor$count_fit <- rep(best_fit$type[best_fit$organ == "count"][1], nrow(factor))
        
# ------------------------------------------------------------------------------------------        
# ------------- Merge the data by genotype / treamtent / dag
# ------------------------------------------------------------------------------------------        

      data <- ddply(root, .(genotype, treatment, dag), summarize,
                      mean_prim = mean(prim_growth, na.rm = T),
                      sd_prim = sd(prim_growth, na.rm = T),
                      mean_lat = mean(lat_growth, na.rm = T),
                      sd_lat = sd(lat_growth, na.rm = T),
                      mean_lat_number = mean(lat_number_growth, na.rm = T),
                      sd_lat_number = sd(lat_number_growth, na.rm = T),
                      mean_tot = mean(tot_growth, na.rm = T),
                      sd_tot = sd(tot_growth, na.rm = T)
        )
        
# ------------------------------------------------------------------------------------------        
# ------------- Compute the modelled data
# ------------------------------------------------------------------------------------------

        for(g in l.genotype){
          for(t in l.treatment){
            temp <- data[data$genotype == g & data$treatment == t & data$dag == min(data$dag),]
            l.dag <- data$dag[data$genotype == g & data$treatment == t] 
              
            # ------------------------------------------------------------------------------------------
            # Model for the primary            
            if(input$fitting == "Quadratic" ||(input$fitting == "Find best" && best_fit$type[best_fit$organ == "prim"] == "Quadratic")){
              data$prim_model[data$genotype == g & data$treatment == t] <- 
                factor$prim_a[factor$genotype == g & factor$treatment == t] + (factor$prim_b[factor$genotype == g & factor$treatment == t] * l.dag)^2
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Linear")){
              data$prim_model[data$genotype == g & data$treatment == t] <- 
                factor$prim_a[factor$genotype == g & factor$treatment == t] + factor$prim_b[factor$genotype == g & factor$treatment == t] * (l.dag)
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "prim"] == "Exponential")){
              a <- factor$prim_a[factor$genotype == g & factor$treatment == t]
              b <- factor$prim_b[factor$genotype == g & factor$treatment == t]
              data$prim_model[data$genotype == g & data$treatment == t] <- a + exp(b*(l.dag))
            }
            
            # Model for the laterals            
            if(input$fitting == "Quadratic" ||(input$fitting == "Find best" && best_fit$type[best_fit$organ == "lat"] == "Quadratic")){
              data$lat_model[data$genotype == g & data$treatment == t] <- 
                factor$lat_a[factor$genotype == g & factor$treatment == t] + (factor$lat_b[factor$genotype == g & factor$treatment == t] * l.dag)^2 
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "lat"] == "Linear")){
              data$lat_model[data$genotype == g & data$treatment == t] <- 
                factor$lat_a[factor$genotype == g & factor$treatment == t] + factor$lat_b[factor$genotype == g & factor$treatment == t] * (l.dag)
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "lat"] == "Exponential")){
              a <- factor$lat_a[factor$genotype == g & factor$treatment == t]
              b <- factor$lat_b[factor$genotype == g & factor$treatment == t]
              data$lat_model[data$genotype == g & data$treatment == t] <- a + exp(b*(l.dag))
            }
            
            # Model for the lateral count           
            if(input$fitting == "Quadratic" ||(input$fitting == "Find best" && best_fit$type[best_fit$organ == "count"] == "Quadratic")){
              data$lat_number_model[data$genotype == g & data$treatment == t] <- 
                factor$lat_number_a[factor$genotype == g & factor$treatment == t] + (factor$lat_number_b[factor$genotype == g & factor$treatment == t] * l.dag)^2  
            }else if(input$fitting == "Linear" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "count"] == "Linear")){
              data$lat_number_model[data$genotype == g & data$treatment == t] <- 
                factor$lat_number_a[factor$genotype == g & factor$treatment == t] + factor$lat_number_b[factor$genotype == g & factor$treatment == t] * (l.dag)
            }else if(input$fitting == "Exponential" || (input$fitting == "Find best" &&  best_fit$type[best_fit$organ == "count"] == "Exponential")){
              a <- factor$lat_number_a[factor$genotype == g & factor$treatment == t]
              b <- factor$lat_number_b[factor$genotype == g & factor$treatment == t]
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
        rootfit$relative_results <<- relative.growth.data
        rootfit$genotype_results <<- genotype.growth.data
        rootfit$condition_results <<- condition.growth.data
        
        return(rootfit)            
      
    })
  
    
#------------------------------------------------------
#------------------------------------------------------
# PLOTS
#------------------------------------------------------
#------------------------------------------------------
    
    # Plot the different growth factors
    factorPlot <- function(){
      
      get_letters <- function(HSD, flev, d, var){
        # Extract labels and factor levels from Tukey post-hoc 
        Tukey.levels <- HSD[[flev]][,4]
        if(length(Tukey.levels) > 1){
          Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
        }else{
          if(Tukey.levels < 0.05) diff <- c("a", "b")
          else diff <- c("a", "b")
          names(diff) <- unique(d[[flev]])
          Tukey.labels <- list(Letters = diff)
        }
        treatment <- names(Tukey.labels[['Letters']])
        
        # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
        # upper quantile and label placement
        boxplot.df <- ddply(d, flev, function (x) max(fivenum(x[[var]])) + 0.2)
        
        # Create a data frame out of the factor levels and Tukey's homogenous group letters
        plot.levels <- data.frame(treatment, labels = Tukey.labels[['Letters']],stringsAsFactors = FALSE)
        
        # Merge it with the labels
        labels.df <- merge(plot.levels, boxplot.df, by.x = 'treatment', by.y = flev, sort = FALSE)
        
        return(labels.df)
      }       
      
      
       rs  <- Results()$growth
       #rs <- read.csv("~/Desktop/ROOT-FIT_factors.csv")
       
       gens <- input$to_plot
       rs$treatment <- factor(rs$treatment)
       treats <- unique(rs$treatment)
       rs <- rs[rs$genotype %in% gens,]
       rs$genotype <- factor(rs$genotype)
       
       if(input$plotting == "Treatment"){
         
         # Get the signification letters (anova)
           labels <- data.frame(treatment = numeric(), labels = character(), V1 = numeric(), genotype = factor())
           for(g in gens){
              temp <- rs[rs$genotype == g,]
              amod <- aov(prim_b ~ treatment, data = temp)
              HSD <- TukeyHSD(amod, ordered = FALSE, conf.level = 0.95)
              tempLab <- get_letters(HSD, 'treatment', temp, 'prim_b')
              tempLab$genotype <- g
              labels <- rbind(labels, tempLab)
           }
          plot1 <- ggplot(rs, aes(factor(treatment), prim_b, fill=treatment)) +
            geom_boxplot() + 
            geom_text(data = labels, aes(x = treatment, y = V1, label = labels)) +
            facet_grid(genotype ~ .) +
            labs(x = "Treatment", y = "Primary growth factor [-]") +
            ggtitle(paste("Primary root growth [", best_fit$type[best_fit$organ == "prim"] ,"]")) + 
            theme_bw()
          
          # Get the signification letters (anova)
          labels <- data.frame(treatment = numeric(), labels = character(), V1 = numeric(), genotype = factor())
          for(g in gens){
            temp <- rs[rs$genotype == g,]
            amod <- aov(lat_b ~ treatment, data = temp)
            HSD <- TukeyHSD(amod, ordered = FALSE, conf.level = 0.95)
            tempLab <- get_letters(HSD, 'treatment', temp, 'lat_b')
            tempLab$genotype <- g
            labels <- rbind(labels, tempLab)
          }          
          plot2 <- ggplot(rs, aes(factor(treatment), lat_b, fill=treatment)) +
            geom_boxplot() + 
            geom_text(data = labels, aes(x = treatment, y = V1, label = labels)) +
            facet_grid(genotype ~ .)+
            labs(x = "Treatment", y = "Lateral growth factor [-]") +
            ggtitle(paste("Lateral root growth [", best_fit$type[best_fit$organ == "lat"] ,"]")) + 
            theme_bw() 
          
          # Get the signification letters (anova)
          labels <- data.frame(treatment = numeric(), labels = character(), V1 = numeric(), genotype = factor())
          for(g in gens){
            temp <- rs[rs$genotype == g,]
            amod <- aov(lat_num_b ~ treatment, data = temp)
            HSD <- TukeyHSD(amod, ordered = FALSE, conf.level = 0.95)
            tempLab <- get_letters(HSD, 'treatment', temp, 'lat_num_b')
            tempLab$genotype <- g
            labels <- rbind(labels, tempLab)
          }           
          plot3 <- ggplot(rs, aes(factor(treatment), lat_num_b, fill=treatment)) +
            geom_boxplot() + 
            geom_text(data = labels, aes(x = treatment, y = V1, label = labels)) +
            facet_grid(genotype ~ .)+
            labs(x = "Treatment", y = "Lateral number factor [-]") +
            ggtitle(paste("Lateral root number [", best_fit$type[best_fit$organ == "count"] ,"]")) + 
            theme_bw()
       }else{
         
         labels <- data.frame(genotype = numeric(), labels = character(), V1 = numeric(), treatment = factor())
         for(t in treats){
           temp <- rs[rs$treatment == t,]
           amod <- aov(prim_b ~ genotype, data = temp)
           HSD <- TukeyHSD(amod, ordered = FALSE, conf.level = 0.95)
           tempLab <- get_letters(HSD, 'genotype', temp, 'prim_b')
           colnames(tempLab)[1] <- "genotype"
           tempLab$treatment <- t
           labels <- rbind(labels, tempLab)
         }          
         plot1 <- ggplot(rs, aes(genotype, prim_b, fill=genotype)) +
           geom_boxplot() + 
           geom_text(data = labels, aes(x = genotype, y = V1, label = labels)) +
           facet_grid(treatment ~ .)+
           labs(x = "Genotype", y = "Primary growth factor [-]") +
           ggtitle(paste("Primary root growth [", best_fit$type[best_fit$organ == "prim"] ,"]")) + 
           theme_bw()   
         
         labels <- data.frame(genotype = numeric(), labels = character(), V1 = numeric(), treatment = factor())
         for(t in treats){
           temp <- rs[rs$treatment == t,]
           amod <- aov(lat_b ~ genotype, data = temp)
           HSD <- TukeyHSD(amod, ordered = FALSE, conf.level = 0.95)
           tempLab <- get_letters(HSD, 'genotype', temp, 'lat_b')
           colnames(tempLab)[1] <- "genotype"
           tempLab$treatment <- t
           labels <- rbind(labels, tempLab)
         }          
         plot2 <- ggplot(rs, aes(factor(genotype), lat_b, fill=genotype)) +
           geom_boxplot() + 
           geom_text(data = labels, aes(x = genotype, y = V1, label = labels)) +
           facet_grid(treatment ~ .)+
           labs(x = "Genotype", y = "Lateral growth factor [-]") +
           ggtitle(paste("Lateral root growth [", best_fit$type[best_fit$organ == "lat"] ,"]")) + 
           theme_bw()   
         
         
         labels <- data.frame(genotype = numeric(), labels = character(), V1 = numeric(), treatment = factor())
         for(t in treats){
           temp <- rs[rs$treatment == t,]
           amod <- aov(lat_num_b ~ genotype, data = temp)
           HSD <- TukeyHSD(amod, ordered = FALSE, conf.level = 0.95)
           tempLab <- get_letters(HSD, 'genotype', temp, 'lat_num_b')
           colnames(tempLab)[1] <- "genotype"
           tempLab$treatment <- t
           labels <- rbind(labels, tempLab)
         }          
         plot3 <- ggplot(rs, aes(factor(genotype), lat_num_b, fill=genotype)) +
           geom_boxplot() + 
           geom_text(data = labels, aes(x = genotype, y = V1, label = labels)) +
           facet_grid(treatment ~ .)+
           labs(x = "Genotype", y = "Lateral number factor [-]") +
           ggtitle(paste("Lateral root number [", best_fit$type[best_fit$organ == "count"] ,"]")) + 
           theme_bw()  
       }
      grid.arrange(plot1, plot2, plot3, nrow=3, ncol=1)
      
      
    }
 
    output$factorPlot <- renderPlot({
      print(factorPlot())
    }, height = 1000)
    
    output$downloadPlot2 <- downloadHandler(
      filename = "rootfit_factors.png",
      content = function(file) {
        png(file, width = input$plot_width, height=input$plot_height)
        factorPlot()
        dev.off()
      })  
    
    
    
    

    # Plot the different growth factors
    relativeFactorPlot <- function(){
      if(input$plotting == "Treatment"){
        
        rs  <- Results()$condition_results
        gens <- input$to_plot
        rs <- rs[rs$genotype %in% gens,]
        rs$treatment <- factor(rs$treatment)
        treats <- factor(sort(unique(rs$treatment)))
        
        rs1 <- data.frame(genotype = factor(rs$genotype), treatment=factor(rs$treatment), val = rs$prim_b, fact = factor(rep("primary",nrow(rs))))
        rs1 <- rbind(rs1, data.frame(genotype = factor(rs$genotype),treatment=factor(rs$treatment), val = rs$lat_b, fact = factor(rep("lateral",nrow(rs)))))
        rs1 <- rbind(rs1, data.frame(genotype = factor(rs$genotype),treatment=factor(rs$treatment), val = rs$lat_num_b, fact = factor(rep("number",nrow(rs)))))
        
        # Get the signification letters (anova)
        rs1$letters <- ""
        rs1$xcoord <- 0
        for(g in gens){
          i <- 0
          for(t in treats){
            i <- i +1
            temp <- rs1[rs1$genotype == g & rs1$treatment == t,]
            amod <- aov(val ~ fact, data = temp)
            tuk <- glht(amod, linfct = mcp(fact = "Tukey"))
            tuk.cld <- cld(tuk)$mcletters$Letters
            
            rs1$letters[rs1$genotype == g & rs1$treatment == t & rs1$fact == "primary"] <- tuk.cld["primary"]
            rs1$xcoord[rs1$genotype == g & rs1$treatment == t & rs1$fact == "primary"] <- i - 0.3
            rs1$letters[rs1$genotype == g & rs1$treatment == t & rs1$fact == "lateral"] <- tuk.cld["lateral"]
            rs1$xcoord[rs1$genotype == g & rs1$treatment == t & rs1$fact == "lateral"] <- i
            rs1$letters[rs1$genotype == g & rs1$treatment == t & rs1$fact == "number"] <- tuk.cld["number"]
            rs1$xcoord[rs1$genotype == g & rs1$treatment == t & rs1$fact == "number"] <- i + 0.3
          }
        }
        
        lines <- c(1:(length(treats)-1))+0.5
        
        plot1 <-  ggplot(rs1, aes(x=factor(treatment), y=val, fill=fact)) +
          geom_boxplot(width = 0.8) + 
          geom_text(aes(x=xcoord, y=max(val), label=letters, group=fact)) + 
          facet_grid(genotype ~ .) +
          labs(x = "Treatment", y = "Root growth relatvie factors [-]") +
          ggtitle(paste("Root growth relatvie factors [", best_fit$type[best_fit$organ == "prim"] ,"]")) + 
          geom_vline(xintercept = lines) + 
          theme_bw()   
        
      }else{

        rs  <- Results()$genotype_results
        gens <- input$to_plot
        rs <- rs[rs$genotype %in% gens,]
        rs$treatment <- factor(rs$treatment)
        treats <- factor(sort(unique(rs$treatment)))
        
        rs1 <- data.frame(genotype = factor(rs$genotype), treatment=factor(rs$treatment), val = rs$prim_b, fact = factor(rep("primary",nrow(rs))))
        rs1 <- rbind(rs1, data.frame(genotype = factor(rs$genotype),treatment=factor(rs$treatment), val = rs$lat_b, fact = factor(rep("lateral",nrow(rs)))))
        rs1 <- rbind(rs1, data.frame(genotype = factor(rs$genotype),treatment=factor(rs$treatment), val = rs$lat_num_b, fact = factor(rep("number",nrow(rs)))))
        
        # Get the signification letters (anova)
        rs1$letters <- ""
        rs1$xcoord <- 0
        for(t in treats){
          i <- 0
          for(g in gens){
            i <- i +1
            temp <- rs1[rs1$genotype == g & rs1$treatment == t,]
            amod <- aov(val ~ fact, data = temp)
            tuk <- glht(amod, linfct = mcp(fact = "Tukey"))
            tuk.cld <- cld(tuk)$mcletters$Letters
            
            rs1$letters[rs1$genotype == g & rs1$treatment == t & rs1$fact == "primary"] <- tuk.cld["primary"]
            rs1$xcoord[rs1$genotype == g & rs1$treatment == t & rs1$fact == "primary"] <- i - 0.3
            rs1$letters[rs1$genotype == g & rs1$treatment == t & rs1$fact == "lateral"] <- tuk.cld["lateral"]
            rs1$xcoord[rs1$genotype == g & rs1$treatment == t & rs1$fact == "lateral"] <- i
            rs1$letters[rs1$genotype == g & rs1$treatment == t & rs1$fact == "number"] <- tuk.cld["number"]
            rs1$xcoord[rs1$genotype == g & rs1$treatment == t & rs1$fact == "number"] <- i + 0.3
          }
        }
        
        lines <- c(1:(length(gens)-1))+0.5        
        
        plot1 <-  ggplot(rs1, aes(x=factor(genotype), y=val, fill=fact)) +
          geom_boxplot(width = 0.8) + 
          geom_text(aes(x=xcoord, y=max(val), label=letters, group=fact)) + 
          facet_grid(treatment ~ .) +
          labs(x = "Genotype", y = "Root growth relative factors [-]") +
          ggtitle(paste("Root growth relatvie factors [", best_fit$type[best_fit$organ == "prim"] ,"]")) + 
          geom_vline(xintercept = lines) + 
          theme_bw()         
      }
      plot1
    }    

    output$relativeFactorPlot <- renderPlot({
      print(relativeFactorPlot())
    })

    output$downloadPlot3 <- downloadHandler(
      filename = "rootfit_relative_factors.png",
      content = function(file) {
        png(file, width = input$plot_width, height=input$plot_height)
        relativeFactorPlot()
        dev.off()
      })  
    
    
    
    
    modelPlot <- function(){      
      
      if(input$runROOTFIT == 0){return()}
      
      data  <- Results()$data
      gens <- as.vector(input$to_plot)
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
          theme_bw() + 
          theme(text = element_text(size=input$plot_font))
        
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
      grid.arrange(plot1, plot2, plot3, plot4, nrow=4, ncol=1)
      
        
    }
    
    output$modelPlot <- renderPlot({
      print(modelPlot())
    }, height = 1000)
    
    output$downloadPlot <- downloadHandler(
      filename = "rootfit_fitting.png",
      content = function(file) {
        png(file, width = input$plot_width, height=input$plot_height)
        modelPlot()
        dev.off()
    })       
    
    
    
    
#------------------------------------------------------
#------------------------------------------------------
# TABLES
#------------------------------------------------------
#------------------------------------------------------ 
    
    output$results <- renderTable({
      if (is.null(input$data_file) || input$runROOTFIT == 0) { return()}
      rs <- Results()$data
      rs
    })
    
    output$fitting_results <- renderTable({
      if (is.null(input$data_file) || input$runROOTFIT == 0) { return()}
      rs <- Results()$fitting_results 
      rs
    })

    output$factors <- renderTable({
      if (is.null(input$data_file) || input$runROOTFIT == 0) { return()}
      rs <- Results()$growth
      rs
    })

    output$relfactorsgen <- renderTable({
      if (is.null(input$data_file) || input$runROOTFIT == 0) { return()}
      rs <- head(Results()$genotype_results)
      rs
    })
    
    output$relfactorstr <- renderTable({
      if (is.null(input$data_file) || input$runROOTFIT == 0) { return()}
      rs <- head(Results()$condition_results)
      rs
    })    
    
#------------------------------------------------------
#------------------------------------------------------
# DOWNLOAD TABLES
#------------------------------------------------------
#------------------------------------------------------ 
    
    output$downloadData <- downloadHandler(
      filename = function() {"ROOT-FIT_results.csv"},
      content = function(file) {
        write.csv(Results()$data, file)
      }
    )

    output$downloadFactors <- downloadHandler(
      filename = function() {"ROOT-FIT_factors.csv"},
      content = function(file) {
        write.csv(Results()$growth, file)
      }
    )
    
    output$downloadFactorsRelGen <- downloadHandler(
      filename = function() {"ROOT-FIT_factors_relative_genotype.csv"},
      content = function(file) {
        write.csv(Results()$genotype_results, file)
      }
    )
      
    output$downloadFactorsRelTr <- downloadHandler(
      filename = function() {"ROOT-FIT_factors_relative_treatment.csv"},
      content = function(file) {
        write.csv(Results()$condition_results, file)
      }
    )    
    
  }
)