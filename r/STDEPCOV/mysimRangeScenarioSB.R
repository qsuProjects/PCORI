mysim <- function(seed.init, comb, n_subjects, gm, lim, w, md, rv, output_location, source_location, simrange, trad_cens, admin_cens, ns, mediangoal, rho) {
  ###
  #lim: limits
  #w: shape parameter (nu) for weibull defining the g function
  #md: median survival desired - together with nu it fixes the scale parameter lambda
  #rv: covariates type: 2 normal or 1 normal 1 binary
  #comb: unique iteration number for all combinations between evaluated
  ###
  #   print(paste0(seed.init, comb, n_subjects, gm, lim, w, md, rv, output_location, source_location, range(simrange), trad_cens[1,], admin_cens[1,], ns, mediangoal, collapse = "   \n")) 
  #   print(seed.init+comb^2*100000)
  #do simulate using simrange
  #   print(match.call(expand.dots = F))
  #   print(range(simrange))
  #   print(simrange)
  #   results_time <- system.time({
  l_ply(simrange, .parallel = T, .fun = function(.df, seed.init, n_subjects, gm, lim, w, md, rv, comb, source_location, output_location, trad_cens, admin_cens, ns, mediangoal, rho) {
    DEBUG <- F
    #       .cornum <- WORKER.ID
    if (DEBUG) {  
      #       .df <- 1
      print(.df)
      WORKER.ID <- .df
      print(paste0(output_location, "/log/", WORKER.ID, "_output.txt"))
      source_location <- "~/shexport/PCORI/lib/"
      source_location <- "~/PCORI/lib/"
      n_subjects <- 100
      rv <- 1
      gm <- 1
      lim <- 1
      w <- 1
      md <- 1
      ns <- 1
      seed.init<-switch(lim,3182014,264064,856492)
      comb <- 111111
      mediangoal <- c(35,75,150)[md]
      rho <- -0.8
      
    }
    
    dataset_number <- sprintf("%04d", .df)
    cat(paste0(WORKER.ID, "-", dataset_number, " rep start at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
    ptm0 <- proc.time()
    require(msm)
    require(survival)
    require(mvtnorm)
    require(plyr)
    require(MBESS)
    require(data.table)
    source(paste0(source_location, "/jointly_generate_binary_normal.R"))
    source(paste0(source_location, "/get_censoring_times.R"))
    
    #it indicates teh unique iteration number for each combination of values we are testing
    myseed <- (seed.init+(comb^2)*100000+.df) %% .Machine$integer.max
    set.seed(myseed)
    cat(paste0(WORKER.ID, "-", dataset_number, " seed is: ", myseed, " at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
    #     rho <- c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)[ns]
    
    ### beta coefficients for generating covariates and mean
    if (rv==1){
      #2 Normals
      #true coefficients
      b1 <- 0.02
      b2 <- 0.04
      #mean/sd var1
      mu1 <- 50
      sd1 <- 10
      #mean/sd var2
      mu2 <- 30
      sd2 <- 5
    } else {
      ## z1~N(50,100) and Z2~Bin(1,0.7)
      #1 Normal + 1 Binary
      b1 <- 0.02
      b2 <- -0.5
      #mean/sd var1
      mu1 <- 50
      sd1 <- 10
      #mean/sd var2 (binary only mean necessary = proportion)
      mu2 <- 0.5
    }    
    
    #nu and lambda are the parameters of a weibull distribution
    #I fixed nu (the shape parameter <1 decreasing hazard, =1 constant hazard, >1 increasing hazard)
    #the scale parameter lambda was computed by inverting the formula for the estimated median 
    ## survival time from a weibull with shape nu and estimated median 50 where we don't take 
    ## into account the covariates, i.e., assume zb=1
    nu <- switch(w,2,1,0.5)
    
    ### lambda is different depending on whether we simulate 2 normals or 1 normal and 1 binary
    ## with 2 Normals we need to use the formula where ZB is incorporated but not for
    ## 1 normal + 1 binary
    if (rv==1){
      rv.type <- "2Normals"
      ZBbar <- mu1*b1+mu2*b2
      #ZBbar <- 0
    } else {
      rv.type <- "1Nor1Bin"
      ZBbar <- mu1*b1+mu2*b2
      #       ZBbar <- 0
    }
    
    
    ### in order to get an estimated survival time closer to what we want we need to manipulate
    # the formula and use a slightly different median. The results dependens on what the shape parameter 
    # is as well as the limits set. No matter what we do, we can't have median survival to be 100 if
    # the limit only goes up to 50!
    median <- mediangoal
    lambda <- (log(2)/exp(ZBbar))*median^(-nu)
    
    # the g function is defined as the inverse of the baseline cummulative hazard from
    ## a Weibull with shape nu and scale lambda defined above
    g <- function(x){
      ((1/lambda)*x)^(1/nu)
    }
    g.inv <- function(x){
      lambda*(x^nu)
    }  
    
    # CREATING THE TIME SCALE AND TRANSFORMED TIME SCALE
    t <- 0:350
    t.diff <- (t[-1] - t[1:(length(t) - 1)])[-(length(t) - 1)]
    g.inv.t <- g.inv(t)
    g.inv.t.diff <- (g.inv(t[-1]) - g.inv(t[1:(length(t) - 1)]))[-(length(t) - 1)]
    
    #CREATING THE BOUNDS OF TRUNCATION
    t.max <- switch(lim,300,150,50)
    t.min <- 20
    
    g.inv.t.max <- g.inv(t.max)
    g.inv.t.min <- g.inv(t.min)
    
    #DATA GENERATING PROCESS FOR COVARIATE
    
    #CREATING DATA VECTOR
    ## assumed parameters for z1 and z2 respectively
    z.list <- list()
    if (rv==1){
      # using correlation=0
      cov <- rho*sd1*sd2
      Sigma <- matrix(c(sd1^2,cov,cov,sd2^2), nrow = 2)
      for (i in 1:n_subjects) {
        vars <- rmvnorm(length(t), mean = c(mu1,mu2), sigma=Sigma)
        z1 <- vars[ ,1]
        z2 <- vars[ ,2]
        XB <- exp(b1*z1 + b2*z2)
        z.list[[i]] <- data.frame(id = i, z1 = z1,z2 =  z2, xB = XB)
      }
    } else {
      # z1~N(50,100) and Z2~Bin(1,0.7)
      # 1 Normal + 1 Binary
      
      #determine whether specified correlations are within bounds
      p = mu2
      q = 1-p
      
      #compute bound
      # using floor and * 1000 avoids "not positive definite" warning
      hiBound = floor(dnorm( qnorm(p) ) / sqrt(p*q) * 1000)
      bounded_corr <- sign(rho) * hiBound/1000
      
      for (i in 1:n_subjects) {
        #this function returns binary variables in the first columns followed by normal variables
        vars <- jointly.generate.binary.normal(no.rows=length(t), no.bin=1, no.nor=1,
                                               prop.vec.bin=c(mu2), mean.vec.nor=c(mu1),
                                               var.nor=c(sd1^2), corr.vec=c(bounded_corr))
        z1 <- vars[ ,2] 
        z2 <- vars[ ,1]   
        XB <- exp(b1*z1 + b2*z2)
        z.list[[i]] <- data.frame(id = i, z1 = z1,z2 =  z2, xB = XB)
      }
    }     
    
    #K function applies ACCEPT-REJECT algorithm
    k <- function(x, m, M, rates, t){
      ifelse(x <= m | x >= M, 0, dpexp(x, rates, t))
    }
    
    #define survival time generation function
    gen.y <- function(.x, .g_inv_t,  .g_inv_t_min, .g_inv_t_max) {
      .x1 <- .x$xB
      .d <- ppexp(.g_inv_t_max, .x1, .g_inv_t) - ppexp(.g_inv_t_min, .x1, .g_inv_t)
      .M <- 1 / .d
      .r <- 60
      .count<-0
      #counter of times repeat is run
      while (.count<1000) {
        .count <- .count+1
        .y <- rpexp(.r, .x1, .g_inv_t)
        .u <- runif(.r)
        .t <- .M * (k(.y, .g_inv_t_min, .g_inv_t_max, .x1, .g_inv_t) / .d / dpexp(.y, .x1, .g_inv_t))
        .y <- .y[.u <= .t][1]
        if (!is.na(.y)) {break}
      }
      .y
    }
    
    ## generate times in g_inv scale
    cat(paste0(WORKER.ID, "-", dataset_number, " setup complete at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
    y <- ldply(z.list, gen.y, g.inv.t,  g.inv.t.min, g.inv.t.max)
    cat(paste0(WORKER.ID, "-", dataset_number, " times generated at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
    
    #### check for error in generating the untransformed time values
    ### only continue running the program if generation of the times in the previous step works
    ### if iteration was repeated 1000 times the function breaks and there will be missing values
    success <- ifelse(any(is.na(y[,1]))==TRUE,0,1)
    if (success != 1 ) {cat(paste0(WORKER.ID, "-", dataset_number, " time generation failed for:", dataset_number), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)}
    
    if (success==1) {
      y$g.y <- g(y[ ,1])
      
      #######
      ####### Create four types of censoring 
      
      #file to keep run times
      runtime_S <- proc.time() - ptm0
      ptm1 <- proc.time()
      cat(paste0(WORKER.ID, "-", dataset_number, " generating uncensored data at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      
      ############## No Censoring
      ### everyone has the event
      #CREATING CENSORING INDICATOR
      d_none <- rep(1, n_subjects)
      
      #CREATING DATASET
      data_none <- NULL
      for (i in 1:n_subjects) {
        id.temp <- rep(i, ceiling(y$g.y[i]))
        time.temp <- c(1:ceiling(y$g.y[i]))
        time0.temp <- 0:ceiling(y$g.y[i] - 1)
        d.temp <- c(rep(0, length(time.temp) - 1), d_none[i])
        z.temp <- z.list[[i]][1:ceiling(y$g.y[i]), c(names(z.list[[i]]))]
        data.temp <- cbind(z.temp[,1:3], time.temp, time0.temp, d.temp)
        data_none <- rbind(data_none, data.temp)
      }
      colnames(data_none) <- c('id','z1','z2','t','t0','d_none')
      runtime_none <- proc.time() - ptm1
      
      data_none <- data.frame(data_none,runtime_S=runtime_S[3],runtime_none=runtime_none[3])
      cat(paste0(WORKER.ID, "-", dataset_number, " uncensored data generated at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      
      ############## Random Censoring
      ### people are censored with probability prop.censor
      cat(paste0(WORKER.ID, "-", dataset_number, " generating random censored data at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      
      ptm2 <- proc.time()  
      data_random2 <- NULL
      for(pc in 1:3){
        #CREATING CENSORING INDICATOR
        prop.cen <- switch(pc,0.2,0.5,0.8)
        d_random <- sample(0:1, n_subjects, replace = TRUE, prob = c(prop.cen, 1 - prop.cen))
        
        #CREATING DATASET
        data_random <- NULL
        for (i in 1:n_subjects) {
          time.temp <- c(1:ceiling(y$g.y[i]))
          d.temp <- c(rep(0, length(time.temp) - 1), d_random[i])
          data.temp <- cbind(d.temp)
          data_random <- rbind(data_random, data.temp)
        }
        
        colnames(data_random) <- c(paste("d_random",pc,sep=""))
        data_random2 <- cbind(data_random2,data_random)
      }
      #runtime for 1 computed as the avg time from running the three scenarios (for all 3 prop_censor)
      runtime_r <- (proc.time()-ptm2)[3]/3
      data_random3<-data.frame(data_random2)
      cat(paste0(WORKER.ID, "-", dataset_number, " random censored data generated at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      data_random3[1:10,]
      #Putting data together for admin and random censoring
      # survival times are the same what changes is the censoring indicator
      mydata <- cbind(data_none,data_random3,runtime_random=runtime_r)
      mydata <- data.frame(mydata)
      mydata$write_time <- Sys.time()
      mydata$job_id <- WORKER.ID
      cat(paste0(WORKER.ID, "-", dataset_number, " writing first set of data at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      
      write.table(mydata, paste(output_location, "/data", gm, lim, w, md, rv, ns, "_", dataset_number, ".csv", sep=""), sep=",", row.names = FALSE)
      cat(paste0(WORKER.ID, "-", dataset_number, " first set of data written at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      
      
      ########## Traditional Censoring
      ### we generate a censoring time and then take the minimum of the censored and survival time
      ### censoring distribution follows a weibull where the shape and scale where determined to give
      ### the correct percentage of censoring observations: 20%, 50%, and 80%
      
      ### percent censoring
      #       pc <- 1
      cat(paste0(WORKER.ID, "-", dataset_number, " starting modeled censoring at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      for(pc in 1:3){
        
        ptm3 <- proc.time()
        #get parameters for censoring distribution 
        ##### CHANGE THE NAME OF THIS FILE
        #         trad_cens <- read.csv(paste0(source_location, "/autoresultsTRAD_2014_06_06.csv"),header=T) 
        cens.nu <- trad_cens$trad_nu[trad_cens$pc==pc]
        cens.mediangoal <- trad_cens$trad_mdg[trad_cens$pc==pc]
        #         print("before")
        cat(paste0(WORKER.ID, "-", dataset_number, " generating censor times for ", pc, " at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
        
        y_c <- get.censoring.times(n_subjects, gm, lim, w, md = NA, rv, pc, cens.nu, cens.mediangoal)
        
        cat(paste0(WORKER.ID, "-", dataset_number, " censor times generated for ", pc, " at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
        #         print("after")
        y_1_modeled <- y
        
        d_modeled_1 <- rep(1, n_subjects)
        d_modeled_1[y_c$g.y < y_1_modeled$g.y] <- 0
        y_1_modeled$g.y[d_modeled_1 == 0] <- y_c$g.y[d_modeled_1 == 0]
        
        data_modeled <- NULL  
        for(i in 1:n_subjects) {
          id.temp <- rep(i, ceiling(y_1_modeled$g.y[i]))
          time.temp <- c(1:ceiling(y_1_modeled$g.y[i]))
          time0.temp <- 0:ceiling(y_1_modeled$g.y[i] - 1)
          d.temp <- c(rep(0, length(time.temp) - 1), d_modeled_1[i])
          z.temp <- z.list[[i]][1:ceiling(y_1_modeled$g.y[i]), c(names(z.list[[i]]))]
          data.temp <- cbind(z.temp[,1:3], time.temp, time0.temp, d.temp)
          data_modeled <- rbind(data_modeled, data.temp)
        }
        
        colnames(data_modeled) <- c('id','z1','z2','t','t0','d')
        runtime_t <- proc.time()-ptm3
        data_trad <- data.frame(data_modeled,runtime=runtime_t[3])
        data_trad$write_time <- Sys.time()
        data_trad$job_id <- WORKER.ID
        cat(paste0(WORKER.ID, "-", dataset_number, " writing data for ", pc, " at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
        
        write.table(data_trad, paste(output_location, "/data", gm, lim, w, md, rv, ns, pc, "_", dataset_number, "_trad.csv",sep=""), sep=",",row.names = FALSE)
        
        cat(paste0(WORKER.ID, "-", dataset_number, " data written for ", pc, " at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      }
      
      cat(paste0(WORKER.ID, "-", dataset_number, " modeled censoring complete at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      ############## Administrative Censoring
      ### the correct quantile of time was computed separately for each possible combination of parameters to give
      #### approximately 20%, 50% and 80% of censoring
      cat(paste0(WORKER.ID, "-", dataset_number, " start admin censoring at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      for(pc in 1:3){
        
        ptm4 <- proc.time()
        #get parameters for censoring distribution 
        #         admin_cens <- read.csv(paste0(source_location, "/admin_censor_times_2014_06_03.csv"),header=T) 
        #note that md=1 is simply because of a rerun of a previous program to generate the times 
        cens.quantile <- admin_cens$admin_time[admin_cens$pc==pc]
        
        y_1_admin <- y
        
        admin_time <- rep(cens.quantile,n_subjects)
        
        d_admin_1 <- rep(1, n_subjects)
        d_admin_1[admin_time < y_1_admin$g.y] <- 0
        y_1_admin$g.y[d_admin_1 == 0] <- y$g.y[d_admin_1 == 0]
        
        data_admin <- NULL 
        for(i in 1:n_subjects) {
          id.temp <- rep(i, ceiling(y_1_admin$g.y[i]))
          time.temp <- c(1:ceiling(y_1_admin$g.y[i]))
          time0.temp <- 0:ceiling(y_1_admin$g.y[i] - 1)
          d.temp <- c(rep(0, length(time.temp) - 1), d_admin_1[i])
          z.temp <- z.list[[i]][1:ceiling(y_1_admin$g.y[i]), c(names(z.list[[i]]))]
          data.temp <- cbind(id.temp,z.temp[,2:3],time.temp, time0.temp, d.temp)
          data_admin <- rbind(data_admin, data.temp)
        }   
        colnames(data_admin) <- c('id','z1','z2','t','t0','d')
        runtime_a <- proc.time()-ptm4
        data_admin <- data.frame(data_admin,runtime=runtime_a[3])
        data_admin$write_time <- Sys.time()
        data_admin$job_id <- WORKER.ID
        
        cat(paste0(WORKER.ID, "-", dataset_number, " writing data for ", pc, " at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
        write.table(data_admin, paste(output_location, "/data", gm, lim, w, md, rv, ns, pc, "_", dataset_number, "_admin.csv", sep=""), sep=",", row.names = FALSE)
        cat(paste0(WORKER.ID, "-", dataset_number, " data written for ", pc, " at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      }
      cat(paste0(WORKER.ID, "-", dataset_number, " admin censoring complete at: ", Sys.time(), "\n"), file = paste0(output_location, "/log/", WORKER.ID, "_output.txt"), append = T)
      
      
    } # close success if
  }, seed.init, n_subjects, gm, lim, w, md, rv, comb, source_location, output_location, trad_cens, admin_cens, ns, mediangoal, rho) #close l_ply
  #   }) # close system.time
  #   cat(paste0("Rep complete\n"), file = paste0(output_location, "/output.txt"), append = T)
  #   cat(paste0(unclass(results_time)[3],"\n"), file = paste0(output_location, "/runtime", WORKER.ID, ".txt"), append = T)
  #return(results_time) #close mysim function
}
