# Application-specific functions for sequential Monte Carlo.

# Function names should be kept unaltered, functions contents can be changed as
# long as they are internally consistent with each other.

# filtering
# Author: Luc Coffeng
# Created: September 13, 2018
# Adopted: October 12, 2018, By Peri Kontoroupis

# Function to run simulation model, possibly resulting in a dataset that is not
# yet of the same dimensions as the observed data.
  gen.data.pop <- function(theta.raw = NULL, theta = NULL, theta_trans = NULL) {

    # commented 11/10 by peri I need to check since reults are not logtranhormed
    if(is.null(theta) & is.null(theta.raw))
      stop("gen.data.pop requires either natural or transformed parameter values")

    if(is.null(theta)) {
      theta <- inv.transform.theta(theta_raw = theta.raw, theta_trans = theta_trans)
    }

  # Generate replicate dataset using vectorisation
    # commented 11/10 by peri
    # data.rep <- rsn(n = N.data * nrow(theta),
    #                 xi = rep(theta[, "loc"], each = N.data),
    #                 omega = rep(theta[, "scale"], each = N.data),
    #                 alpha = rep(theta[, "alpha"], each = N.data))

  # Split vector of replicate data into appropriate chunks and return
    # commented 11/10 by peri
    # i <- seq_along(data.rep)
    # split(data.rep, ceiling(i / N.data))

  # # For parallel simulation with more complex models:
    theta_list <- alply(.data = if(is.matrix(theta)) theta else t(theta),
                        .margins = 1,
                        .fun = as.list)


    foreach(i = 1:length(theta_list),
            .inorder = TRUE,
            .errorhandling = "pass",
            .packages = c("XML",
                          "data.table"),
            .export = c("wormsim.source.dir",
                        "wormsim.dir",
                        "input.schema")) %dopar% {

                          # Set paths and load functions (re-execute for parallel sessions)
                          setwd(wormsim.source.dir)
                          source("xml_substitute_functions.r")
                          source("create_xml_functions.r")
                          source("basic_functions.r")
                          source("par_grid_functions.r")
                          source("special_functions.r")

                          seed <- as.integer(runif(1, 0, 1e6))
                          #print(seed)
                          #  setwd( wormsim.dir )

                          # Run WORMSIM for parameter set i and return output
                      xx<-   run_wormsim(wormsim_dir = wormsim.dir,
                                      i_run = i,
                                      input_schema = input.schema,
                                      par_alt = theta_list[[i]],
                                      seed_start = seed,
                                      seed_end = seed)
                           #
                            if (  is.null(xx)==TRUE )

                            {
                            sink("Report.txt", append=TRUE) #open sink file and add output
                            print(xx$intensity)
                            print(theta_list[[i]])
                            print(seed)
                            result=c(i=i,seed=seed,theta_list[[i]][1:2])
                            write.table(result, file=paste("failed.runs.log",".txt", sep=""),
                                        sep="\t", row.names=F,append=TRUE)
                            sink()

                          #
                           }

                        run_wormsim(wormsim_dir = wormsim.dir,
                                    i_run = i,
                                    input_schema = input.schema,
                                    par_alt = theta_list[[i]],
                                    seed_start = seed,
                                    seed_end = seed)
                        }


  }


# Function to generate summary statistics from (replicate) dataset(s). Input is
# expected to be of identical dimensions as output from gen.data.obs.
  gen.summary <- function(x) {

    "prop" = x$prev

  }


# Function to calculate distance between replicate and observed data
  D <- function(y_rep, y_obs) {

  # Calculate summary estimates of observed and replicate data
    data_obs_summary <- gen.summary(y_obs)
    data_rep_summary <- llply(.data = y_rep,
                              .fun = gen.summary)

  # Optional: check validity of replicate data, i.e. disqualify simulations
  # in which the process "failed" for some reason
    penalty <- sapply(y_rep, function(x) {

      if ( length(x$summ) > 1 & (x$summ[year == min(year), get("w+")] == 0)  ) {
      #  the case when a simulation failed, disqualify using Inf to the distance
          Inf
        print("u are here")
      }

    })

  # Calculate and return final distance
    unname(sapply(data_rep_summary, function(x) {
      sum((x - data_obs_summary)^2)
    })) + penalty

  }


  sim_study <- function(x) {

    y1 <- extract_intens_prev(x$intensity[year == 2011 & age > 5])

    y1[age == 10, cases := rmvhyper(cases, study_size[age == "05-9", N])]
    y1[age == 20, cases := rmvhyper(cases, study_size[age == "10-19", N])]
    y1[age == 30, cases := rmvhyper(cases, study_size[age == "20-29", N])]
    y1[age == 40, cases := rmvhyper(cases, study_size[age == "30-39", N])]
    y1[age == 50, cases := rmvhyper(cases, study_size[age == "40-49", N])]
    y1[age == 99, cases := rmvhyper(cases, study_size[age == "50+", N])]

    x <- y1[, lapply(.SD, sum, na.rm=TRUE), by = intens, .SDcols=c("cases")]

    x[, prev := cases / sum(cases)]

    x = data.table(x)[, .(prev,intens)]

  }
### END OF CODE ###
