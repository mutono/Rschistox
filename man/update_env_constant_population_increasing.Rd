% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/schistox_functions.R
\name{update_env_constant_population_increasing}
\alias{update_env_constant_population_increasing}
\title{update env constant population increasing
   update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)
update the population for a given length of time. Here we include deaths and for each death an individual is immediately born.
Interventions are included in this function and larvae are uptaken immediately as worms and egg production follows a monotonically increasing function}
\usage{
update_env_constant_population_increasing(
  num_time_steps,
  humans,
  miracidia,
  cercariae,
  pars,
  mda_info,
  vaccine_info
)
}
\description{
update env constant population increasing
   update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)
update the population for a given length of time. Here we include deaths and for each death an individual is immediately born.
Interventions are included in this function and larvae are uptaken immediately as worms and egg production follows a monotonically increasing function
}
