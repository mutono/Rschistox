library(devtools)
#use_git()

## load all R scripts
#load_all()

## define functions for schistox

## functions


#' Human -------------------------------------------------------------------


#' This struct contains the information about a human individual. This contains age, the pre determined age of death, community they are in,
#'    their gender, predisposition to picking up cercariae, the number of larvae, female and male worms and eggs in the individual along with
#'    a count of total lifetime eggs. Also it has their age dependent contact rate, adherence and access to interventions.



#' @export
Human <- function(age,death_age,gender,predisposition,female_worms,male_worms,
                  eggs,vac_status,age_contact_rate,adherence,access,community,
                  relative_contact_rate,uptake_rate,acquired_immunity,
                  total_worms, larvae,last_uptake){
  age <- as.numeric(age)
  death_age <- as.numeric(death_age)
  gender <- as.integer(gender)
  predisposition <- as.numeric(predisposition)
  female_worms <- as.array(as.integer(female_worms))
  male_worms <- as.array(as.integer(male_worms))
  eggs <- as.integer(eggs)
  vac_status <- as.integer(vac_status)
  age_contact_rate <- as.numeric(age_contact_rate)
  adherence <- as.integer(adherence)
  access <- as.integer(access)
  community <- as.integer(community)
  relative_contact_rate <- as.numeric(relative_contact_rate)
  uptake_rate <- as.numeric(uptake_rate)
  acquired_immunity <- as.numeric(acquired_immunity) # level of acquired immunity
  total_worms <- as.integer(total_worms) # total number of worms over lifetime
  larvae <- as.array(as.integer(larvae))
  last_uptake <- as.integer(last_uptake)


  human <- list("age"=age,"death_age"=death_age, "gender"=gender, "predisposition"=predisposition,
                "female_worms"=female_worms, "male_worms"=male_worms, "eggs"=eggs, "vac_status"=vac_status,
                "age_contact_rate"=age_contact_rate, "adherence"=adherence, "access"=access, "community"=community,
                "relative_contact_rate"=relative_contact_rate, "uptake_rate"=uptake_rate, "acquired_immunity"=acquired_immunity,
                "total_worms"=total_worms, "larvae"=larvae, "last_uptake"=last_uptake)

  return(human)

}




#' Parameters function
#'
#' This function takes in all the parameters for the model
#'
#' @param N human population size
#'
#' @param time_step length of time step (in days)
#' @param N_communities number of communities in the population sharing the same environmental source
#' @param community_probs probability of being in each community
#' @param community_contact_rate contact rate with the environment for each of the commununity
#' @param density_dependent_fecundity decrease in egg production per worm due to high density of worms
#' @param average_worm_lifespan average expectancy of a worm
#' @param max_age maximum age of individual
#' @param initial_worms initial no. of worms
#' @param initial_miracidia  initial no. of miracidia in the environment
#' @param initial_miracidia_days no.of days miracidia will age into cercariae larvae
#' @param init_env_cercariae initial no of cercaria in the environment
#' @param worm_stages number of stages in the worm. Having 1 stage will result to a Gamma distribution
#' @param contact_rate global contact rate for the uptake of larvae from the environment
#' @param max_fec_contact_rate_product product of max fecundity and the contact rate in the population. Setting this to a desired value is often a good way to ensure that the epidemic stays within a reasonable range, as when the max fecundity increases, if the contact rate doesn't decrease appropriately, then the behaviour of the outbreak can be unrealistically difficult to control.
#' @param max_fecundity expected no. of eggs from a single worm
#' @param age_contact_rates contact rate for the uptake of larvae from the environment for the chosen age groups
#' @param ages_for_contacts  age groups for specifying contact rates
#' @param contact_rate_by_age_array <- rep(0,times=max_age+1) array holding contact rate for each age
#' @param mda_adherence proportion of people who adhere to the MDA
#' @param mda_access proportion of people who have access to the MDA
#' @param female_factor factor for altering the contact rate for females, if we choose to have gender-specific behavior which affects contact rate
#' @param male_factor factor for altering the contact rate for males, if we choose to have gender-specific behavior which affects contact rate
#' @param miracidia_maturity no of days after which miracidias will mature to cercariae
#' @param birth_rate rate of birth of humans
#' @param human_cercariae_prop proportion of cercariae which are able to infect humans
#' @param predis_aggregation aggregation parameter for Poisson distributed egg production
#' @param cercariae_survival proportion of cercariae that survive from one time point to the next
#' @param miracidia_survival proportion of miracidia that survive from one time point to the next
#' @param death_prob_by_age probability of dying each year, specified by age
#' @param ages_for_death age ranges for death probailities
#' @param r aggregation parameter for negative binomially distributed egg production
#' @param vaccine_effectiveness efficacy of a vaccine if one is used
#' @param drug_effectiveness efficacy of a drug given during MDA
#' @param spec_ages number of individuals by age group
#' @param ages_per_index how many different ages we include in the spec_ages parameter
#' @param record_frequency how often we should record the prevalence in the population dusing simulation
#' @param use_kato_katz if 1, use Kato-Katz for egg count, if 0, do not use KK
#' @param kato_katz_par parameter for Gamma distribution if KK is used
#' @param heavy_burden_threshold number of eggs at which an individual is said to have a heavy infection
#' @param rate_acquired_immunity rate at which immunity will be acquired for individuals. This will be multiplied by the cumulative nymber of worms people have had during their life to decide the level of immunity acquired
#' @param M0 if a particular formula of egg production is used, this parameter is required and is a proxy for mean worm burden
#' @param human_larvae_maturity_time length of time (in days) after which a cercariae uptake by a human will mature into a worm
#' @param egg_sample_size the proportion of eggs which are sampled from each individual every time we check their burden (between 0 and 1). 1= all eggs in the person are sampled. Typical value fpr a urine sample may be ~1/100
#' @param input_ages input ages for contructing contact array
#' @param input_contact_rates input contact rates
#' @param scenario can be one of "low adult", "moderate adult" or high adult"
#'
#' @export
Parameters <- function(N, time_step, N_communities, community_probs,
                       community_contact_rate, density_dependent_fecundity,
                       average_worm_lifespan, max_age, initial_worms,
                       initial_miracidia, initial_miracidia_days, init_env_cercariae,
                       worm_stages, contact_rate,max_fec_contact_rate_product, max_fecundity, age_contact_rates,
                       ages_for_contacts, contact_rate_by_age_array, mda_adherence,
                       mda_access, female_factor, male_factor, miracidia_maturity,
                       birth_rate, human_cercariae_prop, predis_aggregation, cercariae_survival,
                       miracidia_survival, death_prob_by_age, ages_for_death, r,
                       vaccine_effectiveness, drug_effectiveness, spec_ages, ages_per_index,
                       record_frequency, use_kato_katz, kato_katz_par, heavy_burden_threshold,
                       rate_acquired_immunity, M0, human_larvae_maturity_time, egg_sample_size, input_ages, input_contact_rates,
                       scenario)
{
  N <- as.integer(N)
  time_step <- as.numeric(time_step)
  N_communities <- as.integer(N_communities)
  community_probs <- as.array(as.numeric(community_probs))
  community_contact_rate <-as.array(as.numeric(community_contact_rate))
  density_dependent_fecundity <- as.numeric(density_dependent_fecundity)
  average_worm_lifespan <- as.numeric(average_worm_lifespan)
  max_age <- as.numeric(max_age)
  initial_worms <- as.integer(initial_worms)
  initial_miracidia <- as.integer(initial_miracidia)
  initial_miracidia_days <- as.integer(initial_miracidia_days)
  init_env_cercariae <- as.integer(init_env_cercariae)
  worm_stages <- as.integer(worm_stages)
  contact_rate <-as.numeric(contact_rate)
  max_fec_contact_rate_product <- as.numeric(max_fec_contact_rate_product)
  max_fecundity <- as.numeric(max_fecundity)
  age_contact_rates <- as.array(as.numeric(age_contact_rates))
  ages_for_contacts <- as.array(as.integer(ages_for_contacts))
  contact_rate_by_age_array <- as.array(as.numeric(contact_rate_by_age_array))
  mda_adherence <- as.numeric(mda_adherence)
  mda_access <- as.numeric(mda_access)
  female_factor <- as.numeric(female_factor)
  male_factor <- as.numeric(male_factor)
  miracidia_maturity <- as.integer(miracidia_maturity)
  birth_rate <- as.numeric(birth_rate)
  human_cercariae_prop <- as.numeric(human_cercariae_prop)
  predis_aggregation <- as.numeric(predis_aggregation)
  cercariae_survival <- as.numeric(cercariae_survival)
  miracidia_survival <- as.numeric(miracidia_survival)
  death_prob_by_age <- as.array(as.numeric(death_prob_by_age))
  ages_for_death <- as.array(as.numeric(ages_for_death))
  r <- as.numeric(r)
  vaccine_effectiveness <- as.numeric(vaccine_effectiveness)
  drug_effectiveness <- as.numeric(drug_effectiveness)
  spec_ages<- as.array(as.numeric(spec_ages))
  ages_per_index <- as.integer(ages_per_index)
  record_frequency <- as.numeric(record_frequency)
  use_kato_katz <- as.integer(use_kato_katz) # if 0, then don't use KK, if 1, use KK
  kato_katz_par <- as.numeric(kato_katz_par)
  heavy_burden_threshold <- as.integer(heavy_burden_threshold)
  rate_acquired_immunity <- as.numeric(rate_acquired_immunity)
  M0 <- as.numeric(M0) # mean worm burden
  human_larvae_maturity_time <- as.integer(human_larvae_maturity_time)
  egg_sample_size <- as.numeric(egg_sample_size)
  egg_production_distribution <- as.character(egg_production_distribution)
  input_contact_rates <- as.array(as.numeric(input_contact_rates))
  input_ages <- as.array(as.numeric(input_ages))
  scenario <- as.character(scenario)


  parameters <- list("N"=N, "time_step"=time_step, "N_communities" = N_communities,"community_probs" = community_probs,
                     "community_contact_rate" = community_contact_rate, "density_dependent_fecundity"=density_dependent_fecundity,
                     "average_worm_lifespan"=average_worm_lifespan, "max_age"=max_age, "initial_worms"=initial_worms,
                     "initial_miracidia"=initial_miracidia, "initial_miracidia_days"=initial_miracidia_days, "init_env_cercariae"=init_env_cercariae,
                     "worm_stages"=worm_stages, "contact_rate"=contact_rate, "max_fec_contact_rate_product"=max_fec_contact_rate_product,
                     "max_fecundity"=max_fecundity, "age_contact_rates" = age_contact_rates, "ages_for_contacts" = ages_for_contacts,
                     "contact_rate_by_age_array" = contact_rate_by_age_array, "mda_adherence"=mda_adherence, "mda_access"=mda_access,
                     "female_factor"=female_factor, "male_factor"=male_factor, "miracidia_maturity"=miracidia_maturity,
                     "birth_rate"=birth_rate, "human_cercariae_prop"=human_cercariae_prop, "predis_aggregation"=predis_aggregation,
                     "cercariae_survival"=cercariae_survival, "miracidia_survival"= miracidia_survival, "death_prob_by_age" = death_prob_by_age,
                     "ages_for_death" = ages_for_death, "r"=r, "vaccine_effectiveness"=vaccine_effectiveness, "drug_effectiveness"=drug_effectiveness,
                     "spec_ages" = spec_ages, "ages_per_index"=ages_per_index, "record_frequency"=record_frequency, "use_kato_katz"=use_kato_katz,
                     "kato_katz_par" = kato_katz_par, "heavy_burden_threshold"=heavy_burden_threshold, "rate_acquired_immunity" = rate_acquired_immunity,
                     "M0" = M0, "human_larvae_maturity_time"=human_larvae_maturity_time, "egg_sample_size" = egg_sample_size,
                     "egg_production_distribution" = egg_production_distribution, "input_ages"= input_ages, "input_contact_rates"=input_contact_rates, "scenario"=scenario)

  return(parameters)
}


#'  Out function -----------------------------------------------------------

#' This function contains the different outputs we are interested in recording. This is the
#' overall population burden, with categories for low, moderate and heavy burdens, along with
#' separate categories for the school age children and adults. Along with these, the time of each
#' result is recorded, so we can subsequently see the prevalence of the outbreak over time.
#' @export
out <- function(population_burden,sac_burden,adult_burden, pop_prev, sac_prev,
                adult_prev, sac_pop, adult_pop, final_ages, recorded_eggs, time){

  return(data.frame(population_burden,sac_burden,adult_burden,
                    pop_prev, sac_prev, adult_prev, sac_pop, adult_pop,
                    final_ages, recorded_eggs, time))
}


#' mda_information ---------------------------------------------------------


#' This function contains the information for the mda, storing the coverage, minimum and maximum age targeted,
#' gender, drug efficacy and the time for the mda to be done
#' @export
mda_information <- function(mda_information, coverage,min_age,
                            max_age, gender, effectiveness, time)
{
  mda_info <- data.frame(mda_information, coverage,min_age,
                         max_age, gender, effectiveness, time)
  return(mda_info)
}


#' vaccine_information -----------------------------------------------------

#' This struct contains the information for the vaccine, storing the coverage, minimum and maximum age targeted,
#' gender, drug efficacy and the time for the vaccine to be done along with how long the vaccine provides protection for
#' @export
vaccine_information <- function(coverage,min_age,max_age,gender,duration,
                                time){
  vacc_info <- data.frame(coverage,min_age,max_age,gender,duration,
                          time)
  return(vacc_info)
}


#' Contact settings --------------------------------------------------------


#' create the age specific contact settings given the scenario


#' This will create age dependent contact rates based on the scenario for simulation which is input. This is either
#'    "low adult", "moderate adult" or "high adult"
#' @export
create_contact_settings <- function(scenario)
{
  scenario <- as.character(scenario)
  if (scenario %in% "low adult")
    pars$age_contact_rates = c(0.01, 1.2, 1, 0.02)

  else if (scenario %in% "moderate adult")
    pars$age_contact_rates = c(0.032, 0.61, 1, 0.06)
  else
    pars$age_contact_rates = c(0.01, 0.61, 1, 0.12)

  return(pars$age_contact_rates)

}


#' Age dependent contact rate ----------------------------------------------


#' function to get age dependent contact rate.
#' the contact rates are taken from the
#' "What is required in terms of mass drug administration to interrupt the transmission
#'     of schistosome parasites in regions of endemic infection?" paper
#' at some point we may change this to be an input from a file instead
#' @export
make_age_contact_rate_array <- function(pars, scenario, input_ages, input_contact_rates)
{
  #' This will make the contact rate array for each age of individual in the population, based either on a scenario basis
  #'    ("low adult", "moderate adult" or "high adult"), through the create_contact_settings function, or through
  #'    specifying an array of age breaks and the desired contact rates for the ages specified by the ages, using the
  #'        input_ages and input_contact_rates variables.
  #'    For example: input_contact_rates = [0.02,0.61, 1,0.06], input_ages= [4,9,15,100] will make 0-4 year olds have contact rate 0.02,
  #'    5-9 will have rate 0.61, 10-15 rate 1 and 16+ 0.06

  if (pars$max_age < 60){
    error("max_age must be greater than 60")
  }
  #if (length(pars$input_ages) == 0) {
    contact_settings = as.array(as.numeric())
    contact_settings = create_contact_settings(scenario)
    x <- rep(contact_settings[4], times=pars$max_age+1)
  #}

  #' initialize an array with the same value for contact rate across all ages
  pars$contact_rate_by_age_array <- rep(contact_settings[4], times=pars$max_age+1)
  #'         pars.contact_rate_by_age_array = pars.contact_rate_by_age_array[1]

  #' then edit the entries for different ages according to values
  for (i in 1:5)
  {
    pars$contact_rate_by_age_array[i] = contact_settings[1]
  }

  if (scenario == "high adult") {
    for (i in 6:12)
    {
      pars$contact_rate_by_age_array[i] = contact_settings[2]

    }
    for (i in 13:21){
      pars$contact_rate_by_age_array[i] = contact_settings[3]
    }

  } else {
    for (i in 6:10){
      pars$contact_rate_by_age_array[i] = contact_settings[2]
    }
    for (i in 11:16){
      pars$contact_rate_by_age_array[i] = contact_settings[3]
    }
    for (i in 1:length(pars$contact_rate_by_age_array)){
      pars$contact_rate_by_age_array[i] = input_contact_rates[i]
    }
    for (i in 1 : length(input_contact_rates))
    {
      if (i == 1){
        for (j in 1:trunc(input_ages[i] + 1))
        {
          pars$contact_rate_by_age_array[j] = input_contact_rates[i]
        }
      } else{
        for (j in trunc(input_ages[i-1]+2):trunc(input_ages[i] + 1)){
          pars$contact_rate_by_age_array[j] = input_contact_rates[i]
        }
      }
    }
  }

  return(pars)
}



#' Age for death of an individual ------------------------------------------


#' function to generate an age for death of an individual

#' This will create the initial human population with randomly chosen age, and gender.
#' Predisposition is taken to be gamma distributed
#' There is also a male and female adjustment to predisposition adjusting for gender specific behaviour
#' In addition to this, it will create the initial miracidia environment vector
#' @export
get_death_age <- function(pars){
  age <- 0
  k <- 1
  p <- pars$death_prob_by_age[k]
  goal_age <- pars$ages_for_death[k]
  x = rgamma(1,1)
  while (x > p)
    age <- age+ 1
  if (age >= goal_age){
    k <- k + 1
    p <- pars$death_prob_by_age[k]
    goal_age <- pars$ages_for_death[k]
  }
  else{
    x = rgamma(1,1)
  }
  death_age <- age + rgamma(1,1)

  return(death_age)
}


#' age population and generating death ages --------------------------------

#' function to age population and generating death ages

#' Step forward the population by a number of steps, where we will go through aging and removing individuals when they
#' pass their age of death. This will generate an age distribution in the population which corresponds to the death_prob_by_age and
#' ages_for_deaths parameters, which specify the probability of dying at each age.
#' @export
generate_ages_and_deaths <- function(num_steps, humans, pars){
  for (i in 1:num_steps){
    x <- array()
    for (i in 1:length(humans)){
      if(humans[i]$age - humans[i]$death_age > 0){
        x[i]
      }
      else{
        humans[i]$age <- humans[i]$age + (pars$time_step/365)
      }
    }
    k <- length(x)
    if (k > 0){
      for (j in seq(from=k,to=1, by=-1)){
        humans <- humans[-x[j]]
        humans <- birth_of_human(human, pars)
      }
    }
  }

  return(humans)
}


#' Initial population ------------------------------------------------------


#' This will create the initial human population with randomly chosen age, and gender.
#' Predisposition is taken to be gamma distributed
#' There is also a male and female adjustment to predisposition adjusting for gender specific behaviour
#' In addition to this, it will create the initial miracidia environment vector
#' @export
create_population <- function(pars){

  if (length(pars$community_probs) != pars$N_communities){
    error("must provide probabilities for membership of each community")
  }
  else if (pars$N_communities > 1){
    community_selection = cumsum(pars$community_probs)/sum(pars$community_probs)
  }
  else {
    community_selection = 1
  }

  humans <- list()
  humans1 <- list()

  #'  initialize and fill the environmental variable

  miracidia<- arrray(integer())
  for (i in 1 : pars$initial_miracidia_days){
    miracidia[i] <- pars$initial_miracidia
  }

  cercariae <- pars$init_env_cercariae

  #'  initialize the Gamma distribution for predisposition selection
  pre <- rgamma(pars$predis_aggregation, 1/pars$predis_aggregation)
  #' select all predispositions
  predisposition = sample(c(pre, pars$N))


  for (i in 1:pars$N){

    f_worms <- rep(0, times=pars$worm_stages)
    f_worms[1] = round(runif(1)*pars$initial_worms)
    m_worms = rep(0, times=pars$worm_stages)
    m_worms[1] = round(runif(1)*pars$initial_worms)

    if (runif(1) > pars$mda_adherence)
    {
      adherence = 0
    }

    else
    {
      adherence = 1
    }

    if (runif(1) > pars$mda_access)
    {
      access = 0
    }

    else
    {
      access = 1
    }

    community = community_selection[community_selection > runif(1)][1]
  }
  # Human(age,gender, predisposition, female_worms, male_worms, eggs, vac_status, age_contact_rate, adherence, access, community)
  death_age <- get_death_age(pars)
  humans1 <-Human(pars$max_age*runif(1), death_age, sample(c(0,1)), predisposition[i],
                  f_worms, m_worms,
                  0, 0, 0, adherence, access, community, 0, 0,0,0, NA, 0)
  humans <- rbind( humans,humans1)

  age = trunc(as.numeric(tail(data.frame(humans)$age,1)))

  contact_rate_age = pars$ages_for_contacts[pars$ages_for_contacts> age][1]


  as.numeric(tail(data.frame(humans)$age_contact_rate,1)) = pars$age_contact_rates[contact_rate_age]
  as.numeric(tail(data.frame(humans)$relative_contact_rate,1)) = as.numeric(tail(data.frame(humans)$age_contact_rate,1)) *  pars$community_contact_rate[as.numeric(tail(data.frame(humans)$community,1))]/
    (max(pars$community_contact_rate) * max(pars$contact_rate_by_age_array))


  if (tail(as.array(data.frame(humans)$gender)[[1]],1) == 0){
    as.numeric(tail(data.frame(humans)$predisposition,1)) = as.numeric(tail(data.frame(humans)$predisposition,1)) * pars$female_factor
  }else if (tail(as.array(data.frame(humans)$gender)[[1]],1) == 1){
    as.numeric(tail(data.frame(humans)$predisposition,1)) = as.numeric(tail(data.frame(humans)$predisposition,1)) * pars$male_factor
  }
  else{
    as.numeric(tail(data.frame(humans)$uptake_rate,1)) = as.numeric(tail(data.frame(humans)$predisposition,1)) * pars$contact_rate * as.numeric(tail(data.frame(humans)$age_contact_rate,1)) *
      pars$community_contact_rate[community]
  }
  humans = generate_ages_and_deaths(20000, humans, pars)
  humans = update_contact_rate(humans,  pars)
  return(humans, miracidia, cercariae)

}






#' Initial human population with an age distribution ------------------------

#' This will create the initial human population with an age distribution
#' specified by the spec_ages variable
#' Predisposition is taken to be gamma distributed. There is also a male and female
#' adjustment to predisposition adjusting for gender specific behaviour
#' In addition to this, it will create the initial miracidia environment vector
#' @export
create_population_specified_ages <- function(pars){

  if (length(pars$community_probs) != pars$N_communities){
    error("must provide probabilities for membership of each community")
  }
  else{
    community_selection <- 1
  }
  community_selection <- cumsum(pars$community_probs)/sum(pars$community_probs)

  ages <- specified_age_distribution(pars)

  humans <- array()
  #'  initialize and fill the environmental variable
  miracidia <- array(as.integer())
  for (i in 1 : pars$initial_miracidia_days){
    miracidia <- as.array(append(miracidia, pars$initial_miracidia))
  }
  cercariae <- pars$init_env_cercariae
  #' initialize the Gamma distribution for predisposition selection
  pre <- rgamma(pars$predis_aggregation, 1/pars$predis_aggregation)

  #' select all predispositions
  predisposition <- runif(1,min=pre, max=pars$N)


  for (i in 1:pars$N){

    f_worms <- array(rep(0, pars$worm_stages))
    f_worms[1] <- trunc(round(runif(1)*pars$initial_worms))
    m_worms <- array(rep(0, pars$worm_stages))
    m_worms[1] <- trunc(round(runif(1)*pars$initial_worms))

    adherence <- 1 - 1*(runif(1) > pars$mda_adherence)

    access <- 1 - 1*(runif(1) > pars$mda_access)


    community <- community_selection[community_selection> runif(1)][1]

    death_age <- get_death_age(pars)
    humans<- append(humans, Human(ages[i]+rand(), death_age, runif(1,min=0,max=1), predisposition[i],
                                 f_worms, m_worms,
                                 0, 0, 0, adherence, access, community, 0, 0,0,0, 0, 0))

    age <- trunc(as.numeric(tail(data.frame(humans)$age,1)))

    contact_rate_age = pars$ages_for_contacts[pars$ages_for_contacts>age][1]


    as.numeric(tail(data.frame(humans)$age_contact_rate,1)) = pars$age_contact_rates[contact_rate_age]
    as.numeric(tail(data.frame(humans)$relative_contact_rate,1)) = as.numeric(tail(data.frame(humans)$age_contact_rate,1)) *  pars$community_contact_rate[as.numeric(tail(data.frame(humans)$community,1))]/
      (max(pars$community_contact_rate) * max(pars$contact_rate_by_age_array))
    as.numeric(tail(data.frame(humans)$predisposition,1)) = as.numeric(tail(data.frame(humans)$predisposition,1)) * ((1-tail(as.array(data.frame(humans)$gender)[[1]],1))* pars$female_factor) + (tail(as.array(data.frame(humans)$gender)[[1]],1)* pars$male_factor)
    as.numeric(tail(data.frame(humans)$uptake_rate,1)) = as.numeric(tail(data.frame(humans)$predisposition,1)) * pars$contact_rate * as.numeric(tail(data.frame(humans)$age_contact_rate,1)) *
      pars$community_contact_rate[community]

  }
  humans <- generate_ages_and_deaths(20000, humans, pars)
  humans <- update_contact_rate(humans,  pars)
  return(humans, miracidia, cercariae)



}



#' update contact rates ----------------------------------------------------

#' function to update the contact rate of individuals in the population. This is necessary
#' as over time when people age, they will move through different age groups which have
#' different contact rates

#' @export
update_contact_rate <- function(humans,  pars){


  age = min(length(pars$contact_rate_by_age_array) - 1, (trunc(humans$age)))
  humans$age_contact_rate =  pars$contact_rate_by_age_array[age+1]
  humans$relative_contact_rate = humans$age_contact_rate *  pars$community_contact_rate[humans$community]/
    (max(pars$community_contact_rate) * max(pars$contact_rate_by_age_array))

  humans$uptake_rate = humans$predisposition * pars$contact_rate * humans$age_contact_rate *
    pars$community_contact_rate[humans$community]

  return(humans)

}

  #' Cercariae uptake --------------------------------------------------------




#' cercariae_uptake(humans, cercariae, miracidia, pars)

#' uptake cercariae into humans, whilst updating cercariae with miracidia.
#' Uptaken cercariae immediately become worms in this formulation
#' @export
cercariae_uptake <-function(humans, cercariae, miracidia, pars){
  #cercariae = 0
  humans = sample(humans)
  k = length(humans)
  #= assign larvae which have been in the environment for 40 days to become infective.
  #then delete those larvae from the environmental larvae =#

    if (length(miracidia) > (pars$miracidia_maturity / pars$time_step)){
      cercariae = cercariae+(miracidia[1] * pars$human_cercariae_prop)
      splice(miracidia, 1)
    }
   else{
     splice(miracidia, 1)
   }




  #' if there are still infective larvae in the environment,
 #' we will uptake from choose from Poisson
#'  distribution. otherwise, just uptake 0. =#
    #'   if cercariae > 0
    #' calculate the rate of the poisson distribution
    humans$last_uptake = humans$last_uptake+pars$time_step
  if (humans$uptake_rate > 0)
  {
    pois_rate  = max(humans$uptake_rate * (1-pars$rate_acquired_immunity * humans$total_worms) *  cercariae / k, 0)

  }
   else{#' reduce the rate according to the effectiveness of the vaccine (if any is given)
     pois_rate = pois_rate * (1 - (humans$vac_status > 0) * pars$vaccine_effectiveness)
   }


  #' choose from the Poisson distribution
  uptake = rpois(1,pois_rate)
  if (uptake > 0){
    humans$last_uptake = -1
  }
  else{
    n = sample(rbinom(uptake,1, 0.5))

  }
  #println(uptake)

  humans$female_worms[1] =humans$female_worms[1]+ n
  humans$male_worms[1] = humans$male_worms[1]+ (uptake - n)
  humans$total_worms = humans$total_worms+ uptake

  #' reduce the infective larvae by the number of larvae uptaken
  cercariae = cercariae-uptake
  cercariae = max(0, cercariae)

  return(humans, cercariae, miracidia)


}
  #' cercariae uptake human larvae -------------------------------------------



#'uptake cercariae into humans, whilst updating cercariae with matured miracidia.
#'Uptaken cercariae become larvae within humans, rather than immediately into worms with this function.
#' @export
cercariae_uptake_with_human_larvae <-function(humans, cercariae, miracidia, pars){

  cercariae = 0
  humans = sample(humans)
  k = length(humans)
  #' assign larvae which have been in the environment for 40 days to become infective.
  #' then delete those larvae from the environmental larvae =#

  if (length(miracidia) > (pars$miracidia_maturity / pars$time_step))
  {
    cercariae =cercariae + (miracidia[1] * pars$human_cercariae_prop)
    splice(miracidia, 1)
  }
  else{
    splice(miracidia, 1)
  }


  #' if there are still infective larvae in the environment,
 #' we will uptake from choose from Poisson
#'  distribution. otherwise, just uptake 0. =#
    #'   if cercariae > 0
    #' calculate the rate of the poisson distribution
    humans$last_uptake =humans$last_uptake+ pars$time_step
  if (humans$uptake_rate > 0)
  {
    pois_rate  = max(humans$uptake_rate * (1-pars$rate_acquired_immunity * humans$total_worms) *  cercariae / k, 0)

  }
 else{
   #' reduce the rate according to the effectiveness of the vaccine (if any is given)
   pois_rate = pois_rate * (1 - (humans$vac_status > 0) * pars$vaccine_effectiveness) * (1 - humans$acquired_immunity)

 }

  #' choose from the Poisson distribution
  uptake = rpois(1,pois_rate)

  humans$larvae<- append(humans$larvae, uptake)
  if (uptake > 0)
  {
    humans$last_uptake = 0
  }
 else{
   humans$last_uptake = humans$last_uptake
 }
  #' reduce the infective larvae by the number of larvae uptaken
  cercariae = cercariae-uptake
  cercariae = max(0, cercariae)
  return(humans, cercariae, miracidia)

}
#' @export
enact_maturity_function_false <-  function(h){
  return(h)
}

#' @export
enact_maturity_function_true <-  function(h){
  females = rbinom(1,h.larvae[1], 0.5)
  humans$female_worms[1] =humans$female_worms[1]+ females
  humans$male_worms[1] = humans$male_worms[1] + ( humans$larvae[1] - females)
  humans$total_worms = humans$total_worms+ humans$larvae[1]
  slice(humans$larvae, 1)
  return(humans)
}




  #' Human larvae maturity ---------------------------------------------------



#' This will mature the human larvae into worms after a chosen number of days, which is specified
#' by the human_larvae_maturity_time parameter in the pars struct

  # function human_larvae_maturity(humans, pars)
  #
  #     a = "enact_maturity_function_"
  # #=  loop over humans  =#
  #     @inbounds for h in humans
  #
  # #=  if we there are non-zero larvae over the age of 35 days, then add
  #      these to worms and remove from human_cercariae  =#
  #         fn = a * string(length(h.larvae) > round(pars.human_larvae_maturity_time/pars.time_step, digits = 0))
  #         h = getfield(Schistoxpkg, Symbol(fn))(h)
  #     end
  #
  # #= return arrays  =#
  #     return humans
  # end
#' @export
human_larvae_maturity <- function(humans, pars){


  #' if we there are non-zero larvae over the age of 35 days, then add
  #' these to worms and remove from human_cercariae  =#
    if (length(humans$larvae) > round(pars$human_larvae_maturity_time/time_step, digits = 0)){
      females = rbinom(1, humans$larvae[1], 0.5)
      humans$female_worms[1] =humans$female_worms[1] + females
      humans$male_worms[1] = humans$male_worms[1] + (humans$larvae[1] - females)
      humans$total_worms = humans$total_worms + humans$larvae[1]
      slice(humans$larvae, 1)
    }
   else
   {
     humans$female_worms[1]=humans$female_worms[1]
   }
  #' return arrays
  return(humans)

}

#' Kill miracidia in the environment ---------------------------------------



#' Kill a chosen proportion of miracidia in the environment governed by the
#' miracidia_survival parameter in the pars struct
#' @export
miracidia_death <-function(miracidia, pars){

  #' as miracidia is an array, we need to use . syntax to apply
  #' the functions to each element in the array =#
    #' return rand.(Binomial.(env_miracidia, 1 - env_miracidia_survival_prop))
    if (pars$miracidia_survival <=0){
      print("miracidia_survival_prop must be bigger than 0")
    }
  else
  {
    tail(miracidia,1) = round(tail(miracidia,1) * pars$miracidia_survival)
  }
  return(miracidia)

}
  #' kill cercariae in the environment


  #' function to kill cercariae in the environment

#    cercariae_death!(miracidia, pars)

#' Kill a chosen proportion of cercariae in the environment governed by the
#' cercariae_survival parameter in the pars struct. This parameter governs what
#' proportion of cercariae survive for for one additional day, so if the time step
#' is greater than one, we have to calculate the correct proportion who die over
#' the chosen time step
#' @export
cercariae_death <-function(cercariae, miracidia, pars){

  updated_cercariae = cercariae
  if (cercariae > 0)
  {
    time_step_specific_cerc_death = pars$cercariae_survival
  }

  else if (pars$time_step > 1){
    c = cercariae
    for (i in 1:pars$time_step)
    {
      #' calculate the number of cercariae remaining if there was
      #' daily addition and culling of cercariae to the pool
      c = (c + miracidia[1]/pars$time_step) * pars$cercariae_survival

    }

  }
  c1 = cercariae
  c1 = (c1 +  miracidia[1])
  if (c1 > 0)
  {
    time_step_specific_cerc_death = c/c1
  }
  #' calculate the death rate needed to match the daily version of
  #' adding and culling cercariae

  if (pars$cercariae_survival <= 0)
  {
    print("cercariae_survival_prop must be bigger than 0")
  }

  else
  {
    updated_cercariae = round(cercariae * time_step_specific_cerc_death)

  }

  return(updated_cercariae)


}

  #' Kill worms within human host --------------------------------------------



  #' Worms die have a specified mean life span, and hence a rate of deaths per day .
  #' the number of deaths, or maturing from each stage is dependent
  #' on the number of worm stages and rate of aging through stages and dying.
  #' This p is multiplied by the time scale of the simulation, so if 2 days
  #' pass between consecutive time points, twice as many worms age and die

  #' for human in humans
  #'     println(human.eggs)
  #' end

#    worm_maturity!(humans, pars)

#' function to kill worms within human hosts, and if there is more than one stage for worm life,
#' to update how many worms are in each stage
#' @export
worm_maturity <-function(humans, pars){


  #' probability of aging out of category/ dying
  p = pars$time_step * pars$worm_stages/ (365 * pars$average_worm_lifespan)


  if (humans$last_uptake > 1.5 * 365*pars$average_worm_lifespan){


  for (j in 1:pars$worm_stages)
  {
  humans$female_worms[j] = 0
  humans$male_worms[j] = 0
  }
}
  else{

    #' kill appropriate number of worms in the final stage
    n = humans$female_worms[pars$worm_stages]
  humans$female_worms[pars$worm_stages] = rbinom(1, n, 1-p)

  n = humans$male_worms[pars$worm_stages]
  humans$male_worms[pars$worm_stages] = rbinom(1, n, 1-p)

  #
#'  for aging worms, we do this in reverse order, which ensures the
#'  correct order of aging is respected
  #

  for(j in (pars$worm_stages-1):1)
  {
  #'  choose the number of male and female worms to age from one stage to the next
  aging_females = rbinom(1, humans$female_worms[j], p)
  aging_males = rbinom(1,humans$male_worms[j], p)

  #' add and subtract the number of worms from the appropriate categories
  humans$female_worms[j+1] = humans$female_worms[j+1]+aging_females
  humans$female_worms[j] = humans$female_worms[j]- aging_females
  humans$male_worms[j+1] =humans$male_worms[j+1]+ aging_males
  humans$male_worms[j] = humans$male_worms[j]-aging_males
  }
  }

  return(humans)


}
  #' Calculate worm pairs in humans ------------------------------------------

  #' function calculate_worm_pairs(humans)
  #'     return min(sum(humans.female_worms), sum(humans.male_worms))
  #' end


 #   calculate_worm_pairs(female_worms, male_worms)

#' calculate how many pairs of worms there are in each human host
#' @export
calculate_worm_pairs <- function(female_worms, male_worms)
{
  worm_pairs = as.integer()
  for (i in 1:length(female_worms)){
    append(worm_pairs, min(sum(female_worms[i]), sum(male_worms[i])))
  }



  return(worm_pairs)
  #return min.(sum.(female_worms), sum.(male_worms))
}



  #' Number of eggs produced -------------------------------------------------

 #' function to calculate the number of eggs produced
#'  this is done by choosing from a negative binomial distribution for each worms,
#'  where the mean and aggregation parameters are calculated as in the
#' "Refined stratified-worm-burden models that incorporate specific biological features
#' of human and snail hosts provide better estimates of Schistosoma diagnosis,
#' transmission, and control" paper
 #' for julia the negative binomial describes the number of failures before
#'  the given number of successes
#'  in a collection of independent Bernoulli trials.
#'  we need to specify a probability of success, and a given number of
#'  successes, which are derived
#'  from the mean and aggregation in the function below
  #

    #' inputs

    #' r - aggregation factor for NB distribution


#'    egg_production!(humans, pars)

#' function to produce eggs for individuals, dependent on how many worms they have
#'        and the max fecundity and density dependent fecundity of the population
#' @export
egg_production <-  function(humans, pars)
{
  male_worms = humans$male_worms
  female_worms = humans$female_worms
  worm_pairs = calculate_worm_pairs(female_worms, male_worms)

  for (i in 1 : length(humans)){
    wp = worm_pairs[i]
    eggs = 0
  }
  #        wp = calculate_worm_pairs(humans[i])

  if (wp > 0){
    wp = max(wp, 1e-10)
    mean_eggs = pars$max_fecundity * wp *
      exp(- pars$density_dependent_fecundity  * wp) * pars$time_step

  }

  #' if we have a positive number of worms, then make calculation,
#'  otherwise the number of eggs is trivially 0 =#
    #         if worm_pairs > 0

    #println(female_worms[i][1])
    #' calculate the mean number of eggs we would expect
    #    mean_eggs = pars.max_fecundity * wp *
    #            exp(- pars.density_dependent_fecundity  * female_worms[i][1])

else{
  mean_eggs = pars$max_fecundity * wp *
    exp(- pars$density_dependent_fecundity  * wp) * pars$time_step

}

  if(pars$egg_production_distribution == "Poisson")
  {
    eggs = rpois(1, mean_eggs)
  }

  else if (pars$egg_production_distribution == "NegBin"){
    NB_r = pars$r * wp

    #' calculate the probability of a success
    p = NB_r/(NB_r+mean_eggs)

    #' choose from NB
    eggs = sample(NegativeBinomial(NB_r,p))
  }

  else{
    print("egg_production_distribution must be 'Poisson' or 'NegBin'")

  }
    # mean_eggs = 0.5*(pars.max_fecundity * M0)/(M0 + M) * (1-exp(-M0*log(2)*M)) * M
  #' calculate the number of successes

  # println(eggs)
  humans[i]$eggs = eggs

  return(humans)

}
  #' egg production increasing


#'    egg_production_increasing!(humans, pars)

#' function to produce eggs for individuals, dependent on how many worms they have
#'        and the max fecundity and density dependent fecundity of the population
#' @export
egg_production_increasing <-  function(humans, pars){

  male_worms = (pars$male_worms)#.(humans)
  female_worms = (pars$female_worms)#.(humans)
  worms = sum(male_worms) + sum(female_worms)
  for (i in 1 : length(humans)){

  #        wp = calculate_worm_pairs(humans[i])
  M = worms[i]

  #' if we have a positive number of worms, then make calculation,
 #' otherwise the number of eggs is trivially 0 =#

    #' calculate the mean number of eggs we would expect


    mean_eggs = max(1E-10, 0.5*(pars$max_fecundity * pars$M0)/(pars$M0 + M) * (1-exp(-pars$M0*log(2)*M)) * M)
  #' calculate the number of successes


  #' calculate the probability of a success


  #' choose from NB
  #eggs = rand(NegativeBinomial(NB_r,p))
  eggs = rpois(1, mean_eggs)
  # println(eggs)
  humans[i]$eggs = eggs
  }
  return(humans)

}

  #' Miracidia production ----------------------------------------------------


  #' function

#'    miracidia_production!(humans)

#' release eggs from individuals into the environment as miracidia. Release is relative to the contact rate with
#' the environment for each individual.
#' @export
miracidia_production <- function(humans){
  released_eggs = 0
  released_eggs =released_eggs + (humans$eggs * humans$relative_contact_rate)

  return(round(released_eggs))

}



#' Death of humans ---------------------------------------------------------
#' @export
death_of_human <- function(humans){

  for (i in length(humans):1){
    if (humans[i]$age > humans[i]$death_age)
    {
      splice(humans, i)
    }


  }

  return(humans)

}
#' Birth of humans


#' add an individual to the population
#' @export
birth_of_human <-  function(humans, pars){

  if ((length(pars$community_probs) != pars$N_communities)){
    print("must provide probabilities for membership of each community")

  }
   else{
     community_selection = 1
   }

  if(pars$N_communities > 1){
    community_selection = cumsum(pars$community_probs)/sum(pars$community_probs)

  }


  pre = rgamma(pars$predis_aggregation, 1/pars$predis_aggregation)
  predisp = sample(pre)[1]
  f_worms = rep(0, pars$worm_stages)
  m_worms = rep(0, pars$worm_stages)


  if (runif(1) > pars$mda_adherence){
    adherence = 0
  }
  else{
    adherence = 1
  }


  if (runif(1) > pars$mda_access){
    access = 0
  }

  else{
    access = 1
  }

  community = filter(community_selection > runif(1))[1]

  # Human(age,gender, predisposition, female_worms, male_worms, eggs, vac_status, age_contact_rate, adherence, access, community)
  death_age = get_death_age(pars)


  append(humans, Human(0, death_age, sample(c(0,1),1), predisp,
                      f_worms, m_worms,
                      0, 0, 0, adherence, access, community,0,0,0,0, "", 0))


  tail(humans$age_contact_rate,1) = pars$contact_rate_by_age_array[1]



  tail(humans$relative_contact_rate,1) = tail(humans$age_contact_rate,1) *  pars$community_contact_rate[tail(humans$community,1)]/
    (max(pars$community_contact_rate) * max(pars$contact_rate_by_age_array))


  if (tail(humans$gender,1) == 0){
    tail(humans$predisposition) = tail(humans$predisposition,1) * pars$female_factor

  }
  else{
    tail(humans$predisposition,1) = tail(humans$predisposition,1) * pars$male_factor
}

  tail(humans$uptake_rate,1) = tail(humans$predisposition,1) * pars$contact_rate * tail(humans$age_contact_rate,1) *
    pars$community_contact_rate[community]

  return(humans)

}

##' Administer drugs


#' function to administer drug to a specific variable (e.g. female_worms or eggs).
#' input the variable, the indices to apply to and the effectiveness of treatment


#'    administer_drug(humans, indices, drug_effectiveness)

#'administer mda drugs to chosen individuals in the population. If they adhere to the drugs, then they
#'reduce male and female worms with a given efficacy alongside removing eggs
#' @export
administer_drug <-  function(humans, indices, drug_effectiveness){

 for (i in 1:length(indices)){
   index = indices[i]
   p = 1 - (drug_effectiveness * humans[index]$adherence)
   humans[index]$female_worms = rbinom(1,humans[index]$female_worms,p)
   humans[index]$male_worms = rbinom(1,humans[index]$male_worms,p)
   #         @inbounds humans[index].human_cercariae = rand.(Binomial.(humans[index].human_cercariae,p))
   #humans[index].eggs = 0
 }

  return(humans)

}
#' Administer vaccine ------------------------------------------------------


#' function to administer drug to a specific variable (e.g. female_worms or eggs).
#' input the variable, the indices to apply to and the effectiveness of treatment

#'    administer_vaccine(humans, indices, vaccine_effectiveness, vaccine_duration)

#'administer vaccine to chosen individuals in the population.
#'reduce male and female worms with a given efficacy alongside removing eggs and
#'adding to their vaccine status signifying that they will have increased immunity for a chosen period of time
#' @export
administer_vaccine <-  function(humans, indices, vaccine_effectiveness, vaccine_duration){


    for (i in 1:length(indices)){
    index = indices[i]
    p = 1 - (vaccine_effectiveness)
    humans[index]$female_worms = rbinom(1,humans[index]$female_worms,p)
    humans[index]$male_worms = rbinom(1, humans[index]$male_worms, p)
    humans[index]$eggs = 0
    humans[index]$vac_status = vaccine_duration * humans[index]$adherence

  }
  return(humans)
}

  #' Mass drug administration ------------------------------------------------


  #' function for mass drug administration
  #' currently there is no correlation between individuals chosen each time

  #'  mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)

#'administer mda in the population. This includes choosing individuals between specified ages,
#'having a certain level of coverage and taking access and adherence into consideration
#' @export
mda <-  function(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender){

  ages = humans$age
  gender = humans$gender
  access = humans$access

  #' find index of people with correct ages for the mda treatment
 # in_gender = in(mda_gender)

  x = filter(((ages <= min_age_mda &ages <max_age_mda) & gender==mda_gender & (access == 1)))


  k = length(x)


  #' randomly permute the indices
  y = sample(x)

  #' only take as many as are indicated by the coverage
  y = y[1:trunc(round(k*mda_coverage))]

  #' update female and male worms, human cercariae and eggs

  humans = administer_drug(humans, y, mda_effectiveness)
  #   else

  #   end

  return(humans)
}


  #' Update MDA --------------------------------------------------------------

  #' function to update the mda information


#    update_mda(mda_info, mda_round)

#' update when the next mda will take place
#' @export
update_mda <-  function(mda_info, mda_round){

  i = min(mda_round + 1, length(mda_info))
  mda_coverage = mda_info[i]$coverage
  min_age_mda =  mda_info[i]$min_age
  max_age_mda =  mda_info[i]$max_age
  mda_effectiveness =  mda_info[i]$effectiveness
  mda_gender = mda_info[i]$gender
  if (mda_round == length(mda_info)){
    next_mda_time = Inf
  }
  else{
    next_mda_time = mda_info[i]$time
  }

  return(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender)

}

#' create MDA --------------------------------------------------------------


#' function to create a set of mda's which will be performed regularly
#' first_mda_time specifies when this will first occur in years,
#'last_mda_time is the final mda in this block
#'regularity is how often to perform the mda in years.
#'specify the proportion of pre SAC, SAC and adults at each of these time points
#'also specify genders for these differect age groups, along with the effectiveness of mda


#    create_mda(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
#            last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness)

#' function to create a set of mda's which will be performed regularly
#' @param first_mda_time specifies when mda will first be administered in years
#' @param last_mda_time is the final mda in this block
#' @param regularity is how often to perform the mda in years
#' @param pre_SAC_prop is the proportion of pre SAC given treatment at each of the time points
#' @param SAC_prop is the proportion of SAC given treatment at each of the time points
#' @param adult_prop is the proportion of adults given treatment at each of the time points
#        also specify genders for these different age groups, along with the effectiveness of mda
#' @export
create_mda <-  function(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
                      last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness){



  mda_info <- as.array()
  mda_time = first_mda_time
  while (mda_time < last_mda_time){
    append(mda_info, mda_information(pre_SAC_prop, 0, 4, pre_SAC_gender, mda_effectiveness, mda_time))
    append(mda_info, mda_information(SAC_prop, 4, 16, SAC_gender, mda_effectiveness, mda_time))
    append(mda_info, mda_information(adult_prop, 16, 110, adult_gender, mda_effectiveness, mda_time))
    mda_time = mda_time+ regularity
  }

  return(mda_info)

}

  #' add vaccination to population -------------------------------------------

  #' function to add vaccination to population

  # function vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
  #     vaccine_gender, vaccine_duration, vaccine_round)

  #     x = []
  #     #= find index of people with correct ages for the mda treatment =#
  #     for i in 1 : length(humans)
  #         if (humans[i].gender in vaccine_gender) & (min_age_vaccine <= humans[i].age <= max_age_vaccine)
  #             push!(x, i)
  #         end
  #     end


  # #=  if this is the first mda round, then treat entirely at random
  #     find how many people are eligible for the treatment  =#

  #     #if mda_round == 0
  #         k = length(x)


  # #= randomly permute the indices =#
  #         y = shuffle(x)

  # #= only take as many as are indicated by the coverage  =#
  #         y = y[1:trunc(Int, round(k*vaccine_coverage))]



  #         humans = administer_vaccine(humans, y, vaccine_effectiveness, vaccine_duration)
  #  #   end

  #     #println("output = ", female_worms, male_worms, human_cercariae, eggs)
  #     return humans
  #     #return x
  # end


  #' update vaccine information ----------------------------------------------


  #' function to update vaccine information =#
  # function update_vaccine(vaccine_info, vaccine_round)

  #     i = min(vaccine_round + 1, size(vaccine_info)[1])
  #     vaccine_coverage = vaccine_info[i].coverage
  #     min_age_vaccine =  vaccine_info[i].min_age
  #     max_age_vaccine =  vaccine_info[i].max_age
  #     vaccine_duration = vaccine_info[i].duration
  #     vaccine_gender = vaccine_info[i].gender
  #     if vaccine_round === size(vaccine_info)[1]
  #         next_vaccine_time = Inf
  #     else
  #         next_vaccine_time = vaccine_info[min(vaccine_round + 1, size(vaccine_info)[1])].time
  #     end
  #     return vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender

  # end


  #' vaccine decay -----------------------------------------------------------


#    vac_decay!(humans, pars)

#' decrease vaccination status for each person by 1 each day
#' @export
vac_decay <-  function(humans, pars)
{

  humans$vac_status = humans$vac_status-pars$time_step

  return(humans)
}



  #' kato_katz eggs ----------------------------------------------------------


#    kato_katz(eggs, gamma_k)
#
#' calculate number of eggs using kato katz method. Gamma_k is a gamma distribution with shape and scale
#' defined by pars.kato_katz_par
#' @export
kato_katz <-  function(eggs, gamma_k){
  gamma1 = runif(gamma_k)[1]
  gamma2 = runif(gamma_k)[1]
  pois1 = rpois(1, 1*gamma1*eggs)
  pois2 = rpois(1, 1*gamma2*eggs)
  return(floor(0.5*(pois1 + pois2)))

}



  #' count number of eggs ----------------------------------------------------

 #  count_eggs(humans)

#' count the total number of eggs in the human population
#' @export
count_eggs <-  function(humans){
  eggs = 0
  eggs = eggs+humans$eggs

  return(eggs)

}



  #' get prevalences ---------------------------------------------------------


#    get_prevalences!(humans, time, pars)

#' calculate the desired prevalences in the human population, and store them in an out struct
#' @export
get_prevalences <-  function(humans, time, pars){


  pop_burden = array()
  sac_burden = array()
  adult_burden = array()
  pop_prev = 0
  sac_prev = 0
  adult_prev = 0
  sac_pop = 0
  adult_pop = 0
  recorded_eggs = c()
  final_ages = c()
  gamma_k  = rgamma(pars$kato_katz_par, 1/kato_katz_par)
  num_humans = length(humans)
  append(final_ages, humans$age)
  # final_eggs = kato_katz(eggs[i], gamma_k)
  # take a sample of the eggs of the individual. For urine sample this may be ~ 1/100 and is defined
  # by the egg_sample_size parameter. for use of kato katz method, this is used for stool samples and hence the sample
  #size is different to urine sample
  sampled_eggs1 = rbinom(1,trunc(humans$eggs/pars$time_step), pars$egg_sample_size)
  sampled_eggs2 = rbinom(1,trunc(humans$eggs/pars$time_step), pars$egg_sample_size)
  sampled_eggs = 0.5 * (sampled_eggs1 + sampled_eggs2)
  final_eggs = sampled_eggs * (1-pars$use_kato_katz) + kato_katz(sampled_eggs, gamma_k) * pars$use_kato_katz

  append(recorded_eggs, final_eggs)
  if (humans$age > 5 & humans$age < 15){
    sac_pop = sac_pop + 1
  }

  if (humans$age > 15){
    adult_pop = adult_pop + 1
  }
  if (final_eggs > pars$heavy_burden_threshold){
    pop_burden[3] = pop_burden[3] + 1
    pop_burden[2] = pop_burden[2] + 1
    pop_burden[1] = pop_burden[1] + 1
    pop_prev = pop_prev + 1
  }

  else if (humans$age < 5  & humans$age < 15){
    sac_burden[3] = sac_burden[3] + 1
    sac_burden[2] = sac_burden[2] + 1
    sac_burden[1] = sac_burden[1] + 1
    sac_prev = sac_prev + 1
  }
  else if (humans$age > 15){
    adult_burden[3] = adult_burden[3] + 1
    adult_burden[2] = adult_burden[2] + 1
    adult_burden[1] = adult_burden[1] + 1
    adult_prev = adult_prev + 1
  }

  else if (final_eggs > 4){
    pop_burden[2] = pop_burden[2] + 1
    pop_burden[1] = pop_burden[1] + 1
    pop_prev = pop_prev + 1
  }

  else if (humans$age > 5 & humans$age < 15){
    sac_burden[2] = sac_burden[2] + 1
    sac_burden[1] = sac_burden[1] + 1
    sac_prev = sac_prev + 1
  }
  else if (humans$age > 15){
    adult_burden[2] = adult_burden[2] + 1
    adult_burden[1] = adult_burden[1] + 1
    adult_prev = adult_prev + 1
  }

  else if (final_eggs > 0)
  {
    pop_burden[1] = pop_burden[1] + 1
    pop_prev = pop_prev + 1
  }

  else if (humans$age > 5 & humans$age < 15){
    sac_burden[1] = sac_burden[1] + 1
    sac_prev = sac_prev +  1
  }

  if (humans$age > 15)
  {
    adult_burden[1] = adult_burden[1] + 1
    adult_prev = adult_prev + 1
  }


  output = out( round(100 * pop_burden/num_humans, digits = 2),
                round(100 *sac_burden/sac_pop, digits = 2),
                round(100 *adult_burden/adult_pop, digits = 2),
                round.(100 *pop_prev / num_humans, digits = 2),
                round(100 *sac_prev / sac_pop, digits = 2),
                round(100 *adult_prev / adult_pop, digits = 2),
                sac_pop, adult_pop,  final_ages, recorded_eggs,
                time)

  return(output)


}
#' save population to file -------------------------------------------------


#' save the enironment variables in a specified file
#' @export
save_population_to_file <- function(filename, humans, miracidia, cercariae, pars){
    data <- tibble("humans"= humans,  "miracidia"=miracidia, "cercariae"=cercariae, "pars"= pars)
    write_csv(data, filename)

  }




#' load_population_from_file(filename)

#' load the environmental variables saved in the specified file
#' @export

load_population_from_file <- function(filename){
  d = read_csv(filename)
  humans = d[, "humans"][1]
  miracidia = d[, "miracidia"][1]
  cercariae = d[, "cercariae"][1]
  pars = d[, "pars"][1]

  return(humans,miracidia,cercariae, pars)

}



#' generate a distribution for ages ----------------------------------------

#' function to generate a distribution for ages based on a specified demography

#' generate population numbers for each age in

#' @export
generate_age_distribution <-  function(pars)
{
  number_per_age = c()
  for (i in 1:(length(pars$spec_ages) * pars$ages_per_index))
  {
    index = trunc( (i-1)/pars$ages_per_index) + 1
    append(number_per_age, pars$spec_ages[index])

  }
  cumsum_spec_ages = cumsum(number_per_age)/sum(number_per_age)
  return(cumsum_spec_ages)

}



#' construct the set of ages -----------------------------------------------

#' function to construct the set of ages, with size N

#' specified_age_distribution(pars)

#' output ages according to a specified age distribution
#' @export
specified_age_distribution <-function(pars)
{
  cumsum_spec_ages = generate_age_distribution(pars)
  ages = c()
  for (i in 1:pars$N)
  {
    x = runif(1)
    k = findall(cumsum_spec_ages > x)[1]
    append(ages, k-1)
  }

  return(ages)

}




  #' update environment to equilibrium ---------------------------------------


  #'  update_env_to_equilibrium(num_time_steps, humans, miracidia, cercariae, pars)

#' update the population for a given length of time. Here we do not age the population or include birth, deaths or interventions.
#' @export
update_env_to_equilibrium <-  function(num_time_steps, humans, miracidia, cercariae, pars){

  sim_time = 0
  record_time = 0
  record <- c()

  for (j in 1:num_time_steps){
    if (sim_time >= record_time){
      a = get_prevalences(humans, sim_time, pars)
      append(record, a)
      record_time = record_time + pars$record_frequency
    }

  }



  sim_time =  sim_time + pars$time_step/365

  humans =   egg_production(humans, pars)

  humans =  worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

  #'  uptake larvae into humans from the environment
  humans <- cercariae_uptake(humans, cercariae, miracidia, pars)[1]
  cercariae <- cercariae_uptake(humans, cercariae, miracidia, pars)[2]
  miracidia <- cercariae_uptake(humans, cercariae, miracidia, pars)[3]
  #'  kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
  #'  kill cercariae in the environment at specified death rate
  cercariae =   cercariae_death(cercariae, miracidia, pars)

  return(humans, miracidia, cercariae, record)
}

  #' update the pop for a given length of time



#'update_env_to_equilibrium_human_larvae(num_time_steps, humans, miracidia, cercariae, pars)

#'update the population for a given length of time. Here we do not age the population or include birth, deaths or interventions and for this
#'function larvae are uptaken from the environment into a larvae category in the humans, rather than immediately becoming worms
#' @export
update_env_to_equilibrium_human_larvae <-  function(num_time_steps, humans, miracidia, cercariae, pars){

  sim_time = 0
  record_time = 0
  record <- c()

  for (j in 1:num_time_steps){
    if (sim_time >= record_time){

    a = get_prevalences(humans, sim_time, pars)
    append(record, a)
    record_time = record_time + pars$record_frequency
  }


  sim_time = sim_time+ pars$time_step/365

  humans =   egg_production(humans, pars)

  humans =  worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

  #' uptake larvae into humans from the environment
  humans <- cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)
  cercariae <- cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)[2]
  miracidia <- cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)[3]

  humans = human_larvae_maturity(humans, pars)
  #' kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
  #' kill cercariae in the environment at specified death rate
  cercariae =   cercariae_death(cercariae, miracidia, pars)


  }
  return(humans, miracidia, cercariae, record)
}


#' update environment to equilibrium ---------------------------------------


#' update_env_to_equilibrium_increasing(num_time_steps, humans, miracidia, cercariae, pars)

#' update the population for a given length of time. Here we do not age the population or include birth, deaths or interventions and for this
#' function larvae are uptaken from the environment immediately to worms and eggs are produced using a monotonically increasing function

#' @export
update_env_to_equilibrium_increasing <- function(num_time_steps, humans, miracidia, cercariae, pars){

  sim_time = 0
  record_time = 0
  record <-c()

  for (j in 1:num_time_steps){


  if (sim_time >= record_time){

  a = get_prevalences(humans, sim_time, pars)
  append(record, a)
  record_time = record_time+ pars$record_frequency

}
  sim_time = sim_time + pars$time_step/365

  humans =   egg_production_increasing(humans, pars)

  humans =  worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

#'  uptake larvae into humans from the environment
  humans = cercariae_uptake(humans, cercariae, miracidia, pars)
  cercariae= cercariae_uptake(humans, cercariae, miracidia, pars)
   miracidia = cercariae_uptake(humans, cercariae, miracidia, pars)

#'  kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
#'  kill cercariae in the environment at specified death rate
  cercariae =   cercariae_death(cercariae, miracidia, pars)

}
  return(humans, miracidia, cercariae, record)


}
  #' update env constant population ------------------------------------------



  #'  update_env_constant_population(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)

#'update the population for a given length of time. Here we include deaths and for each death an individual is immediately born.
 #'   Interventions are included in this function and larvae are immediately uptaken as worms

#' @export
update_env_constant_population <- function(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){


  update_contact_death_rates = 1/5
  sim_time = 0
  record_time = pars$record_frequency
  record = c()
  print_time = 0

  if (length(mda_info)[1] > 0)
  {

  mda_round = 0
  mda_gender = mda_info[1]$gender
  mda_coverage = mda_info[1]$coverage
  min_age_mda =  mda_info[1]$min_age
  max_age_mda =  mda_info[1]$max_age
  mda_effectiveness =  mda_info[1]$effectiveness
  next_mda_time = mda_info[1]$time
  }
  else{
    next_mda_time = Inf
  }


  if (size(vaccine_info)[1] > 0){

  vaccine_round = 0
  vaccine_coverage = vaccine_info[1]$coverage
  vaccine_gender = vaccine_info[1]$gender
  min_age_vaccine =  vaccine_info[1]$min_age
  max_age_vaccine =  vaccine_info[1]$max_age
  next_vaccine_time = vaccine_info[1]$time
  vaccine_duration = vaccine_info[1]$duration
  }
  else{
    next_vaccine_time = Inf
  }


  for (j in 1:num_time_steps){


  if (sim_time >= update_contact_death_rates){
    humans = update_contact_rate(humans,  pars)
    update_contact_death_rates = update_contact_death_rates + (1/5)

  }


  if (sim_time >= record_time){

  a =  get_prevalences(humans, sim_time, pars)
  append(record, a)
  record_time = record_time + pars$record_frequency
}

  sim_time = sim_time+ pars$time_step/365

  humans$age =humans$age + pars$time_step/365

}

  humans =  egg_production(humans, pars)

  humans =    worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

  humans = vac_decay(humans, pars)

  humans =   death_of_human(humans)

  if (length(humans) < pars$N){

  for (k in 1:(pars$N - length(humans))){

  humans = birth_of_human(humans, pars)
  }
  }


  #'  uptake larvae into humans from the environment
  #'  uptake larvae into humans from the environment
  humans = cercariae_uptake(humans, cercariae, miracidia, pars)[1]
  cercariae = cercariae_uptake(humans, cercariae, miracidia, pars)[2]
  miracidia = cercariae_uptake(humans, cercariae, miracidia, pars)[3]


  #' check if we are at a point in time in which an mda is scheduled to take place
  if (sim_time >= next_mda_time){


  #' perform mda
  humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
  #' update information for the next round of mda
  mda_round = mda_round + 1
  mda_coverage =
    update_mda(mda_info, mda_round)[1]
  min_age_mda =
    update_mda(mda_info, mda_round)[2]
  max_age_mda =
    update_mda(mda_info, mda_round)[3]
 mda_effectiveness =
    update_mda(mda_info, mda_round)[4]
 next_mda_time =
   update_mda(mda_info, mda_round)[5]
 mda_gender =
   update_mda(mda_info, mda_round)[6]
  }




  #' check if we are at a point in time in which a vaccine is scheduled to take place
  if (sim_time >= next_vaccine_time){



  #' perform vaccination
  humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                     vaccine_gender, vaccine_duration, vaccine_round)
  #' update information for the next round of vaccination
  vaccine_round = vaccine_round + 1
  vaccine_coverage =
    update_vaccine(vaccine_info, vaccine_round)[1]
  min_age_vaccine =
    update_vaccine(vaccine_info, vaccine_round)[2]
  max_age_vaccine =
    update_vaccine(vaccine_info, vaccine_round)[3]
  next_vaccine_time=
    update_vaccine(vaccine_info, vaccine_round)[4]
  vaccine_gender =
    update_vaccine(vaccine_info, vaccine_round)[5]
  }



  #'  kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
  #'  kill cercariae in the environment at specified death rate
  cercariae =  cercariae_death(cercariae, miracidia, pars)


  return(humans, miracidia, cercariae, record)
}



  #' update env_sonstant population ------------------------------------------



#'    update_env_constant_population_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)

#'update the population for a given length of time. Here we include deaths and for each death an individual is immediately born.
#'Interventions are included in this function and larvae are uptaken as larvae in the humans.
#' @export
update_env_constant_population_human_larvae <-  function(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){

  update_contact_death_rates = 1/5
  sim_time = 0
  record_time = pars$record_frequency
  record = c()
  print_time = 0

  if (size(mda_info)[1] > 0){

  mda_round = 0
  mda_gender = mda_info[1]$gender
  mda_coverage = mda_info[1]$coverage
  min_age_mda =  mda_info[1]$min_age
  max_age_mda =  mda_info[1]$max_age
  mda_effectiveness =  mda_info[1]$effectiveness
  next_mda_time = mda_info[1]$time
  }
  else{
    next_mda_time = Inf
  }



  if (size(vaccine_info)[1] > 0){

  vaccine_round = 0
  vaccine_coverage = vaccine_info[1]$coverage
  vaccine_gender = vaccine_info[1]$gender
  min_age_vaccine =  vaccine_info[1]$min_age
  max_age_vaccine =  vaccine_info[1]$max_age
  next_vaccine_time = vaccine_info[1]$time
  vaccine_duration = vaccine_info[1]$duration
  }
  else{
    next_vaccine_time = Inf
  }


  for (j in 1:num_time_steps){

  if (sim_time >= update_contact_death_rates){
  humans = update_contact_rate(humans,  pars)
  update_contact_death_rates =update_contact_death_rates + 1/5
  }

  if (sim_time >= record_time){

  a = get_prevalences(humans, sim_time, pars)
  append(record, a)
  record_time = record_time+ pars$record_frequency
}

  sim_time = sim_time + pars$time_step/365

  humans$age = humans$age + pars$time_step/365

  }

  humans =  egg_production(humans, pars)

  humans =  worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

  humans = vac_decay(humans, pars)

  humans = death_of_human(humans)

  if (length(humans) < pars$N){


  for (k in 1:(pars$N - length(humans))){
    humans = birth_of_human(humans, pars)
  }
  }


  #'  uptake larvae into humans from the environment
  #'  uptake larvae into humans from the environment

  humans = cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)[1]
  cercariae = cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)[2]
  miracidia = cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)[3]

  humans = human_larvae_maturity(humans, pars)

  #' check if we are at a point in time in which an mda is scheduled to take place
  if (sim_time >= next_mda_time){

  #' perform mda
  humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
  #' update information for the next round of mda
  mda_round =  mda_round + 1
  mda_coverage = update_mda(mda_info, mda_round)[1]
  min_age_mda = update_mda(mda_info, mda_round)[2]
  max_age_mda = update_mda(mda_info, mda_round)[3]
  mda_effectiveness = update_mda(mda_info, mda_round)[4]
  next_mda_time = update_mda(mda_info, mda_round)[5]
  mda_gender = update_mda(mda_info, mda_round)[6]
}


  #' check if we are at a point in time in which a vaccine is scheduled to take place
  if (sim_time >= next_vaccine_time){

  #' perform vaccination
  humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                     vaccine_gender, vaccine_duration, vaccine_round)
  #' update information for the next round of vaccination
  vaccine_round = vaccine_round + 1
  vaccine_coverage =
    update_vaccine(vaccine_info, vaccine_round)[1]
  min_age_vaccine =
    update_vaccine(vaccine_info, vaccine_round)[2]
  max_age_vaccine =
    update_vaccine(vaccine_info, vaccine_round)[3]
  next_vaccine_time =
    update_vaccine(vaccine_info, vaccine_round)[4]
  vaccine_gender =
    update_vaccine(vaccine_info, vaccine_round)[5]
}



  #'  kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
  #'  kill cercariae in the environment at specified death rate
  cercariae =  cercariae_death(cercariae, miracidia, pars)


  return(humans, miracidia, cercariae, record)

}
  #' update env constant population increasing -------------------------------



#'    update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)

#'update the population for a given length of time. Here we include deaths and for each death an individual is immediately born.
#'Interventions are included in this function and larvae are uptaken immediately as worms and egg production follows a monotonically increasing function
#' @export
update_env_constant_population_increasing<- function(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){


  update_contact_death_rates = 1/5
  sim_time = 0
  record_time = pars$record_frequency
  record = c()
  print_time = 0

  if (size(mda_info)[1] > 0)
  {


  mda_round = 0
  mda_gender = mda_info[1]$gender
  mda_coverage = mda_info[1]$coverage
  min_age_mda =  mda_info[1]$min_age
  max_age_mda =  mda_info[1]$max_age
  mda_effectiveness =  mda_info[1]$effectiveness
  next_mda_time = mda_info[1]$time
  }
  else{
    next_mda_time = Inf
  }




  if (length(vaccine_info)[1] > 0){
  vaccine_round = 0
  vaccine_coverage = vaccine_info[1]$coverage
  vaccine_gender = vaccine_info[1]$gender
  min_age_vaccine =  vaccine_info[1]$min_age
  max_age_vaccine =  vaccine_info[1]$max_age
  next_vaccine_time = vaccine_info[1]$time
  vaccine_duration = vaccine_info[1]$duration
  }
  else{
    next_vaccine_time = Inf
  }


  for (j in 1:num_time_steps){

  if (sim_time >= update_contact_death_rates){

  humans = update_contact_rate(humans,  pars)
  update_contact_death_rates = update_contact_death_rates + 1/5
  }

  if (sim_time >= record_time){

  a = get_prevalences(humans, sim_time, pars)
  append(record, a)
  record_time = record_time + pars$record_frequency
  }

  sim_time = sim_time + pars$time_step/365

  humans$age = humans$age + pars$time_step/365

  }

  humans = egg_production_increasing(humans, pars)

  humans = worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

  humans = vac_decay(humans, pars)

  humans = death_of_human(humans)

  if (length(humans) < pars$N){

  for (k in 1:(pars$N - length(humans))){
  humans = birth_of_human(humans, pars)
  }
  }


  #'  uptake larvae into humans from the environment
  #'  uptake larvae into humans from the environment
  humans = cercariae_uptake(humans, cercariae, miracidia, pars)[1]
  cercariae = cercariae_uptake(humans, cercariae, miracidia, pars)[2]
  miracidia = cercariae_uptake(humans, cercariae, miracidia, pars)[3]


  #' check if we are at a point in time in which an mda is scheduled to take place
  if (sim_time >= next_mda_time){

  #' perform mda
  humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
  #' update information for the next round of mda
  mda_round = mda_round + 1

  mda_coverage = update_mda(mda_info, mda_round)[1]
  min_age_mda = update_mda(mda_info, mda_round)[2]
  max_age_mda = update_mda(mda_info, mda_round)[3]
  mda_effectiveness = update_mda(mda_info, mda_round)[4]
  next_mda_time = update_mda(mda_info, mda_round)[5]
  mda_gender = update_mda(mda_info, mda_round)[6]

  }


  #' check if we are at a point in time in which a vaccine is scheduled to take place
  if (sim_time >= next_vaccine_time){


  #' perform vaccination
  humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                     vaccine_gender, vaccine_duration, vaccine_round)
  #' update information for the next round of vaccination
  vaccine_round = vaccine_round + 1

  vaccine_coverage = update_vaccine(vaccine_info, vaccine_round)[1]
  min_age_vaccine = update_vaccine(vaccine_info, vaccine_round)[2]
  max_age_vaccine = update_vaccine(vaccine_info, vaccine_round)[3]
  next_vaccine_time = update_vaccine(vaccine_info, vaccine_round)[4]
  vaccine_gender = update_vaccine(vaccine_info, vaccine_round)[5]

  }



  #'  kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
  #'  kill cercariae in the environment at specified death rate
  cercariae =   cercariae_death(cercariae, miracidia, pars)


  return(humans, miracidia, cercariae, record)

}

  #' update env births and deaths --------------------------------------------



#'    update_env_no_births_deaths(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)

#'update the population for a given length of time. Here we do not include births or deaths and individuals do not age
#'Interventions are included in this function and larvae are uptaken immediately as worms
#' @export
update_env_no_births_deaths <-  function(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){


  update_contact_death_rates = 1/5
  sim_time = 0
  record_time = pars$record_frequency
  record = c()
  print_time = 0

  if (length(mda_info)[1] > 0){
  mda_round = 0
  mda_gender = mda_info[1]$gender
  mda_coverage = mda_info[1]$coverage
  min_age_mda =  mda_info[1]$min_age
  max_age_mda =  mda_info[1]$max_age
  mda_effectiveness =  mda_info[1]$effectiveness
  next_mda_time = mda_info[1]$time
  }
  else{
    next_mda_time = Inf
  }


  if (length(vaccine_info)[1] > 0){

  vaccine_round = 0
  vaccine_coverage = vaccine_info[1]$coverage
  vaccine_gender = vaccine_info[1]$gender
  min_age_vaccine =  vaccine_info[1]$min_age
  max_age_vaccine =  vaccine_info[1]$max_age
  next_vaccine_time = vaccine_info[1]$time
  vaccine_duration = vaccine_info[1]$duration
  }
  else{
    next_vaccine_time = Inf
  }

  for (j in 1:num_time_steps){


  if (sim_time >= update_contact_death_rates){

  humans = update_contact_rate(humans,  pars)
  update_contact_death_rates = update_contact_death_rates + 1/5
  }

  if (sim_time >= record_time){

  a = get_prevalences(humans, sim_time, pars)
  append(record, a)
  record_time = record_time + pars$record_frequency
  }

  sim_time = sim_time +  pars$time_step/365



  humans =   egg_production(humans, pars)

  humans =  worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

  humans = vac_decay(humans, pars)

  #'  uptake larvae into humans from the environment  =#
  #'  uptake larvae into humans from the environment  =#
  humans = cercariae_uptake(humans, cercariae, miracidia, pars)[1]
  cercariae = cercariae_uptake(humans, cercariae, miracidia, pars)[2]
  miracidia = cercariae_uptake(humans, cercariae, miracidia, pars)[3]



  #' check if we are at a point in time in which an mda is scheduled to take place
  if (sim_time >= next_mda_time){


  #' perform mda
  humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
  #' update information for the next round of mda
  mda_round = mda_round + 1
  mda_coverage = update_mda(mda_info, mda_round)[1]
  min_age_mda = update_mda(mda_info, mda_round)[2]
  max_age_mda = update_mda(mda_info, mda_round)[3]
  mda_effectiveness = update_mda(mda_info, mda_round)[4]
  next_mda_time = update_mda(mda_info, mda_round)[5]
  mda_gender = update_mda(mda_info, mda_round)[6]
  }
  }


  #' check if we are at a point in time in which a vaccine is scheduled to take place
  if (sim_time >= next_vaccine_time){


  #' perform vaccination
  humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                     vaccine_gender, vaccine_duration, vaccine_round)
  #' update information for the next round of vaccination
  vaccine_round = vaccine_round + 1
  vaccine_coverage = update_vaccine(vaccine_info, vaccine_round)[1]
  min_age_vaccine = update_vaccine(vaccine_info, vaccine_round)[2]
  max_age_vaccine = update_vaccine(vaccine_info, vaccine_round)[3]
  next_vaccine_time = update_vaccine(vaccine_info, vaccine_round)[4]
  vaccine_gender = update_vaccine(vaccine_info, vaccine_round)[5]
  }



  #'  kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
  #'  kill cercariae in the environment at specified death rate
  cercariae =   cercariae_death(cercariae, miracidia, pars)


  return (humans, miracidia, cercariae, record)
  }


  #' update env_no_births_deaths_human larvae --------------------------------


 #'   update_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)

#'update the population for a given length of time. Here we do not include births or deaths and individuals do not age
#'Interventions are included in this function and larvae are uptaken as larvae in humans
#' @export
update_env_no_births_deaths_human_larvae <-  function (num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){



  update_contact_death_rates = 1/5
  sim_time = 0
  record_time = pars$record_frequency
  record = c()
  print_time = 0

  if (length(mda_info)[1] > 0){

  mda_round = 0
  mda_gender = mda_info[1]$gender
  mda_coverage = mda_info[1]$coverage
  min_age_mda =  mda_info[1]$min_age
  max_age_mda =  mda_info[1]$max_age
  mda_effectiveness =  mda_info[1]$effectiveness
  next_mda_time = mda_info[1]$time
  }
  else{

    next_mda_time = Inf
  }


  if (length(vaccine_info)[1] > 0){

  vaccine_round = 0
  vaccine_coverage = vaccine_info[1]$coverage
  vaccine_gender = vaccine_info[1]$gender
  min_age_vaccine =  vaccine_info[1]$min_age
  max_age_vaccine =  vaccine_info[1]$max_age
  next_vaccine_time = vaccine_info[1]$time
  vaccine_duration = vaccine_info[1]$duration
  }
  else{
    next_vaccine_time = Inf
  }

  for (j in 1:num_time_steps){

  if (sim_time >= update_contact_death_rates){

  humans = update_contact_rate(humans,  pars)
  update_contact_death_rates = update_contact_death_rates + 1/5
  }

  if (sim_time >= record_time){

  a = get_prevalences(humans, sim_time, pars)
  append(record, a)
  record_time = record_time + pars$record_frequency
  }

  sim_time = sim_time + pars$time_step/365



  humans =   egg_production(humans, pars)

  humans =  worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

  humans = vac_decay(humans, pars)

  #'  uptake larvae into humans from the environment
  #' uptake larvae into humans from the environment
  humans = cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)[1]
  cercariae = cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)[2]
  miracidia = cercariae_uptake_with_human_larvae(humans, cercariae, miracidia, pars)[3]

  humans = human_larvae_maturity(humans, pars)


  #' check if we are at a point in time in which an mda is scheduled to take place
  if (sim_time >= next_mda_time){


  #' perform mda
  humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
  #' update information for the next round of mda
  mda_round =  mda_round + 1
  mda_coverage = update_mda(mda_info, mda_round)[1]
  min_age_mda = update_mda(mda_info, mda_round)[2]
  max_age_mda = update_mda(mda_info, mda_round)[3]
  mda_effectiveness = update_mda(mda_info, mda_round)[4]
  next_mda_time = update_mda(mda_info, mda_round)[5]
  mda_gender = update_mda(mda_info, mda_round)[6]


  }


  #' check if we are at a point in time in which a vaccine is scheduled to take place
  if (sim_time >= next_vaccine_time){

  #' perform vaccination
  humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                     vaccine_gender, vaccine_duration, vaccine_round)
  #' update information for the next round of vaccination
  vaccine_round = vaccine_round + 1
  vaccine_coverage = update_vaccine(vaccine_info, vaccine_round)[1]
  min_age_vaccine = update_vaccine(vaccine_info, vaccine_round)[2]
  max_age_vaccine = update_vaccine(vaccine_info, vaccine_round)[3]
  next_vaccine_time = update_vaccine(vaccine_info, vaccine_round)[4]
  vaccine_gender = update_vaccine(vaccine_info, vaccine_round)[5]

}



  #'  kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
  #'  kill cercariae in the environment at specified death rate
  cercariae =   cercariae_death(cercariae, miracidia, pars)

}
  return(humans, miracidia, cercariae, record)
}


  #' update env no. births and deaths increasing -----------------------------


#'    update_env_no_births_deaths_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)

#'update the population for a given length of time. Here we do not include births or deaths and individuals do not age
#'Interventions are included in this function and larvae are uptaken as immediately as worms and egg production is monotonically increasing
#' @export
update_env_no_births_deaths_increasing <-  function(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info){

  update_contact_death_rates = 1/5
  sim_time = 0
  record_time = pars$record_frequency
  record = c()
  print_time = 0

  if (length(mda_info)[1] > 0){

  mda_round = 0
  mda_gender = mda_info[1]$gender
  mda_coverage = mda_info[1]$coverage
  min_age_mda =  mda_info[1]$min_age
  max_age_mda =  mda_info[1]$max_age
  mda_effectiveness =  mda_info[1]$effectiveness
  next_mda_time = mda_info[1]$time
  }
  else{

    next_mda_time = Inf
  }


  if (length(vaccine_info)[1] > 0){

  vaccine_round = 0
  vaccine_coverage = vaccine_info[1]$coverage
  vaccine_gender = vaccine_info[1]$gender
  min_age_vaccine =  vaccine_info[1]$min_age
  max_age_vaccine =  vaccine_info[1]$max_age
  next_vaccine_time = vaccine_info[1]$time
  vaccine_duration = vaccine_info[1]$duration
  }
  else{
    next_vaccine_time = Inf
  }

  for (j in 1:num_time_steps){


  if (sim_time >= update_contact_death_rates){

  humans = update_contact_rate(humans,  pars)
  update_contact_death_rates = update_contact_death_rates + 1/5
  }

  if (sim_time >= record_time){

  a = get_prevalences(humans, sim_time, pars)
  append(record, a)
  record_time = record_time + pars$record_frequency
  }

  sim_time = sim_time + pars$time_step/365



  humans =   egg_production_increasing(humans, pars)

  humans =  worm_maturity(humans, pars)

  append(miracidia, miracidia_production(humans))

  humans = vac_decay(humans, pars)

  #'  uptake larvae into humans from the environment
  #'  uptake larvae into humans from the environment
  humans = cercariae_uptake(humans, cercariae, miracidia, pars)[1]
  cercariae = cercariae_uptake(humans, cercariae, miracidia, pars)[2]
  miracidia = cercariae_uptake(humans, cercariae, miracidia, pars)[3]
  }

  #' check if we are at a point in time in which an mda is scheduled to take place
  if (sim_time >= next_mda_time){


  #' perform mda
  humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
  #' update information for the next round of mda
  mda_round = mda_round + 1

  mda_coverage = update_mda(mda_info, mda_round)[1]
  min_age_mda = update_mda(mda_info, mda_round)[2]
  max_age_mda = update_mda(mda_info, mda_round)[3]
  mda_effectiveness = update_mda(mda_info, mda_round)[4]
  next_mda_time = update_mda(mda_info, mda_round)[5]
  mda_gender = update_mda(mda_info, mda_round)[6]

  }


  #' check if we are at a point in time in which a vaccine is scheduled to take place
  if (sim_time >= next_vaccine_time){



  #' perform vaccination
  humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
                     vaccine_gender, vaccine_duration, vaccine_round)
  #' update information for the next round of vaccination
  vaccine_round = vaccine_round + 1
  vaccine_coverage = update_vaccine(vaccine_info, vaccine_round)
  min_age_vaccine = update_vaccine(vaccine_info, vaccine_round)
  max_age_vaccine = update_vaccine(vaccine_info, vaccine_round)
  next_vaccine_time = update_vaccine(vaccine_info, vaccine_round)
  vaccine_gender = update_vaccine(vaccine_info, vaccine_round)

  }



  #'  kill miracidia in the environment at specified death rate
  miracidia =  miracidia_death(miracidia, pars)
  #'  kill cercariae in the environment at specified death rate
  cercariae =   cercariae_death(cercariae, miracidia, pars)


  return(humans, miracidia, cercariae, record)
  }


  #' Store prevalence and sac prevalence -------------------------------------

  #' when we run multiple simulations, we store them in an array. This function will store the prevalence and sac prevalence


#    collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)

#'collect multiple prevalences within the population and store in appropriate arrays
#' @export
collect_prevs <-  function(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run){

  if (run == 1){

  for (i in 1 : length(record)){

  append(times, record[i]$time)
  append(prev, record[i]$pop_prev)
  append(sac_prev, record[i]$sac_prev)
  append(high_burden, record[i]$population_burden[3])
  append(high_burden_sac, record[i]$sac_burden[3])
  append(adult_prev, record[i]$adult_prev)
  append(high_adult_burden, record[i]$adult_burden[3])
  }}
  else{
    for (i in 1 : length(record)){
  append(prev[i], record[i]$pop_prev)
  append(sac_prev[i], record[i]$sac_prev)
  append(high_burden[i], record[i]$population_burden[3])
  append(high_burden_sac[i], record[i]$sac_burden[3])
  append(adult_prev[i], record[i]$adult_prev)
  append(high_adult_burden[i], record[i]$adult_burden[3])
    }
  }

  return (times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden)
}


  #' repeat simulations ------------------------------------------------------

  #' repeat simulations where we allow mdas and vaccination, but keep the population the same by adding a birth for every death

#    run_repeated_sims_no_population_change_human_larvae(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

#'run multiple simulations where the population is aged, but each death is replaced by a birth and larvae are uptaken by humans as larvae
#' @export
run_repeated_sims_no_population_change_human_larvae <-  function (filename, num_time_steps, mda_info, vaccine_info, num_repeats){



  times = c()
  prev = c()
  sac_prev = c()
  high_burden = c()
  high_burden_sac = c()
  adult_prev = c()
  high_adult_burden = c()


  for (run in 1:num_repeats){



  humans = load_population_from_file(filename)[1]
  miracidia = load_population_from_file(filename)[2]
  cercariae = load_population_from_file(filename)[3]
  pars = load_population_from_file(filename)[4]


  humans =
    update_env_constant_population_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[1]
  miracidia =
    update_env_constant_population_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2]
  cercariae =
    update_env_constant_population_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[3]
  record =
    update_env_constant_population_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[4]



  times = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[1]
  prev = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[2]
  sac_prev = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[3]
  high_burden= collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[4]
  high_burden_sac = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[5]
  adult_prev = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[6]
  high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,  high_burden_sac, adult_prev, high_adult_burden, record, run)[7]


  }

  return(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden)
}


  #' repeated simulations ----------------------------------------------------



#'    run_repeated_sims_no_population_change(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

#'run multiple simulations where aging of the population is not included, larvae are uptaken by humans as worms
#' @export
run_repeated_sims_no_population_change_increasing <-  function (filename, num_time_steps, mda_info, vaccine_info, num_repeats){


  times = c()
  prev = c()
  sac_prev = c()
  high_burden = c()
  high_burden_sac =c()
  adult_prev = c()
  high_adult_burden = c()


  for (run in 1:num_repeats){



  humans = load_population_from_file(filename)[1]
  miracidia = load_population_from_file(filename)[2]
  cercariae= load_population_from_file(filename)[3]
  pars = load_population_from_file(filename)[4]



  humans =
    update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[1]
  miracidia =
    update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2]
  cercariae =
    update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[3]
  record =
    update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[4]


  times = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[1]
  prev = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[2]
  sac_prev = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[3]
  high_burden = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[4]
  high_burden_sac = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[5]
  adult_prev = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[6]
  high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[7]


  }

  return(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden)
}


  #' repeat simulations ------------------------------------------------------


  #' repeat simulations where we allow mdas and vaccination, but keep the population the same by adding a birth for every death

#    run_repeated_sims_no_births_deaths(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

#'run multiple simulations where aging of the population is not included and larvae are uptaken by humans as worms
#' @export
run_repeated_sims_no_births_deaths <-  function (filename, num_time_steps, mda_info, vaccine_info, num_repeats){


  times = c()
  prev = c()
  sac_prev = c()
  high_burden = c()
  high_burden_sac =c()
  adult_prev = c()
  high_adult_burden = c()


  for (run in 1:num_repeats){


  humans = load_population_from_file(filename)[1]
  miracidia = load_population_from_file(filename)[2]
  cercariae = load_population_from_file(filename)[3]
  pars = load_population_from_file(filename)[4]



  humans =
    update_env_no_births_deaths(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[1]
  miracidia =
    update_env_no_births_deaths(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2]
  cercariae =
    update_env_no_births_deaths(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[3]
  record =
    update_env_no_births_deaths(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[4]



  times = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[1]
  prev = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[2]
  sac_prev = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[3]
  high_burden = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[4]
  high_burden_sac = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[5]
  adult_prev = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[6]
  high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[7]

  }

  return(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden)

}



  #' repeat simulations ------------------------------------------------------


#    run_repeated_sims_no_births_deaths_human_larvae(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

#' run multiple simulations where aging of the population is not included and larvae are uptaken by humans as larvae
#' @export
run_repeated_sims_no_births_deaths_human_larvae <-  function(filename, num_time_steps, mda_info, vaccine_info, num_repeats){

  times = c()
  prev = c()
  sac_prev = c()
  high_burden = c()
  high_burden_sac =c()
  adult_prev = c()
  high_adult_burden = c()


  for (run in 1:num_repeats){


  humans = load_population_from_file(filename)[1]
  miracidia = load_population_from_file(filename)[2]
  cercariae = load_population_from_file(filename)[3]
  pars = load_population_from_file(filename)[4]

  humans = update_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[1]
  miracidia = update_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2]
  cercariae = update_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[3]
  record = update_env_no_births_deaths_human_larvae(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[4]



  times = collect_prevs(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden, record, run)[1]
  prev = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[2]
  sac_prev = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[3]
  high_burden = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[4]
  high_burden_sac = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[5]
  adult_prev = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[6]
  high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[7]

  }

  return(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden)
}


  #' repeat simulations ------------------------------------------------------

  #' repeat simulations where we allow mdas and vaccination, but keep the population the same by adding a birth for every death


 #   run_repeated_sims_no_births_deaths_human_larvae(filename, num_time_steps, mda_info, vaccine_info, num_repeats)

#'run multiple simulations where aging of the population is not included and larvae are uptaken by humans as worms, and egg production is
#'monotonically increasing
#' @export
run_repeated_sims_no_births_deaths_increasing <-  function (filename, num_time_steps, mda_info, vaccine_info, num_repeats){


  times = c()
  prev = c()
  sac_prev = c()
  high_burden = c()
  high_burden_sac =c()
  adult_prev = c()
  high_adult_burden = c()


  for (run in 1:num_repeats){



  humans = load_population_from_file(filename)[1]
  miracidia = load_population_from_file(filename)[2]
  cercariae = load_population_from_file(filename)[3]
  pars = load_population_from_file(filename)[4]


  humans =
    update_env_no_births_deaths_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[1]
  miracidia =
    update_env_no_births_deaths_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[2]
  cercariae =
    update_env_no_births_deaths_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[3]
  record =
    update_env_no_births_deaths_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)[4]




  times = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[1]
  prev = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[2]
  sac_prev = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[3]
  high_burden = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[4]
  high_burden_sac = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[5]
  adult_prev = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[6]
  high_adult_burden = collect_prevs(times, prev, sac_prev, high_burden,high_burden_sac, adult_prev, high_adult_burden, record, run)[7]



  }

  return(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden)
}


  # et al -------------------------------------------------------------------

  #
  # function update_env(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)
  #
  #
  #     update_contact_death_rates = 1/5
  #     sim_time = 0
  #     record_time = pars.record_frequency
  #     record = []
  #     print_time = 0
  #
  #     if size(mda_info)[1] > 0
  #         mda_round = 0
  #         mda_gender = mda_info[1].gender
  #         mda_coverage = mda_info[1].coverage
  #         min_age_mda =  mda_info[1].min_age
  #         max_age_mda =  mda_info[1].max_age
  #         mda_effectiveness =  mda_info[1].effectiveness
  #         next_mda_time = mda_info[1].time
  #     else
  #         next_mda_time = Inf
  #     end
  #
  #
  #     if size(vaccine_info)[1] > 0
  #         vaccine_round = 0
  #         vaccine_coverage = vaccine_info[1].coverage
  #         vaccine_gender = vaccine_info[1].gender
  #         min_age_vaccine =  vaccine_info[1].min_age
  #         max_age_vaccine =  vaccine_info[1].max_age
  #         next_vaccine_time = vaccine_info[1].time
  #         vaccine_duration = vaccine_info[1].duration
  #     else
  #         next_vaccine_time = Inf
  #     end
  #
  #     for j in 1:num_time_steps
  #
  #         if sim_time >= update_contact_death_rates
  #             humans = update_death_rate(humans, death_rate_per_time_step)
  #             humans = update_contact_rate(humans,  contact_rates_by_age)
  #             update_contact_death_rates += 1/5
  #         end
  #
  #         if sim_time >= record_time
  #             a = get_prevalences!(humans, sim_time, pars)
  #             push!(record, a)
  #             record_time += pars.record_frequency
  #         end
  #
  #         sim_time += pars.time_step/365
  #
  #         for h in humans
  #             h.age += pars.time_step/365
  #
  #         end
  # #=  mature larvae within humans  =#
  #         humans = human_cercariae_maturity(humans, time_step)
  #
  #         humans = egg_production(humans, max_fecundity, r, density_dependent_fecundity)
  #
  #         humans = worm_maturity(humans, worm_stages, average_worm_lifespan, time_step)
  #
  #         miracidia = miracidia_production(humans, miracidia)
  #
  #         humans = vac_decay!(humans)
  #
  # # larvae = larvae_production(d, larvae)
  #         humans = death_of_human(humans, time_step)
  #
  # #=  uptake larvae into humans from the environment  =#
  #         humans, cercariae, miracidia = cercariae_uptake(humans, cercariae, miracidia, pars)
  #
  # #= check if we are at a point in time in which an mda is scheduled to take place =#
  #         if sim_time >= next_mda_time
  #
  # #= perform mda =#
  #             humans = mda(humans, mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, mda_gender)
  # #= update information for the next round of mda =#
  #             mda_round += 1
  #             mda_coverage, min_age_mda, max_age_mda, mda_effectiveness, next_mda_time, mda_gender =
  #                             update_mda(mda_info, mda_round)
  #
  #         end
  #
  #
  # #= check if we are at a point in time in which a vaccine is scheduled to take place =#
  #         if sim_time >= next_vaccine_time
  #
  # #= perform vaccination =#
  #             humans = vaccinate(humans, vaccine_coverage, min_age_vaccine, max_age_vaccine, vaccine_effectiveness,
  #                 vaccine_gender, vaccine_duration, vaccine_round)
  # #= update information for the next round of vaccination =#
  #             vaccine_round += 1
  #             vaccine_coverage, min_age_vaccine, max_age_vaccine, next_vaccine_time, vaccine_gender =
  #                                         update_vaccine(vaccine_info, vaccine_round)
  #         end
  #
  #
  #
  # #=  kill miracidia in the environment at specified death rate =#
  #         miracidia = miracidia_death(miracidia, pars)
  #
  # #=  kill cercariae in the environment at specified death rate =#
  #         cercariae = cercariae_death(cercariae, env_cercariae_death_rate, time_step)
  #
  #
  #         l = rand(Binomial(length(humans), birth_rate))[1]
  #         if l > 0
  #             for i in 1:l
  #                 humans = birth_of_human(humans, pars)
  #             end
  #         end
  #     end
  #     return humans, miracidia, cercariae, record
  # end



