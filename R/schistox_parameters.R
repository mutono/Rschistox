## schistox

### parameters file

# set working directory

## import required files
#source("schistox_functions.R")

# define parameters
N <- 1000 # human population size
time_step<- 20 # length of time step (in days)
N_communities <- 1 #number of communities in the population sharing the same environmental source
community_probs <- 1.0 # probability of being in each community
community_contact_rate <- 1.0 # contact rate with the environment for each of the commununity
density_dependent_fecundity <- 0.0007 # for S. mansoni [Toor et al JID paper SI] #decrease in egg production per worm due to high density of worms
#density_dependent_fecundity <- 0.0006 # for S. haematobium [Toor et al JID paper SI]
average_worm_lifespan <- 5.7 # years for S. mansoni [Toor et al JID paper SI] # average expectancy of a worm
#average_worm_lifespan = 4 # years for S. haematobium [Toor et al JID paper SI]
max_age <- 100 #maximum age of individual
initial_worms <- 10 # initial no. of worms
# miracidia_maturity = 21 # for S. haemotobium
miracidia_maturity <- 24 # for S. mansoni # no of days after which miracidias will mature to cercariae
initial_miracidia <- time_step*5000*N # initial no. of miracidia in the environment
initial_miracidia_days <- round(miracidia_maturity/time_step) # no.of days miracidia will age into cercariae larvae
init_env_cercariae <- time_step*5000*N # initial no of cercaria in the environment
worm_stages <- 1 # number of stages in the worm. Having 1 stage will result to a Gamma distrinution
egg_production_distribution <- "NegBin" # Distribution for egg production, either "Poisson" or "NegBin"
# contact_rate <- 0.1 #global contact rate for the uptake of larvae from the environment
input_contact_rates <- c(0.01, 1.2, 1, 0.02) # input contact rates
age_contact_rates <- input_contact_rates/sum(input_contact_rates) # contact rate for the uptake of larvae from the environment for the chosen age groups
input_ages <- c(4, 9, 15, 100) # input ages for contructing contact array
ages_for_contacts <- input_ages  # age groups for specifying contact rates
contact_rate_by_age_array <- rep(0,times=max_age+1) #array holding contact rate for each age
mda_adherence <- 1 #proportion of people who adhere to the MDA
mda_access <- 1 # proportion of people who have access to the MDA
last_uptake <- 0 # last uptake of MDA
egg_multiplier <- 100 #??
sd_decrease <- 1 #??
female_factor <- 1 # factor for altering the contact rate for females, if we choose to have gender-specific behavior which affects contact rate
male_factor <- 1 # factor for altering the contact rate for males, if we choose to have gender-specific behavior which affects contact rate
birth_rate <- 28*time_step/(1000*365) # rate of birth of humans
human_cercariae_prop <- 1 # proportion of cercariae which are able to infect humans
cercariae_survival <- 0.05 # proportion of cercariae that survive from one time point to the next
miracidia_survival <-0.05 # proportion of miracidia that survive from one time point to the next
death_prob_by_age <- c(0.0656, 0.0093, 0.003, 0.0023, 0.0027, 0.0038, 0.0044, 0.0048,
                      0.0053, 0.0065, 0.0088, 0.0106, 0.0144, 0.021, 0.0333, 0.0529,
                      0.0851, 0.1366, 0.2183, 0.2998 , 0.3698, 1) # probability of dying each year, specified by age
ages_for_death <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                   65, 70, 75, 80, 85, 90, 95, 100, 110) # age ranges for death probailities
r <- 0.03 # aggregation parameter for negative binomially distributed egg production
vaccine_effectiveness <- 0.86 # efficacy of a vaccine if one is used
drug_effectiveness <- 0.86 # efficacy of a drug given during MDA
spec_ages <- c(7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
               2636, 1970, 1468, 1166, 943, 718, 455, 244) #number of individuals by age group
ages_per_index <- 5 #how many different ages we include in the spec_ages parameter
record_frequency <- 1/24 #how often we should record the prevalence in the population during simulation
use_kato_katz <- 0 #if 1, use Kato-Katz for egg count, if 0, do not use KK
kato_katz_par <- 0.87 #parameter for Gamma distribution if KK is used
scenario <- "moderate adult" #can be one of "low adult", "moderate adult" or high adult"

## main parameters that we change
predis_aggregation <- 0.24 # 0.24 for high prev settings; 0.04 for low prev settings # From "The design of schistosomiasis monitoring and evaluation programmes:
#The importance of collecting adult data to inform treatment strategies for Schistosoma mansoni". aggregation for predisposition of individuals to uptake laevale. This is chosen from a Gamma distribution with mean 1 for each individual and set for life. If high, the aggregation is low, meaning individuals have roughly the same predisposition. If low, larvae become concentrated in a few individuals.
max_fecundity <- 20 # expected no. of eggs from a single worm
max_fec_contact_rate_product <- 2 #product of max fecundity and the contact rate in the population. Setting this to a desired value is often a good way to ensure that the epidemic stays within a reasonable range, as when the max fecundity increases, if the contact rate doesn't decrease appropriately, then the behaviour of the outbreak can be unrealistically difficult to control.
contact_rate <- max_fec_contact_rate_product / max_fecundity #global contact rate for the uptake of larvae from the environment
M0 <- 20 #if a particular formula of egg production is used, this parameter is required and is a proxy for mean worm burden
rate_acquired_immunity <- 0 # rate at which immunity will be acquired for individuals. This will be multiplied by the cumulative nymber of worms people have had during their life to decide the level of immunity acquired
human_larvae_maturity_time <- 30 # length of time (in days) after which a cercariae uptake by a human will mature into a worm
egg_sample_size <- 1/100 #the proportion of eggs which are sampled from each individual every time we check their burden (between 0 and 1). 1= all eggs in the person are sampled. Typical value for a urine sample may be ~1/100
heavy_burden_threshold <- 50 #number of eggs at which an individual is said to have a heavy infection



pars <- Parameters(N, time_step, N_communities, community_probs,
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

## update parameters
pars <- make_age_contact_rate_array(pars, scenario, input_ages, input_contact_rates)
