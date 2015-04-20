\name{CombIncrease_next}
\alias{CombIncrease_next}
\alias{print.CombIncrease_next}
\title{Combination determination with logistic model}
\description{
\code{CombIncrease_next} is used to determine the next or recommended combination in a phase I combination clinical trial using the design proposed by Riviere et al. entitled "A Bayesian dose-finding design for drug combination clinical trials based on the logistic model".
}
\usage{
CombIncrease_next(ndose_a1, ndose_a2, target, target_min, target_max, 
prior_tox_a1, prior_tox_a2, final, pat_incl, dose_adm1, dose_adm2, tite=FALSE, 
toxicity, time_full=0, time_tox=0, time_follow=0, c_e=0.85, c_d=0.45, 
c_stop=0.95, n_min)
}
\arguments{ 
  \item{ndose_a1}{Number of dose levels for agent 1.}
  \item{ndose_a2}{Number of dose levels for agent 2.}
  \item{target}{Toxicity (probability) target (for dose allocation).}
  \item{target_min}{Minimum of the targeted toxicity interval (for dose recommendation).}
  \item{target_max}{Maximum of the targeted toxicity interval (for dose recommendation).}
  \item{prior_tox_a1}{A vector of initial guesses of toxicity probabilities associated with the doses of agent 1. Must be of length \code{ndose_a1}.}
  \item{prior_tox_a2}{A vector of initial guesses of toxicity probabilities associated with the doses of agent 2. Must be of length \code{ndose_a2}.}
  \item{final}{A boolean with value TRUE if the trial is finished and the recommended combination for further phases should be given, or FALSE (default value) if the combination determination is performed for the next cohort of patients.}
  \item{pat_incl}{Current number of patients included.}
  \item{dose_adm1}{A vector indicating the dose levels of agents 1 administered to each patient included in the trial. Must be of length \code{pat_incl}.}
  \item{dose_adm2}{A vector indicating the dose levels of agents 2 administered to each patient included in the trial. Must be of length \code{pat_incl}.}
  \item{tite}{A boolean indicating if the toxicity is considered as a time-to-event outcome (TRUE), or as a binary outcome (default value FALSE).}
  \item{toxicity}{A vector of observed toxicities (DLTs) for each patient included in the trial. Must be of length \code{pat_incl}. This argument is used/required only if tite=FALSE.}
  \item{time_full}{Full follow-up time window. This argument is used only if tite=TRUE.}  
  \item{time_tox}{A vector of times-to-toxicity for each patient included in the trial. If no toxicity was observed for a patient, must be filled with +Inf. Must be of length \code{pat_incl}. This argument is used/required only if tite=TRUE.}
  \item{time_follow}{A vector of follow-up times for each patient included in the trial. Must be of length \code{pat_incl}. This argument is used/required only if tite=TRUE.} 
  \item{c_e}{Probability threshold for dose-escalation. The default value is set at 0.85.}
  \item{c_d}{Probability threshold for dose-deescalation. The default value is set at 0.45.}  
  \item{c_stop}{Probability threshold for early trial termination. The default value is set at 0.95.}
  \item{n_min}{Minimum number of patients to be included before possible early trial termination.}
}
\value{
An object of class "CombIncrease_next" is returned, consisting of determination of the next combination and estimations. Objects generated by \code{CombIncrease_next} contain at least the following components: 

\item{n_pat_comb}{Number of patients per combination.}
\item{n_tox_comb}{Number of observed toxicities per combination.}
\item{pi}{Estimated toxicity probabilities (if the start-up ended).}
\item{ptox_inf}{Estimated probabilities that the toxicity probability is inferior to \code{target} (if the start-up ended).}
\item{ptox_inf_targ}{Estimated probabilities of underdosing, i.e. to be inferior to \code{target_min} (if the start-up ended).}
\item{ptox_targ}{Estimated probabilities to be in the targeted interval [\code{target_min},\code{target_max}] (if the start-up ended).}
\item{ptox_sup_targ}{Estimated probabilities of overdosing, i.e. to be superior to \code{target_max} (if the start-up ended).}
\item{startup_in}{Start-up phase is ended or not.}
\item{(cdose1, cdose2)}{NEXT RECOMMENDED COMBINATION.}
\item{cohort}{Cohort size.}
\item{pat_incl}{Number of patients included.}
\item{target}{Toxicity target.}
\item{[target_min, target_max]}{Targeted toxicity interval.}
\item{prior_tox_a1}{Prior toxicity probabilities for agent 1.}
\item{prior_tox_a2}{Prior toxicity probabilities for agent 2.}
\item{n_min}{Minimum number of cohorts to stop the trial.}
\item{c_e}{Escalation threshold.}
\item{c_d}{Deescalation threshold.}
\item{c_stop}{Stopping threshold.}
\item{tite}{Type of outcome for toxicity (time-to-event or binary).}
\item{time_full}{If toxicity is a time-to-event, full follow-up time is also reminded.}
}
\references{
Riviere, M-K., Yuan, Y., Dubois, F., and Zohar, S. (2014). A Bayesian dose-finding design for drug combination clinical trials based on the logistic model. Pharmaceutical Statistics.
}
\author{Jacques-Henri Jourdan and Marie-Karelle Riviere-Jourdan \email{eldamjh@gmail.com}}
\seealso{
  \code{\link{CombIncrease_sim}}.
}
\examples{
prior_a1 = c(0.12, 0.2, 0.3, 0.4, 0.5)
prior_a2 = c(0.2, 0.3, 0.4)
toxicity1 = c(0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1)
dose1 = c(1,1,1,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4)
dose2 = c(1,1,1,2,2,2,3,3,3,2,2,2,1,1,1,1,1,1)
t_tox = c(rep(+Inf,8),2.9,+Inf,4.6,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,5.2)
follow = c(rep(6,15), 4.9, 3.1, 1.3)

next1 = CombIncrease_next(ndose_a1=5, ndose_a2=3, target=0.3, target_min=0.20, 
target_max=0.40, prior_tox_a1=prior_a1, prior_tox_a2=prior_a2, final=FALSE, 
pat_incl=18, dose_adm1=dose1, dose_adm2=dose2, tite=FALSE, toxicity=toxicity1, 
n_min=6)

next1

next2 =CombIncrease_next(ndose_a1=5, ndose_a2=3, target=0.30, target_min=0.20, 
target_max=0.40, prior_tox_a1=prior_a1, prior_tox_a2=prior_a2, final=FALSE, 
pat_incl=18, dose_adm1=dose1, dose_adm2=dose2, tite=TRUE, time_full=6, 
time_tox=t_tox, time_follow=follow, n_min=6)

next2
}