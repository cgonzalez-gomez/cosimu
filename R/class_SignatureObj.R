#' @include class_ModalityObj.R class_SubmodalityObj.R class_ProbabilityObj.R
#' @import R6
NULL

#' @title Global Signature object
SignatureObj <- R6::R6Class("SignatureObj", #nolint
                            public = list(
                              #' @description
                              #' Gets the modality object
                              #'
                              #' @return `ModalityObj` object
                              get_modality_obj = function() {
                                return(private$mod_obj)
                              },
                              #' @description
                              #' Gets the sub-modality object
                              #'
                              #' @return `SubmodalityObj` object
                              get_submodality_obj = function() {
                                return(private$submod_obj)
                              },
                              #' @description
                              #' Gets the probability object
                              #'
                              #' @return `SubmodalityObj` object
                              get_probability_obj = function() {
                                return(private$perc_obj)
                              },
                              #' @description
                              #' Gets the log2 fold change vector
                              #'
                              #' @return log2 fold change vector
                              get_lfc_vect = function() {
                                return(private$lfc_vect)
                              },
                              #' @description
                              #' Gets the base expression
                              #'
                              #' @return base_expression vector
                              get_base_expression = function() {
                                return(private$base_expression)
                              }
                            ),
                            private = list(
                              mod_obj = NULL,
                              submod_obj = NULL,
                              perc_obj = NULL,
                              lfc_vect = NULL,
                              p_up = NULL,
                              p_down = NULL,
                              submod_tr = NULL,
                              perc_tr = NULL,
                              qf_list = NULL,
                              base_expression = NULL,
                              # @description
                              # Generates the fold-change values for each entity.
                              #
                              # @param qf_list nested list containing all the quantile
                              # functions associated to each sub-modality in each modality
                              # @param perc_vect
                              # @param submod_vect
                              # @param mod_vect
                              # @param nb_ent
                              # @return fold change vector.
                              lfc_simulation = function(qf_list, perc_vect, submod_vect, mod_vect,
                                                        nb_ent, entity_id) {
                                lfc_vect <- purrr::map_dbl(seq(1, nb_ent), function(i) {
                                  m <- mod_vect[i]
                                  sm <- submod_vect[i]
                                  p <- perc_vect[i]
                                  qf <- qf_list[[m]][[sm]]
                                  return(qf(p))
                                })
                                names(lfc_vect) <- entity_id
                                return(lfc_vect)
                              }
                            )
)
#' @title
#' Primary signature object
#'
#' @description
#' Object composed by three layers:
#' `PrimaryModalityObj`,`PrimarySubmodalityObj`
#' and `PrimaryProbabilityObj`.
#' @export
PrimarySignatureObj <- R6::R6Class("PrimarySignatureObj", #nolint
                                   inherit = SignatureObj,
                                   public = list(
                                     #' @description
                                     #' Initialize `PrimarySignatureObj` object.
                                     #'
                                     #' @param nb_ent : Integer representing the number of entities
                                     #' (genes or transcripts).
                                     #' @param p_up : Proportion of up-regulated entities.
                                     #' @param p_down : Proportion of down-regulated entities.
                                     #' @param prop_sm_up : Vector of proportions of the
                                     #' sub-modalities of the up-regulated entities.
                                     #' @param prop_sm_down : Vector of proportions of the
                                     #' sub-modalities of the down-regulated entities.
                                     #' @param qf_vect_up : List of quantile functions associated
                                     #' to each sub-modality of the up-regulated modality.
                                     #' @param qf_vect_nr : List of quantile functions associated
                                     #' to each sub-modality of the non-deregulated modality.
                                     #' @param qf_vect_down : List of quantile functions associated
                                     #' to each sub-modality of the down-deregulated modality.
                                     #' @param entity_id (optional) : Character vector with the entity names.
                                     #' @param base_expression (optional) : Boolean set to `TRUE` if their is a
                                     #' need to correct the base expression. Default `FALSE`.
                                     #' @param initial_base (optional) : Double vector with the normalized base
                                     #' expression values to be corrected and control simulated bias.
                                     #' @param up_means (optional) : Vector of means of the up sub-modalities use
                                     #' if `base_expression = TRUE`, to apply the base expression correction.
                                     #' @param nr_means (optional) : Vector with the mean of the non-deregulated
                                     #' use if `base_expression = TRUE`, to apply the base expression correction.
                                     #' @param down_means (optional) : Vector of means of the down sub-modalities use
                                     #' if `base_expression = TRUE`, to apply the base expression correction.
                                     #' Default `NA`.
                                     initialize = function(nb_ent, p_up, p_down, prop_sm_up,
                                                           prop_sm_down, qf_vect_up, qf_vect_nr, qf_vect_down,
                                                           entity_id = NULL, base_expression = FALSE, initial_base = NA, up_means = NA,
                                                           nr_means = NA, down_means = NA) {
                                       private$mod_obj <- PrimaryModalityObj$new(nb_ent = nb_ent,
                                                                                 p_up = p_up, p_down = p_down, entity_id = entity_id)
                                       private$submod_obj <- PrimarySubmodalityObj$new(
                                         mod_obj = private$mod_obj,
                                         prop_sm_up = prop_sm_up,
                                         prop_sm_down = prop_sm_down, entity_id = entity_id)
                                       private$perc_obj <- PrimaryProbabilityObj$new(nb_ent = nb_ent,
                                                                                     entity_id = entity_id)
                                       private$qf_list <- list(qf_vect_up, qf_vect_nr, qf_vect_down)
                                       private$lfc_vect <- private$lfc_simulation(qf_list = private$qf_list,
                                                                                  perc_vect = private$perc_obj$get_values(),
                                                                                  submod_vect = private$submod_obj$get_values(),
                                                                                  mod_vect = private$mod_obj$get_values(),
                                                                                  nb_ent = nb_ent,
                                                                                  entity_id = entity_id)
                                       if(base_expression){
                                         stopifnot(length(initial_base) == nb_ent)
                                         sm_vect <- private$submod_obj$get_values()
                                         m_vect <- private$mod_obj$get_values()
                                         list_means <- list(up_means,nr_means,down_means)
                                         private$base_expression <- purrr::map_dbl(seq(1, nb_ent),
                                                                                   function(i) {
                                                                                     return(2 * initial_base[i] /
                                                                                              (1 + 2^list_means[[m_vect[i]]][sm_vect[i]]))
                                                                                   }
                                         )
                                       }
                                     }
                                   )
)
#' @title
#' Secondary signature object
#'
#' @description
#' Object composed by three layers:
#' `SecondaryModalityObj`,`SecondarySubmodalityObj`
#' and `SecondaryProbabilityObj`.
#' @export
SecondarySignatureObj <- R6::R6Class("SecondarySignatureObj", #nolint
                                     inherit = SignatureObj,
                                     public = list(
                                       #' @description
                                       #' Initialize `SecondarySignatureObj` object.
                                       #'
                                       #' @param prim_sig_obj : `PrimarySignatureObj` used as base for the
                                       #' simulation of the connected object.
                                       #' @param mod_transition : Modality transition "default"
                                       #' or "customized".
                                       #' @param submod_transition : Sub-modality transition
                                       #' "cop", "independent", "stochastic" or "deterministic".
                                       #' @param proba_transition : Probability transition
                                       #' "independent", "cop" or "deterministic".
                                       #' @param connectivity_score (optional) : Double between -1 and 1
                                       #' representing the "connectivity" between the two perturbations.
                                       #' Used to generate the transition matrix if the modality transition is not
                                       #' "customized".
                                       #' @param nr_noise (optional) : Double between 0 and 0.5 representing the
                                       #' noise induced by the non-deregulated entities during the "default"
                                       #' transition. Default `NA`.
                                       #' @param transition_mat (optional) : transition matrix used if
                                       #' `mod_transition` = "customized". Default `NA`.
                                       #' @param prop_sm_up (optional) : vector of proportions of the
                                       #' sub-modalities of the up-regulated entities. Default `NA`.
                                       #' @param prop_sm_down (optional) : vector of proportions of the
                                       #' sub-modalities of the down-deregulated entities. Default `NA`.
                                       #' @param qf_vect_up (optional) : list of quantile functions associated
                                       #' to each sub-modality of the up-regulated modality. Default `NA`.
                                       #' @param qf_vect_nr (optional) : list of quantile functions associated
                                       #' to each sub-modality of the non-deregulated modality. Default `NA`.
                                       #' @param qf_vect_down (optional) : list of quantile functions associated
                                       #' to each sub-modality of the down-deregulated modality. Default `NA`.
                                       #' @param entity_id (optional) : character vector with the entity names.
                                       #' Default `NULL`.
                                       #' @param copula_submod (optional) : 'Gauss','Plackett' or 'Frank' copula
                                       #' used when the sub-modality transition is "cop". Default `NA`.
                                       #' @param rho_submod (optional) : correlation to be reached by the
                                       #' sub-modality "cop" transition.
                                       #' @param eps_submod (optional) : Error between the expected correlation
                                       #' (`rho_submod`) and the one estimated with the copula function.
                                       #' Default `NA`.
                                       #' @param optim_method_submod (optional) : Optimization method used
                                       #' during the sub-modality "cop" transition. See `optim` for more
                                       #' details. Default "Brent".
                                       #' @param copula_prob (optional) : 'Gauss','Plackett' or 'Frank' copula
                                       #' used when the probability transition is "cop". Default `Frank`.
                                       #' @param theta_prob : parameter of the copula for the "cop" transition. If:
                                       #' - cop = "Frank": theta [-Inf, Inf]; theta → Inf comonotonicity ;
                                       #' theta → -Inf countermonotonicity; theta → 0 independence copula.
                                       #' - cop = "Plackett": theta (0, Inf) \ {1}; theta → Inf comonotonicity ;
                                       #' theta → 0 countermonotonicity; theta → 1 independence copula.
                                       #' - cop = "Gauss": theta [-1,1] ; theta → 1 comonotonicity ;
                                       #' theta → -1 countermonotonicity; theta → 0 independence copula
                                       #' Default 10.
                                       #' @param nbins_prob (optional): number of bins used to generate the
                                       #' probability layer by the probability "cop" transition. Default 1e3.
                                       #' @param base_expression (optional) : Boolean set to `TRUE` if their is a
                                       #' need to correct the base expression. Default `FALSE`.
                                       #' @param up_means (optional) : vector of means of the up sub-modalities use
                                       #' if `base_expression = TRUE`, to apply the base expression correction.
                                       #' @param nr_means (optional) : vector with the mean of the non-deregulated
                                       #' use if `base_expression = TRUE`, to apply the base expression correction.
                                       #' @param down_means (optional) : vector of means of the down sub-modalities
                                       #' use if `base_expression = TRUE`, to apply the base expression correction.
                                       #' @param ncpus (optional) : number of cpus for parallel calculation used in
                                       #' 'cop' probability transition. Default 1.
                                       initialize = function(prim_sig_obj,
                                                             mod_transition, submod_transition,
                                                             proba_transition, qf_vect_up = NA, qf_vect_nr = NA, qf_vect_down = NA,
                                                             connectivity_score = NA, transition_mat = NA,
                                                             nr_noise = NA, prop_sm_up = NA, prop_sm_down = NA,
                                                             rho_submod = NA, copula_submod = NA, eps_submod = NA,
                                                             optim_method_submod = "Brent", copula_prob = "Frank",
                                                             theta_prob = 10, nbins_prob = 1e3,
                                                             base_expression = FALSE,
                                                             up_means = NA, nr_means = NA, down_means = NA, ncpus = 1) {
                                         private$mod_obj <- SecondaryModalityObj$new(
                                           prim_mod_obj = prim_sig_obj$get_modality_obj(),
                                           connectivity_score = connectivity_score,
                                           transition_mat = transition_mat,
                                           mod_transition = mod_transition, nr_noise = nr_noise)
                                         nb_ent <- private$mod_obj$get_nb_ent()
                                         private$submod_obj <- SecondarySubmodalityObj$new(
                                           mod_obj = private$mod_obj,
                                           prim_mod_obj = prim_sig_obj$get_modality_obj(),
                                           prim_submod_obj = prim_sig_obj$get_submodality_obj(),
                                           submod_transition = submod_transition, prop_sm_up = prop_sm_up,
                                           prop_sm_down = prop_sm_down,
                                           rho = rho_submod, copula = copula_submod, eps = eps_submod,
                                           optim_method = optim_method_submod
                                         )
                                         private$perc_obj <- SecondaryProbabilityObj$new(
                                           prim_perc_obj = prim_sig_obj$get_probability_obj(),
                                           proba_transition = proba_transition,
                                           copula = copula_prob, theta = theta_prob, nbins = nbins_prob,
                                           ncpus = ncpus)
                                         private$qf_list <- list(qf_vect_up, qf_vect_nr, qf_vect_down)
                                         private$lfc_vect <- private$lfc_simulation(qf_list = private$qf_list,
                                                                                    perc_vect = private$perc_obj$get_values(),
                                                                                    submod_vect = private$submod_obj$get_values(),
                                                                                    mod_vect = private$mod_obj$get_values(),
                                                                                    nb_ent = nb_ent,
                                                                                    entity_id = names(prim_sig_obj$get_lfc_vect()))
                                         
                                         if (base_expression) {
                                           sm_vect <- private$submod_obj$get_values()
                                           m_vect <- private$mod_obj$get_values()
                                           prim_m_vect <- prim_sig_obj$get_modality_obj()$get_values()
                                           prim_sm_vect <- prim_sig_obj$get_submodality_obj()$get_values()
                                           list_means <- list(up_means, nr_means, down_means)
                                           prim_base_exp <- prim_sig_obj$get_base_expression()
                                           private$base_expression <- purrr::map_dbl(seq(1, nb_ent),
                                                                                     function(i) {
                                                                                       m <- m_vect[i]
                                                                                       prim_m <- prim_m_vect[i]
                                                                                       if(m == prim_m) {
                                                                                         return(prim_base_exp[i])
                                                                                       }else{
                                                                                         return(2 * prim_base_exp[i] /
                                                                                                  (1 + 2^list_means[[m_vect[i]]][sm_vect[i]]))
                                                                                       }
                                                                                     })
                                         }
                                       },
                                       #' @description
                                       #' Get the vector with the expected and the effective connectivity scores.
                                       #' For all transitions different to "cop" they should be identical.
                                       #'
                                       #' @return list with connectivities
                                       get_connectivities = function() {
                                         return(private$mod_obj$get_connectivities())
                                       }
                                     )
)
