# @include class_ModalityObj.R
NULL

#' @title Global sub-modality object
#' Parent Class for the suubmodality objects, defines the common parameteres
#' and their get functions.
#' Note : NR modality is supposed to have only one submodality.
SubmodalityObj <- R6::R6Class("SubmodalityObj", #nolint
                              public = list(
                                #' @description
                                #' Gets the list of the  sub-modalities proportions of a specific modality
                                #' (or each modality)
                                #'
                                #' @param mod (optional) : Modality for which the  sub-modalities proportions
                                #' should be returned. Default 'NA' returns the vector for all
                                #' the modalities.
                                #' @return  sub-modalities proportions vector or list of all the sub-modalities
                                #' propotions vectors.
                                get_prop_list = function(mod = NA) {
                                  if (is.na(mod)) {
                                    return(private$submod_prop_by_mod)
                                  }else{
                                    return(private$submod_prop_by_mod[[mod]])
                                  }
                                },
                                #' @description
                                #' Gets the a  sub-modalities vector
                                #'
                                #' @return a named vector of integers (factors) corresponding to the
                                #' sub-modality of each entity.
                                get_values = function() {
                                  return(private$submod_vect)
                                },
                                #' @description
                                #' Gets the number of  sub-modalities in a specific modality
                                #' (or each modality)
                                #'
                                #' @param mod (optional) : Modality for which the number of sub-modalities
                                #' should be returned. Default 'NA' returns the value for all
                                #' the modalities.
                                #' @return an integer equal to the number of entities or a vector of values.
                                get_nb_submod_by_mod = function(mod = NA) {
                                  if (is.na(mod)) {
                                    return(private$nb_submod_by_mod)
                                  }else{
                                    return(private$nb_submod_by_mod[[mod]])
                                  }
                                }
                              ),
                              private = list(
                                submod_prop_by_mod = NULL,
                                submod_vect = NULL,
                                nb_submod_by_mod = NULL
                              )
)
#' @title
#' Primary sub-modality object
#'
#' @description
#' Object used as second layer for the independent ("non connected") signature
#' @export
PrimarySubmodalityObj <- R6::R6Class("PrimarySubmodalityObj", #nolint
                                     inherit = SubmodalityObj,
                                     public = list(
                                       #' @description
                                       #' Initialize `PrimarySubmodalityObj` object.
                                       #'
                                       #' @param mod_obj : `ModalityObjPrimary` corresponding to the first
                                       #' layer of the independent signature.
                                       #' @param prop_sm_up : Vector of proportions of the
                                       #'  sub-modalities of the up regulated entities.
                                       #' @param prop_sm_down : Vector of proportions of the
                                       #'  sub-modalities of the down regulated entities.
                                       #' @param entity_id (optional) : Character vector with the entity names.
                                       #' Default `NULL`.
                                       initialize = function(mod_obj, prop_sm_up, prop_sm_down,
                                                             entity_id = NULL) {
                                         #check proportion sum
                                         if ((sum(prop_sm_up) != 1 & !any(is.na(prop_sm_up)))|
                                             (sum(prop_sm_down) != 1 & !any(is.na(prop_sm_down)))) {
                                           stop("The sum of one of your proportions vectors isn't equal to 1.")
                                         }
                                         private$submod_prop_by_mod <- list(prop_sm_up, 1, prop_sm_down)
                                         private$nb_submod_by_mod <- list(length(prop_sm_up),
                                                                          1, length(prop_sm_down))
                                         private$submod_vect <- private$submod_simulation(
                                           mod_vect = mod_obj$get_values(),
                                           probas = private$submod_prop_by_mod,
                                           nb_submod_by_mod = private$nb_submod_by_mod)
                                         names(private$submod_vect) <- entity_id
                                       }
                                     ),
                                     private = list(
                                       submod_simulation = function(mod_vect, probas, nb_submod_by_mod) {
                                         submod_vect <- vapply(mod_vect, function(x) {
                                           switch(x,
                                                  sample(x = seq(1, nb_submod_by_mod[[x]]),
                                                         size = 1, prob = probas[[x]]),
                                                  as.integer(1),
                                                  sample(x = seq(1, nb_submod_by_mod[[x]]),
                                                         size = 1, prob = probas[[x]])
                                           )
                                         }, numeric(1))
                                         return(submod_vect)
                                       }
                                     )
)
#' @title
#' Secondary sub-modality object
#'
#' @description
#' Object used as second layer for the connected signature
#' @export
SecondarySubmodalityObj <- R6::R6Class("SecondarySubmodalityObj", #nolint
                                       inherit = SubmodalityObj,
                                       public = list(
                                         #' @description
                                         #' Initialize `SecondarySubmodalityObj` object.
                                         #'
                                         #' @param mod_obj : `ModalityObjSecondary` corresponding to the first
                                         #' layer of the connected signature.
                                         #' @param prim_mod_obj : `ModalityObjPrimary` used as base for the
                                         #' simulation of the connected object.
                                         #' @param prim_submod_obj : `PrimarySubmodalityObj` used as base for the
                                         #' simulation of the connected object.
                                         #' @param submod_transition : Choice of the sub-modality transition
                                         #' "cop", "independent", "stochastic" or "deterministic".
                                         #' @param prop_sm_up (optional) : Vector of proportions of the
                                         #'  sub-modalities of the up regulated entities. Default `NA`.
                                         #' @param prop_sm_down (optional) : Vector of proportions of the
                                         #'  sub-modalities of the down regulated entities. Default `NA`.
                                         #' @param copula (optional) : 'Gauss','Plackett' or 'Frank' copula
                                         #' used when the sub-modality transition is "cop". Default "Frank".
                                         #' @param rho (optional) : Correlation to be reached by the
                                         #' sub-modality "cop" transition.
                                         #' @param eps (optional) : Error between the expected correlation
                                         #' (`rho`) and the one estimated with the copula function.
                                         #' Default "1e-3.
                                         #' @param optim_method (optional) : Optimization method used
                                         #' during the sub-modality "cop" transition. See `optim` for more
                                         #' details. Default "Brent".
                                         initialize = function(mod_obj, prim_mod_obj, prim_submod_obj,
                                                               submod_transition, prop_sm_up = NA, prop_sm_down = NA,
                                                               rho = NA, copula = "Frank", eps = 1e-3, optim_method = "Brent") {
                                           private$submod_transition <- submod_transition
                                           private$rho <- rho
                                           if (any(is.na(c(prop_sm_up, prop_sm_down)))) {
                                             private$submod_prop_by_mod <-
                                               prim_submod_obj$get_prop_list()
                                             private$nb_submod_by_mod <-
                                               prim_submod_obj$get_nb_submod_by_mod()
                                           }else{
                                             private$submod_prop_by_mod <- list(prop_sm_up, 1,
                                                                                prop_sm_down)
                                             private$nb_submod_by_mod <- list(length(prop_sm_up), 1,
                                                                              length(prop_sm_down))
                                           }
                                           prim_probas <- prim_submod_obj$get_prop_list()
                                           prim_nb_submod_hy_mod <- prim_submod_obj$get_nb_submod_by_mod()
                                           prim_submod_vect <- prim_submod_obj$get_values()
                                           private$submod_vect <- private$submod_simulation(
                                             mod_vect = mod_obj$get_values(),
                                             prim_mod_vect = prim_mod_obj$get_values(),
                                             prim_submod_vect = prim_submod_vect,
                                             nb_ent = mod_obj$get_nb_ent(), probas = private$submod_prop_by_mod,
                                             nb_submod_by_mod = private$nb_submod_by_mod,
                                             submod_transition = submod_transition,
                                             prim_nb_submod_hy_mod = prim_nb_submod_hy_mod,
                                             rho = rho, prim_probas = prim_probas,
                                             copula = copula, eps = eps, optim_method = optim_method)
                                           names(private$submod_vect) <- names(prim_submod_vect)
                                         },
                                         #' @description
                                         #' Get the optimization results. See `optim` method.
                                         #'
                                         #' @return list with optimization results
                                         get_optim_res = function() {
                                           return(private$optim_res)
                                         },
                                         #' @description
                                         #' Get the vector with the expected and the effective correlation (rho)
                                         #' used for "cop" transition.
                                         #'
                                         #' @return list with two elements expected rho and effective one
                                         get_rho = function() {
                                           return(private$rho)
                                         },
                                         #' @description
                                         #' Get a list of the transition matrix between all the sub-modalities
                                         #' of two modalities.
                                         #'
                                         #' @return list of transition matrix
                                         get_list_transition_mat = function() {
                                           return(private$list_transition_mat)
                                         }
                                       ),
                                       private = list(
                                         submod_transition = NULL,
                                         rho = NULL,
                                         optim_res = NULL,
                                         list_transition_mat = NULL,
                                         # @description
                                         # Generates a nested list containing the sub-modality information of all
                                         # entities (in all the modalities) for the second perturbation
                                         # by calling
                                         # `mapped_submod_simulation` for the "stochastic" and "deterministic"
                                         # transitions.
                                         #
                                         # @param probas List of sub-modalities probabilities of each modalitie.
                                         # @param nb_ent Number of entities.
                                         # @param prim_mod_vect Vector containing the modality of each
                                         # entity (gene/transcript) for the first perturbation.
                                         # @param mod_vect Vector containing the modality of each entity
                                         # (gene/transcript) for the second perturbation.
                                         # @param prim_submod_vect List of vectors containing the sub-modality of each
                                         # entity (gene/transcript) for the first perturbation.
                                         # entities between  sub-modalities for the first perturbation.
                                         # @param nb_submod_by_mod Vector containing the number of sub-modalities
                                         # in each modality for the second perturbation.
                                         # @param submod_transition Transition function
                                         # @return sub-modality vector
                                         submod_simulation = function(mod_vect,
                                                                      prim_mod_vect, prim_submod_vect, nb_ent,
                                                                      probas, nb_submod_by_mod, submod_transition, rho,
                                                                      prim_probas, copula, eps, optim_method, prim_nb_submod_hy_mod) {
                                           submod_vect <- rep(1, nb_ent)
                                           no_nr <- which(mod_vect != 2)
                                           if (submod_transition == "independent") {
                                             submod_vect[no_nr] <- vapply(mod_vect[no_nr], function(x){
                                               sample(x = seq(1, nb_submod_by_mod[[x]]),
                                                      size = 1, prob = probas[[x]])
                                             }, numeric(1))
                                           }else if (submod_transition == "cop") {
                                             list_transition_mat <- purrr::map2(rep(seq(1, 3), each = 2),
                                                                                rep(c(1, 3), 3), function(i, j) {
                                                                                  if (i == 2) return(matrix(probas[[j]], nrow = 1))
                                                                                  return(private$copula_tr(margin1 = prim_probas[[i]],
                                                                                                           copula = copula,
                                                                                                           margin2 = probas[[j]], rho = rho, eps = eps,
                                                                                                           optim_method = optim_method))
                                                                                }
                                             )
                                             asso_mat <- matrix(seq(1, 6), ncol = 2, byrow = T)
                                             asso_mat <- cbind(asso_mat[, 1], NA, asso_mat[, 2])
                                             submod_vect[no_nr] <- vapply(no_nr, function(i) {
                                               mod <- mod_vect[i]
                                               prim_mod <- prim_mod_vect[i]
                                               prim_submod <- prim_submod_vect[i]
                                               n <- nb_submod_by_mod[[mod]]
                                               sample(seq(1, n), size = 1,
                                                      prob = list_transition_mat[[asso_mat[prim_mod, mod]]][prim_submod, ])
                                             }, numeric(1))
                                             private$list_transition_mat <- list_transition_mat
                                           }else if (submod_transition == "deterministic") {
                                             submod_vect[no_nr] <- vapply(no_nr, function(i) {
                                               if (prim_mod_vect[i] != 2) {
                                                 return(prim_submod_vect[i])
                                               }else{
                                                 m <- mod_vect[i]
                                                 sample(x = seq(1, nb_submod_by_mod[[m]]),
                                                        size = 1, prob = probas[[m]])
                                               }
                                             }, numeric(1))
                                             
                                           }else if (submod_transition == "stochastic") {
                                             list_transition_mat <- purrr::map2(rep(seq(1, 3), each = 2),
                                                                                rep(c(1, 3), 3), function(i, j) {
                                                                                  if (i == 2) return(matrix(probas[[j]], nrow = 1))
                                                                                  prim_n <- prim_nb_submod_hy_mod[[i]]
                                                                                  n <- nb_submod_by_mod[[j]]
                                                                                  if (n != prim_n) {
                                                                                    warning("The number of  sub-modalities in the vectors are not
              the same, this could cause some troubles during the
              stochastic transition")
                                                                                  }
                                                                                  mat <- suppressMessages(t(purrr::map_dfc(seq(1, prim_n),
                                                                                                                           ~dbinom(x = seq(1, n), size = n, prob = .x / n))))
                                                                                  return(mat)
                                                                                })
                                             asso_mat <- matrix(seq(1, 6), ncol = 2, byrow = T)
                                             asso_mat <- cbind(asso_mat[, 1], NA, asso_mat[, 2])
                                             submod_vect[no_nr] <- vapply(no_nr, function(i) {
                                               mod <- mod_vect[i]
                                               prim_mod <- prim_mod_vect[i]
                                               prim_submod <- prim_submod_vect[i]
                                               n <- nb_submod_by_mod[[mod]]
                                               sample(seq(1, n), size = 1,
                                                      prob = list_transition_mat[[asso_mat[prim_mod, mod]]][prim_submod, ])
                                             }, numeric(1))
                                             private$list_transition_mat <- list_transition_mat
                                           } else {
                                             stop("Invalide sub-modality transition")
                                           }
                                           return(submod_vect)
                                         },
                                         copula_tr = function(margin1, copula, margin2,
                                                              rho, eps, optim_method) {
                                           distr1 <- cumsum(margin1[-length(margin1)])
                                           distr2 <- cumsum(margin2[-length(margin2)])
                                           support1 <- seq(1, length(margin1))
                                           support2 <- seq(1, length(margin2))
                                           check <- GenOrd::corrcheck(list(distr1, distr2),
                                                                      support = list(support1, support2))
                                           rho_min <- check[[1]][1, 2]
                                           rho_max <- check[[2]][1, 2]
                                           if (rho > rho_max | rho < rho_min) {
                                             warning(paste0("The correlation (rho=", rho, ") is not
            inside the attainable interval [", rho_min, ",", rho_max, "]
            rho will be replace by an extreme value"))
                                             rho <- (rho > rho_max) * rho_max + (rho < rho_min) * rho_min
                                           }
                                           cor2theta_obj <- private$cor2theta(margin1 = margin1, margin2 = margin2,
                                                                              support1 = support1, support2 = support2, copula = copula,
                                                                              rho = rho, eps = eps, optim_method = optim_method)
                                           theta <- cor2theta_obj$theta
                                           private$rho <- c(private$rho, cor2theta_obj$rho)
                                           joint_probas <- private$theta2probas(margin1 = margin1,
                                                                                margin2 = margin2, support1 = support1, support2 = support2,
                                                                                copula = copula, param = theta)
                                           joint_probas <- matrix(joint_probas, byrow = TRUE,
                                                                  ncol = length(margin2))
                                           probas <- joint_probas / rowSums(joint_probas)
                                           return(probas)
                                         },
                                         # @description
                                         # Estimation of the copula parameter
                                         # Function finding the value of the copula parameter inducing the target
                                         # linear correlation rho between two discrete variables with assigned
                                         # marginal distributions. Source https://tinyurl.com/ASTB-D-19-00215
                                         # @param margin1 Probabilities of the modalities 'UP', 'NR' and 'DOWN'
                                         # for the first perturbation.
                                         # @param margin2 Probabilities of the modalities 'UP', 'NR' and 'DOWN'
                                         # for the second perturbation.
                                         # @param support1 Ranks of the modalities for the first perturbation.
                                         # @param support2 Ranks of the modalities for the second perturbation.
                                         # @param copula Copula function "Gauss", "Plackett" or "Frank".
                                         # @param rho Target correlation.
                                         # @param eps Error between the expected connectivity score and the
                                         # correlation score estimated with the copula function.
                                         # @param optim_method see `optim`.
                                         cor2theta = function(margin1, margin2, support1, support2, copula,
                                                              rho, eps, optim_method) {
                                           thetaval <- switch(copula,
                                                              # Theorically limits for Gauss are -1 and 1
                                                              Gauss = c(-0.9999999, rho, 0, rho, 0.9999999),
                                                              Frank = c(-100, -1, 0, 1, 100),
                                                              Plackett = c(0, 0.5, 1, 2, 100))
                                           res_optim <- optim(par = thetaval[2 + 2 * ifelse(rho > 0, 1, 0)],
                                                              lower = thetaval[1], upper = thetaval[5], method = optim_method,
                                                              control = list("abstol" = eps),
                                                              fn = function(thetat) {
                                                                rhot <- private$theta2cor(margin1 = margin1, margin2 = margin2,
                                                                                          support1 = support1, support2 = support2, copula = copula,
                                                                                          param = thetat)
                                                                return(abs(rho - rhot))
                                                              }
                                           )
                                           private$optim_res <- c(private$optim_res, res_optim)
                                           # control of the convergence
                                           if (res_optim$convergence != 0) {
                                             # !!!! warning or stop for the submod ?
                                             stop(paste0("Convergence wasn't reach. The `optim`",
                                                         " convergence code was ", res_optim$convergence,
                                                         "; and the message: ", res_optim$message))
                                           }
                                           return(list(theta = res_optim$par,
                                                       rho = private$theta2cor(margin1 = margin1, margin2 = margin2,
                                                                               support1 = support1, support2 = support2, copula = copula,
                                                                               param = res_optim$par)
                                           )
                                           )
                                         },
                                         # Auxiliary function. Source https://tinyurl.com/ASTB-D-19-00215
                                         #
                                         # @param copula Copula function "Gauss", "Plackett" or "Frank".
                                         # @param param Copula parameter.
                                         # @import copula
                                         copfun = function(param, copula) {
                                           switch(copula,
                                                  Gauss = copula::normalCopula(param),
                                                  Frank = copula::frankCopula(param),
                                                  Plackett = copula::plackettCopula(param))
                                         },
                                         # @description
                                         # Function computing the linear correlation between two discrete variables
                                         # with assigned probability mass functions contained in margin1 and margin2
                                         # and finite supports in support1 and support2 linked to a joint distribution
                                         # through:
                                         # \itemize{
                                         #  \item the Gauss copula with parameter -1<=rho<=1
                                         #  \item the Frank copula with parameter kappa (real)
                                         #  \item the Plackett copula with parameter theta>0
                                         # }
                                         # Source https://tinyurl.com/ASTB-D-19-00215
                                         #
                                         # @param margin1 Probabilities of the modalities 'UP', 'NR' and 'DOWN' for the
                                         # first perturbation.
                                         # @param margin2 Probabilities of the modalities 'UP', 'NR' and 'DOWN' for the
                                         # second perturbation.
                                         # @param support1 Ranks of the modalities for the first perturbation.
                                         # @param support2 Ranks of the modalities for the second perturbation.
                                         # @param copula Copula function "Gauss", "Plackett" or "Frank".
                                         # @param param Copula parameter.
                                         # @import copula
                                         theta2cor = function(margin1, margin2, support1, support2, copula, param) {
                                           # select the copula function
                                           cop <- private$copfun(param = param, copula = copula)
                                           # computing the cumulative distribution functions
                                           f1 <- c(0, cumsum(margin1))
                                           f2 <- c(0, cumsum(margin2))
                                           # computing the mixed moment
                                           s <- sum(
                                             purrr::map2_dbl(
                                               rep(seq(1, length(support1)), each = length(support2)),
                                               rep(seq(1, length(support2)), times = length(support1)),
                                               function(i, j) {
                                                 support1[i] * support2[j] * (
                                                   copula::pCopula(c(f1[i + 1], f2[j + 1]), cop) +
                                                     copula::pCopula(c(f1[i], f2[j]), cop) -
                                                     copula::pCopula(c(f1[i], f2[j + 1]), cop) -
                                                     copula::pCopula(c(f1[i + 1], f2[j]), cop))
                                               }
                                             )
                                           )
                                           # computing the marginal means
                                           m1 <- sum(margin1 * support1)
                                           m2 <- sum(margin2 * support2)
                                           # and variances
                                           v1 <- sum(margin1 * support1^2) - m1^2
                                           v2 <- sum(margin2 * support2^2) - m2^2
                                           # computing the linear correlation
                                           corf <- (s - m1 * m2) / sqrt(v1 * v2)
                                           return(corf)
                                         },
                                         
                                         # @description
                                         # Estimation of the linear correlation between two discrete variables
                                         # with assigned marginal distributions using the copula function and
                                         # its parameter.
                                         #
                                         # @param margin1 Probabilities of the modalities 'UP', 'NR' and
                                         # 'DOWN' for the first perturbation.
                                         # @param margin2 Probabilities of the modalities 'UP', 'NR' and
                                         # 'DOWN' for the second perturbation.
                                         # @param support1 Ranks of the modalities for the first perturbation.
                                         # @param support2 Ranks of the modalities for the second perturbation.
                                         # @param copula Copula function "Gauss", "Plackett" or "Frank".
                                         # @param param Copula parameter
                                         # @import copula
                                         theta2probas = function(margin1, margin2, support1, support2,
                                                                 copula, param) {
                                           cop <- private$copfun(param = param, copula = copula)
                                           p_joint <- c()
                                           f1 <- c(0, cumsum(margin1))
                                           f2 <- c(0, cumsum(margin2))
                                           p_joint <- purrr::map2_dbl(
                                             rep(seq(1, length(support1)), each = length(support2)),
                                             rep(seq(1, length(support2)), times = length(support1)),
                                             function(i, j) {
                                               copula::pCopula(c(f1[i + 1], f2[j + 1]), cop) +
                                                 copula::pCopula(c(f1[i], f2[j]), cop) -
                                                 copula::pCopula(c(f1[i], f2[j + 1]), cop) -
                                                 copula::pCopula(c(f1[i + 1], f2[j]), cop)
                                             }
                                           )
                                           if (any(p_joint < 0)) {
                                             if (any(abs(0 - p_joint[p_joint < 0]) > 1e-15)) {
                                               stop("At least one probability was negative : ",
                                                    paste0(as.character(p_joint), sep = "\n"))
                                             }else{
                                               warning("At least one probability was rounded to 0 in a
          sub-modality transition matrix: ",
                                                       paste0(p_joint[p_joint < 0]))
                                               p_joint[p_joint < 0] <- 0
                                             }
                                           }
                                           return(p_joint)
                                         }
                                       )
)
