#' @title Global Probability object
ProbabilityObj <- R6::R6Class("ProbabilityObj", #nolint
                              public = list(
                                #' @description
                                #' Gets the a probability vector
                                #'
                                #' @return a named vector of integers (factors) corresponding to the
                                #' probability of each entity.
                                get_values = function() {
                                  return(private$probability_vect)
                                }
                              ),
                              private = list(
                                probability_vect = NULL
                              )
)
#' @title
#' Primary probability object
#'
#' @description
#' Object used as third layer for the primary signature
#' @export
PrimaryProbabilityObj <- R6::R6Class("PrimaryProbabilityObj", #nolint
                                     inherit = ProbabilityObj,
                                     public = list(
                                       #' @description
                                       #' Initialize `PrimaryProbabilityObj` object.
                                       #'
                                       #' @param nb_ent : Integer representing the number of entities
                                       #' (genes or transcripts)
                                       #' @param entity_id (optional) : character vector with the entity names.
                                       #' Default `NA`.
                                       initialize = function(nb_ent, entity_id = NULL) {
                                         private$probability_vect <- runif(nb_ent)
                                         names(private$probability_vect) <- entity_id
                                       }
                                     )
)
#' @title
#' Secondary probability object
#'
#' @description
#' Object used as third layer for the connected signature
#' @export
SecondaryProbabilityObj <- R6::R6Class("SecondaryProbabilityObj", #nolint
                                       inherit = ProbabilityObj,
                                       public = list(
                                         #' @description
                                         #' Initialize `SecondaryProbabilityObj` object.
                                         #'
                                         #' @param prim_perc_obj : `PrimaryProbabilityObj` used as base for the
                                         #' simulation of the connected object.
                                         #' @param proba_transition : choice of the probability transition
                                         #' "independent", "cop" or "deterministic".
                                         #' @param theta (optional): parameter of the copula for the "cop"
                                         #' transition. If:
                                         #' - cop = "Frank": theta [-Inf, Inf]; theta → Inf comonotonicity ;
                                         #' theta → -Inf countermonotonicity; theta → 0 independence copula.
                                         #' - cop = "Plackett": theta (0, Inf) \ {1}; theta → Inf comonotonicity ;
                                         #' theta → 0 countermonotonicity; theta → 1 independence copula.
                                         #' - cop = "Gauss": theta [-1,1] ; theta → 1 comonotonicity ;
                                         #' theta → -1 countermonotonicity; theta → 0 independence copula
                                         #' @param copula (optional) : 'Gauss','Plackett' or 'Frank' copula
                                         #' used when the probability transition is "cop".
                                         #' @param nbins (optional): number of bins used to generate the
                                         #' probability layer by the probability "cop" transition.
                                         #' @param ncpus (optional) : number of cpus for parallel calculation used in
                                         #' 'cop' probability transition.
                                         initialize = function(prim_perc_obj, proba_transition,
                                                               copula, theta, nbins, ncpus) {
                                           prim_perc_vect <- prim_perc_obj$get_values()
                                           nb_ent <- length(prim_perc_vect)
                                           private$probability_vect <- switch(as.character(proba_transition),
                                                                              independent = runif(nb_ent),
                                                                              deterministic = prim_perc_vect,
                                                                              cop = private$cop_transition(prim_perc_vect = prim_perc_vect,
                                                                                                           copula = copula, theta = theta, nbins = nbins, ncpus = ncpus)
                                           )
                                           names(private$probability_vect) <- names(prim_perc_vect)
                                         }
                                       ),
                                       private = list(
                                         cop_transition = function(prim_perc_vect,
                                                                   copula, theta, nbins, ncpus) {
                                           if (copula == "Gauss" & (theta > 1 | theta < (-1))) {
                                             stop("The definition interval for the theta parameter of the Gauss
          copula is [-1,1] ")
                                           }else if (copula == "Plackett" & theta < 0) {
                                             stop("The definition interval for the theta parameter of the Plackett
          copula is (0, Inf) \ {1} ")
                                           }
                                           cop_fn <- switch(copula,
                                                            Gauss = copula::normalCopula(theta),
                                                            Frank = copula::frankCopula(theta),
                                                            Plackett = copula::plackettCopula(theta))
                                           support <- seq(0, 1, length.out = nbins)
                                           probability_vect <- parallel::mclapply(prim_perc_vect, function(x) {
                                             probs <- copula::dCopula(copula = cop_fn,
                                                                      u = matrix(c(rep(x, nbins), support), ncol = 2))
                                             sample(x = support, size = 1, prob = probs)
                                           }, mc.cores = ncpus)
                                           return(unlist(probability_vect))
                                         }
                                       )
)
