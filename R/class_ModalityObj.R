#' @title Global Modality object
#' @description
#' Parent Class for the modality objects, defines the common parameteres
#' and their get functions.
ModalityObj <- R6::R6Class("ModalityObj", #nolint
                           public = list(
                             #' @description
                             #' Gets the number of entities in the object
                             #'
                             #' @return an integer equal to the number of entities
                             get_nb_ent = function() {
                               return(private$nb_ent)
                             },
                             #' @description
                             #' Gets the proportion of up regulated entities
                             #'
                             #' @return a double equal to the proportion of up regulates
                             #' entities
                             get_p_up = function() {
                               return(private$p_up)
                             },
                             #' @description
                             #' Gets the proportion of down regulated entities
                             #'
                             #' @return a double equal to the proportion of down regulates
                             #' entities
                             get_p_down = function() {
                               return(private$p_down)
                             },
                             #' @description
                             #' Gets the modalities vector
                             #'
                             #' @return a named vector of integers (factors) corresponding to the
                             #' modality of each entity.
                             get_values = function() {
                               return(private$mod_vect)
                             }
                           ),
                           private = list(
                             nb_ent = NULL,
                             p_up = NULL,
                             p_down = NULL,
                             mod_vect = NULL
                           )
)

#' @title
#' Primary modality object
#' @description
#' Object used as first layer for the independent
#' ("non connected") signature
#' @export
PrimaryModalityObj <- R6::R6Class("PrimaryModalityObj", #nolint
                                  inherit = ModalityObj,
                                  public = list(
                                    #' @description
                                    #' Initialize `PrimaryModalityObj` object.
                                    #'
                                    #' @param nb_ent : Integer representing the number of entities
                                    #' (genes or transcripts)
                                    #' @param p_up : Proportion of up regulated entities
                                    #' @param p_down : Proportion of down regulated entities
                                    #' @param entity_id (optional) : Character vector with the entity names.
                                    #' Default `NULL`.
                                    initialize = function(nb_ent, p_up, p_down, entity_id = NULL) {
                                      private$nb_ent <- nb_ent
                                      private$p_up <- p_up
                                      private$p_down <- p_down
                                      probas <- list(p_up = p_up, p_nr = 1 - (p_up + p_down), p_down = p_down)
                                      
                                      private$mod_vect <- sample(seq(1, 3),  size = nb_ent,
                                                                 replace = TRUE, prob = probas)
                                      names(private$mod_vect) <- entity_id
                                    }
                                  )
)
#' @title
#' Secondary modality object
#' @description
#'  used as first layer for the connected signature
#' @export
SecondaryModalityObj <- R6::R6Class("SecondaryModalityObj", #nolint
                                    inherit = ModalityObj,
                                    public = list(
                                      #' @description
                                      #' Initialize `SecondaryModalityObj` object.
                                      #'
                                      #' @param prim_mod_obj : `PrimaryModalityObj` used as base for the
                                      #' simulation of the connected object.
                                      #' @param mod_transition : Modality transition "default"
                                      #' or "customized".
                                      #' @param connectivity_score (optional) : Double between -1 and 1
                                      #' representing the "connectivity" between the two perturbations.
                                      #' Used to generate the transition matrix if the modality transition is not
                                      #' "customized".
                                      #' @param nr_noise (optional) : Double between 0 and 0.5 representing the
                                      #' noise induced by the non-deregulated entities during the "default"
                                      #' transition. Default `NA`.
                                      #' @param transition_mat (optional) : Transition matrix used if
                                      #' `mod_transition` = "customized". Default `NA`.
                                      initialize = function(prim_mod_obj,
                                                            connectivity_score = NA,
                                                            mod_transition = c("default", "customized"),
                                                            transition_mat = NA, nr_noise = NA) {
                                        private$nb_ent <- prim_mod_obj$get_nb_ent()
                                        private$connectivity_score <- connectivity_score
                                        private$mod_transition <- mod_transition
                                        private$nr_noise <- nr_noise
                                        
                                        prim_mod_vect <- prim_mod_obj$get_values()
                                        private$transition_mat <- switch(mod_transition,
                                                                         # change the name of the transition
                                                                         default = private$sym_transition_mat_generator(
                                                                           connectivity_score = connectivity_score, nr_noise = nr_noise),
                                                                         customized = transition_mat
                                        )
                                        private$mod_vect <- vapply(prim_mod_vect, function(x) {
                                          sample(seq(1, 3), size = 1,
                                                 replace = TRUE, prob = private$transition_mat[x, ])
                                        }, numeric(1))
                                        
                                        names(private$mod_vect) <- names(prim_mod_vect)
                                      },
                                      #' @description
                                      #' Get the transition matrix
                                      #'
                                      #' @return transition matrix
                                      get_transition_mat = function() {
                                        return(private$transition_mat)
                                      }
                                    ),
                                    private = list(
                                      connectivity_score = NULL,
                                      mod_transition = NULL,
                                      nr_noise = NULL,
                                      transition_mat = NULL,
                                      # @description
                                      # Generates a matrix with the symmetric transition probabilities
                                      # between modalities.
                                      #
                                      # @param nr_noise: Double between 0 and 0.5 representing the noise induced
                                      # by the non-deregulated entities.
                                      # @param connectivity_score : Double between -1 and 1 representing the
                                      # "connectivity" between the two perturbations.
                                      # @return A matrix representing the symmetric probability
                                      # transition matrix.
                                      sym_transition_mat_generator = function(connectivity_score, nr_noise) {
                                        p_id <- 0.5 * (1 + connectivity_score) * (1 - nr_noise / 2)
                                        p_op <- 0.5 * (1 - connectivity_score) * (1 - nr_noise / 2)
                                        p_act_ext <- nr_noise / 2
                                        p_nn <- 1 - nr_noise
                                        mat <- matrix(c(p_id, p_act_ext, p_op, p_act_ext,p_nn, p_act_ext, p_op,
                                                        p_act_ext, p_id), ncol = 3, byrow = TRUE)
                                        return(mat)
                                      }
                                    )
)
