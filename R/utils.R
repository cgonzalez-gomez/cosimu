#' Infers cosimu distribution parameters from the log2 fold-change expression of
#' a real data vector.
#' 
#' @param lfc_vector vector of log2 fold-change values.
#' @param min.size the minimum size of the modes. Default tos to 0.01, eliminating modes
#' with less than 1% of the distributional area. See `LaplacesDemon` for more details. 
#' @param sym_mode boolean set to `TRUE` to return a symmetrical parametrization
#' for up and down regulated modalities. Default to `FALSE`.
#' @param plot boolean set to `TRUE` to plot the log2 fold change distribution and
#' its modes, in order to get symmetrical distributions for the deregulated
#' modalities. Default to `FALSE`.
#' @param thresh_up threshold to filter the up-regulated modes. Default to 0.1.
#' @param thresh_down threshold to filter the down-regulated modes. Default to -0.1.
#' @param scale scale of the gamma distributions. Needs to be fixed in order to
#' infer the associated shapes. Default to 0.2.
#' @param nr_mean mean of the Gaussian distribution of the NR modality.
#' Default to 0.
#' @param nr_sd standard deviation of the Gaussian distribution of the NR
#' modality. Default to 0.
#' @param save_yaml character string naming a file to save the inferred parameters.
#' Default to `NULL`.
#' @param return boolean set to `TRUE` to return the parameter values. Note that the
#' quantile functions are only determined during the loading process, see `load_dist`. 
#' Default to `FALSE`.
#' 
#' @import LaplacesDemon yaml purrr
#' @export
infer_dist <- function(lfc_vector,min.size = 0.01, sym_mode = FALSE, plot = FALSE,
                        thresh_up = 0.1, thresh_down = -0.1, scale = 0.2, nr_mean = 0,
                        nr_sd = 0.15, save_yaml=NULL, return = FALSE){
  if(is.unimodal(lfc_vector, min.size = min.size)){
    stop("Your distribution is unimodal. Please try with a different vector or
         modify the min.size parameter")
  }
  mean_gamma <- function(shape,scale) shape*scale
  
  modes_res = Modes(lfc_vector,min.size = min.size)
  
  sel.up <- which(modes_res$modes > thresh_up)
  modes.up <- map(modes_res, ~.x[sel.up])
  names(modes.up) <- names(modes_res)
  
  sel.down <- which(modes_res$modes < thresh_down)
  modes.down <- map(modes_res, ~.x[sel.down])
  names(modes.down) <- names(modes_res)
  
  sel.nr <- which(modes_res$modes > thresh_down & modes_res$modes < thresh_up)
  modes.nr <- map(modes_res, ~.x[sel.nr])
  names(modes.nr) <- names(modes_res)
  
  shapes.up <- map_dbl(seq(1,length(modes.up$modes)), function(i){
    if(i==1){
      return(mean(c(max(modes.nr$modes),modes.up$modes[1]))/scale+1)
    }else{
      return(mean(c(modes.up$modes[i-1],modes.up$modes[i]))/scale+1)
    }
  })
  
  dens.up <- round(map_dbl(seq(1,length(modes.up$mode.dens)), function(i){
    if(i==1){
      nr <- which(modes.nr$modes==max(modes.nr$modes))
      return(mean(c(modes.nr$mode.dens[nr],modes.up$mode.dens[1])))
    }else{
      return(mean(c(modes.up$mode.dens[i-1],modes.up$mode.dens[i])))
    }
  }),digits = 3)
  
  shapes.down <- abs(map_dbl(seq(1,length(modes.down$modes)), function(i){
    if(i==1){
      return(mean(c(min(modes.nr$modes),modes.down$modes[1]))/scale+1)
    }else{
      return(mean(c(modes.down$modes[i-1],modes.down$modes[i]))/scale+1)
    }
  }))
  
  dens.down <- round(map_dbl(seq(1,length(modes.down$mode.dens)), function(i){
    if(i==1){
      nr <- which(modes.nr$modes==min(modes.nr$modes))
      return(mean(c(modes.nr$mode.dens[nr],modes.down$mode.dens[1])))
    }else{
      return(mean(c(modes.down$mode.dens[i-1],modes.down$mode.dens[i])))
    }
  }),digits = 3)
  if(!sym_mode){
    if(plot){
      hist(lfc_vector, breaks=100, freq=F)
      abline(v = modes_res$modes, col="red")
      x <- seq(0,8,0.1)
      for(i in 1:length(dens.up)){
        lines(x,dens.up[i]*dgamma(x,shape=shapes.up[i],scale = scale), col="orange")
      }
      
      for(i in 1:length(dens.down)){
        lines(-x,dens.down[i]*dgamma(x,shape=shapes.down[i],scale = scale), col="blue")
      }
    }
    
    parameters <- list(prop_sm_up = dens.up/sum(dens.up),
                       prop_sm_down = dens.down/sum(dens.down),
                       up_means = purrr::map2_dbl(shapes.up,rep(scale,length(shapes.up)),~mean_gamma(.x,scale=.y)),
                       down_means = -purrr::map2_dbl(shapes.down,rep(scale,length(shapes.down)),~mean_gamma(.x,scale=.y)),
                       nr_means = nr_mean,
                       nr_sd = nr_sd,
                       shapes.up = shapes.up,
                       scales.up = rep(scale,length(shapes.up)),
                       shapes.down = shapes.down,
                       scales.down = rep(scale,length(shapes.down)),
                       alpha= sum(dens.up)/sum(c(dens.up,dens.down)))
  }else{
    if(length(modes.up$modes)!=length(modes.down$modes)){
      stop("A symmetrical distribution can't be inferred because the number of modes differ in up and down regulated modalities. Try a different parametrization")
    }
    modes.deg <- rowMeans(data.frame(up= modes.up$modes,down=abs(modes.down$modes)))
    modes.deg[1] <- mean(c(modes.deg[1],abs(mean(modes.nr$modes)))) 
    shapes.deg <- map_dbl(modes.deg, ~.x/scale+1)
    dens.deg<- rowMeans(data.frame(up= dens.up,down=dens.down))
    if (plot){
      hist(lfc_vector, breaks=100, freq=F)
      abline(v = modes_res$modes, col="red")
      x <- seq(0,8,0.1)
      for(i in 1:length(dens.deg)){
        lines(x,dens.deg[i]*dgamma(x,shape=shapes.deg[i],scale = scale), col="orange")
      }
      
      for(i in 1:length(dens.deg)){
        lines(-x,dens.deg[i]*dgamma(x,shape=shapes.deg[i],scale = scale), col="blue")
      }
    }
    relative_prop <-dens.deg
    scales <- rep(scale, length(shapes.deg))
    parameters <- list(prop_sm_up = relative_prop/sum(relative_prop),
                       prop_sm_down = relative_prop/sum(relative_prop), 
                       up_means = purrr::map2_dbl(shapes.deg,scales,~mean_gamma(.x,scale=.y)),
                       down_means = -purrr::map2_dbl(shapes.deg,scales,~mean_gamma(.x,scale=.y)),
                       nr_means = nr_mean,
                       nr_sd = nr_sd,
                       shapes.up = shapes.deg,
                       scales.up = scales,
                       shapes.down = shapes.deg,
                       scales.down = scales,
                       alpha = 0.5)
  }
  if(!is.null(save_yaml)){
    write_yaml(parameters,file = save_yaml)
  }
  if(return){
    return(parameters)
  }
}

#' Loads the YAML file containing the distribution parameters.
#' 
#' This function loads the YAML file generated by `infer_dist` and prepares the
#' parameters, including quantile functions, for use in the cosimu function.
#' 
#' @param yaml_file character string naming a file to load the inferred parameters.
#' 
#' @return list including:
#' - prop_sm_up: Vector of proportions of the sub-modalities of the up-regulated
#' entities.
#' - prop_sm_down: Vector of proportions of the sub-modalities of the
#' down-regulated entities.
#' - up_means: Vector of means of the up sub-modalities 
#' - nr_means: Mean of the NR (sub-)modality
#' - down_means: Vector of means of the down sub-modalities
#' - qf_vect_up: List of quantile functions associated with each sub-modality of 
#' the up-regulated modality.
#' - qf_vect_nr: List of quantile functions associated with each sub-modality of 
#' the non-deregulated modality.
#' - qf_vect_down : List of quantile functions associated to each sub-modality
#' of the down-deregulated modality.
#' - alpha: Fraction of deregulated genes that are up-regulated. The fraction of
#' deregulated genes that are down-regulated is 1-alpha. (p_up = alpha*p_deg and
#' p_down = (1-alpha)*p_deg with 0 < p_deg < 1).
#' @import yaml purrr
#' @export
load_dist <- function(yaml_file){
  parameters <- read_yaml(yaml_file)
  parameters$qf_vect_up <- purrr::map2(
    parameters$shapes.up,parameters$scales.up,
    function(shape,scale) {
      local(function(x) {
        p <- parent.env(environment())
        p$shape <- shape
        p$scale <- scale
        qgamma(x,shape,scale = scale)
      })
    }
  )
  parameters$qf_vect_nr <- list(
    function(x) qnorm(x,mean=  parameters$nr_means,sd = parameters$nr_sd)
  )
  parameters$qf_vect_down <- purrr::map2(
    parameters$shapes.down,parameters$scales.up,
    function(shape,scale) {
      local(function(x) {
        p <- parent.env(environment())
        p$shape <- shape
        p$scale_r <- scale
        -qgamma(x,shape,scale = scale)
      })
    }
  )
  return(parameters)
}
#' Loads the YAML file containing the parameters needed for a simple simulation.
#' 
#' @param yaml_file character string naming a file to load the inferred parameters.
#' The file must have the structure shown in the cosimu_application repository
#' simple_yaml_structure.yaml
#' 
#' @return list of parameters
#' @import yaml purrr
#' @export
load_simple_param <- function(yaml_file){
  parameters <- read_yaml(yaml_file)
  dist <- parameters$dist
  # dist
  qf_vect_up <- purrr::map2(
    dist$shapes.up,dist$scales.up,
    function(shape,scale) {
      local(function(x) {
        p <- parent.env(environment())
        p$shape <- shape
        p$scale <- scale
        qgamma(x,shape,scale = scale)
      })
    }
  )
  qf_vect_nr <- list(
    function(x) qnorm(x,mean= dist$nr_means,sd = dist$nr_sd)
  )
  qf_vect_down <- purrr::map2(
    dist$shapes.down,dist$scales.up,
    function(shape,scale) {
      local(function(x) {
        p <- parent.env(environment())
        p$shape <- shape
        p$scale_r <- scale
        -qgamma(x,shape,scale = scale)
      })
    }
  )
  list_pDEG <- parameters$signatures$primary$list_pDEG
  nb_tech_rep <- parameters$signatures$nb_tech_rep
  primary_param <- map(rep(list_pDEG,
                           each=nb_tech_rep),
                       function(pDEG){
    list(p_up= pDEG[1],p_down=pDEG[2],prop_sm_up=dist$prop_sm_up,
         prop_sm_down=dist$prop_sm_down,
         qf_vect_up=qf_vect_up,qf_vect_down=qf_vect_down,
         qf_vect_nr=qf_vect_nr,up_means= dist$up_means,
         nr_means=dist$nr_means,down_means= dist$down_means)
  })
  
  mod_transition_params <- expand.grid(parameters$signatures$secondary$connectivity_score,
                                       parameters$signatures$secondary$nr_noise)
  colnames(mod_transition_params) <- c("connectivity_score","nr_noise")
  
  if(!is.null(parameters$signatures$secondary$transition_mat)){
    transition_mat = matrix(parameters$signatures$secondary$transition_mat, ncol=3,nrow = 3)
  }else{
    transition_mat=NULL
  }
  secondary_param = map(seq(1,nrow(mod_transition_params)),function(i){ 
    list(mod_transition=parameters$signatures$secondary$mod_transition,
         submod_transition=parameters$signatures$secondary$submod_transition,
         qf_vect_up=qf_vect_up,
         qf_vect_down=qf_vect_down,
         qf_vect_nr=qf_vect_nr,
         connectivity_score= mod_transition_params$connectivity_score[i],
         transition_mat=transition_mat,
         nr_noise= mod_transition_params$nr_noise[i],
         prop_sm_up=dist$prop_sm_up,
         prop_sm_down=dist$prop_sm_down,
         proba_transition=parameters$signatures$secondary$proba_transition, 
         rho_submod= parameters$signatures$secondary$rho_submod,
         copula_submod = parameters$signatures$secondary$copula_submod,
         eps_submod = parameters$signatures$secondary$eps_submod,
         optim_method_submod = parameters$signatures$secondary$optim_method_submod,
         up_means= dist$up_means,
         nr_means=dist$nr_means,
         down_means= dist$down_means,
         copula_prob= parameters$signatures$secondary$copula_prob,
         theta_prob=parameters$signatures$secondary$theta_prob,
         nbins_prob=parameters$signatures$secondary$nbins_prob)
  }
  )
  idx_primary_base <- rep(seq(1,length(primary_param)), each=length(secondary_param))
  idx_secondary_base <- rep(seq(1,length(secondary_param)),times=length(primary_param))
  
  ## Global parametrization data frame
  gb_param_df <- data.frame(primary_base= idx_primary_base,
                            secondary_base=idx_secondary_base)
  gb_param_df <- cbind(gb_param_df,mod_transition_params, 
                       do.call(
                         rbind,rep(list_pDEG, 
                                   each=nb_tech_rep*length(secondary_param))
                       )
  )
  colnames(gb_param_df)[c(5,6)] <- c("p_up","p_down")
  gb_param_df$secondary_id <- paste0("secondary_",seq(1,nrow(gb_param_df)))
  gb_param_df$primary_id <- paste0("primary_",gb_param_df$primary_base)
  return(list(entity_id = parameters$signatures$entity_id,
              base_expression = parameters$signatures$base_expression,
              initial_base = parameters$signatures$initial_base,
              primary_param= primary_param,
              secondary_param =secondary_param,
              gb_param_df =gb_param_df))
}
