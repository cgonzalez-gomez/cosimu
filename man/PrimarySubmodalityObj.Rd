% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_SubmodalityObj.R
\name{PrimarySubmodalityObj}
\alias{PrimarySubmodalityObj}
\title{Primary sub-modality object}
\description{
Object used as second layer for the independent ("non connected") signature
}
\section{Super class}{
\code{\link[cosimu:SubmodalityObj]{cosimu::SubmodalityObj}} -> \code{PrimarySubmodalityObj}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{PrimarySubmodalityObj$new()}}
\item \href{#method-clone}{\code{PrimarySubmodalityObj$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="cosimu" data-topic="SubmodalityObj" data-id="get_nb_submod_by_mod">}\href{../../cosimu/html/SubmodalityObj.html#method-get_nb_submod_by_mod}{\code{cosimu::SubmodalityObj$get_nb_submod_by_mod()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="cosimu" data-topic="SubmodalityObj" data-id="get_prop_list">}\href{../../cosimu/html/SubmodalityObj.html#method-get_prop_list}{\code{cosimu::SubmodalityObj$get_prop_list()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="cosimu" data-topic="SubmodalityObj" data-id="get_values">}\href{../../cosimu/html/SubmodalityObj.html#method-get_values}{\code{cosimu::SubmodalityObj$get_values()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize \code{PrimarySubmodalityObj} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PrimarySubmodalityObj$new(mod_obj, prop_sm_up, prop_sm_down, entity_id = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mod_obj}}{: \code{ModalityObjPrimary} corresponding to the first
layer of the independent signature.}

\item{\code{prop_sm_up}}{: Vector of proportions of the
sub-modalities of the up regulated entities.}

\item{\code{prop_sm_down}}{: Vector of proportions of the
sub-modalities of the down regulated entities.}

\item{\code{entity_id}}{(optional) : Character vector with the entity names.
Default \code{NULL}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PrimarySubmodalityObj$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
