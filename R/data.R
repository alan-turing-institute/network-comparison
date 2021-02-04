#' Protein-protein interaction (PPI) networks for 5 microorganisms
#'
#' A dataset containing the protein-protein interaction networks for the
#' following 5 microorganisms
#' \itemize{
#'  \item EBV
#'    \itemize{
#'      \item Common name: Epstein Barr virus
#'      \item Scientific name: Human gammaherpesvirus 4
#'      \item TaxonomyID: 10376
#'    }
#'  \item ECL
#'    \itemize{
#'      \item Common name: E.coli
#'      \item Scientific name: Escherichia coli
#'      \item TaxonomyID: 562
#'    }
#'  \item HSV-1
#'    \itemize{
#'      \item Common name: Herpes simplex virus type 1
#'      \item Scientific name: Human alphaherpesvirus 1
#'      \item TaxonomyID: 10298
#'    }
#'  \item KSHV
#'    \itemize{
#'      \item Common name: Karposi's Sarcoma-Associated Herpesvirus
#'      \item Scientific name: Human gammaherpesvirus 8
#'      \item TaxonomyID: 37296
#'    }
#'  \item VZV
#'    \itemize{
#'      \item Common name: Varicella zonster virus
#'      \item Scientific name: Human alphaherpesvirus 3
#'      \item TaxonomyID: 10335
#'    }
#' }
#'
#' @format A list of \code{igraph} objects.
#' @source \strong{PPI data (EBV, HSV-1, KSHV, VZV):} Fossum E, Friedel CC, Rajagopala SV, Titz B, Baiker A, Schmidt T, et al. (2009) Evolutionarily Conserved Herpesviral Protein Interaction Networks. PLoS Pathog 5(9): e1000570. \url{https://doi.org/10.1371/journal.ppat.1000570}. Data from Table S2 in the supporting information.
#' @source \strong{PPI data (ECL):} Peregrín-Alvarez JM, Xiong X, Su C, Parkinson J (2009) The Modular Organization of Protein Interactions in Escherichia coli. PLoS Comput Biol 5(10): e1000523. \url{https://doi.org/10.1371/journal.pcbi.1000523}
#' @source \strong{Taxonomy ground truth:} NCBI taxonomy database. \url{https://www.ncbi.nlm.nih.gov/taxonomy}
#' @encoding UTF-8
"virusppi"





#' World trade networks from 1985–2014
#'
#' The world trade data set consists of bilateral trade flows between countries for the years 1985–2014 [Feenstra et al., 2005, Division, 2015]. This set is a subset of a dataset that spams the years 1962–2014.
#' 
#'
#' \itemize{
#'  \item wtnets:  List of \code{igraph} objects providing the world trade networks from 1985–2014.
#'  \item Counts:  Pre-computed graphlet counts for the networks in 'wtnets'.
#'  }
#'  
#' @format A list of two elements. The first element, 'wtnets', is a list of \code{igraph} objects providing the world trade networks from 1985–2014. The second element, 'Counts', is a list of pre-computed graphlet counts for the aforementioned networks.
#' @source \strong{World trade networks:} Feenstra RC,Lipsey RE, Deng H, Ma AC, and Mo H. (2005) World trade flows: 1962-2000. Technical report, National Bureau of Economic Research
#' @encoding UTF-8
"worldtradesub"

