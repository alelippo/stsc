rm(list=ls())

require(KEGGREST)
require(KEGGgraph)
require(igraph)

data.path <- getwd()
tryCatch({
  data.path <- dirname(rstudioapi::getSourceEditorContext()$path)  
}, error = function(e){
  args <- commandArgs()
  if(length(args) > 1){
    data.path <- args[2]
  }
})

log.path <- file.path(data.path, 'log.txt')
pathways.path <- file.path(data.path, 'pathways')
adjacency.path <- file.path(data.path, 'adjacency_matrices')
adjacency.matrices.corrected.path <- file.path(data.path, 'adjacency_matrices_corrected')

load.kegg.map <- function(map) {
  print('Loading map...')
  kgml_file <- keggGet(map, 'kgml')
  if (! file.exists(file.path(pathways.path, paste0(map,'_pathway.xml')))) {
    dir.create(pathways.path, recursive = TRUE)
    write(kgml_file, file = file.path(pathways.path, paste0(map,'_pathway.xml')))
  }
  
  kegg.pathway <- parseKGML(kgml_file)
  graph <- KEGGpathway2Graph(kegg.pathway, genesOnly = FALSE, expandGenes=FALSE)
  return(graph)
}

from.node.to.ids <- function(graphnel) {
  print('Converting nodes to ids...')
  graphnel.nodes <- nodes(graphnel)
  node.ids <- c()
  for (graphnel.node in graphnel.nodes){
    names <- paste(graphnel@nodeData@defaults$KEGGNode$nodes[[graphnel.node]]@name, collapse = ';')
    node.ids <- c(node.ids, names)
  }
  names(node.ids) <- graphnel.nodes
  
  print(paste0('Found nodes ids for: ', length(node.ids)))
  return(node.ids)
}

get.cleaned.igraph <- function(graphnel, node.ids) {
  print('Cleaning graph...')
  igraph <- igraph::graph_from_graphnel(graphnel)
  undefined <- names(which(node.ids == 'undefined'))
  to.be.deleted <- c(undefined)
  for (node in to.be.deleted){
    igraph <- igraph::delete_vertices(igraph, node)
  }
  return(igraph)
}

update.cache <- function(node.ids, pathway) {
  print('Updating cache...')
  kegg.cache.path <- file.path(data.path, "kegg_cache.rds")
  if (file.exists(kegg.cache.path)) {
    kegg.cache <- readRDS(kegg.cache.path)
  } else {
    kegg.cache <- list()
  }
  
  node.id <- 1
  while (node.id <= length(node.ids)) {
    nodes <- unlist(strsplit(node.ids[node.id], ';'))
    for (node in nodes) {
      if (node != 'undefined' && ! node %in% names(kegg.cache)) {
        print(paste0(node.id, '/', length(node.ids)))
        print(paste0('Download information for: ', node))
        kegg.data <- tryCatch({
          keggGet(node)[[1]]
        }, error = function(e) {
          message("Request error: ", e$message)
          Sys.sleep(30)
          NULL
        })
        if (!is.null(kegg.data)) {
          type <- names(kegg.data$ENTRY)
          if (type == 'CDS' || type == 'miRNA' || type == 'ncRNA' || type == 'KO' || type == 'Tight') {
            names <- ifelse(is.null(kegg.data$SYMBOL), node, strsplit(kegg.data$SYMBOL[1], ", ")[[1]])
            kegg.cache[[node]] <- names[1]
          } else if (type == 'Compound' || type == 'Glycan' || type == 'Drug') {
            names <- ifelse(is.null(kegg.data$NAME), node, gsub(";$", "", kegg.data$NAME[1]))
            kegg.cache[[node]] <- names
          } else if (type == 'Pathway') {
            names <- kegg.data$NAME[1]
            kegg.cache[[node]] <- names
          } else {
            log.warning <- paste('>>> Unknown type', type, 'for node', node, '(pathway:', pathway, ') <<<')
            warning(log.warning)
            cat(log.warning, file = log.path, append = TRUE, sep = "\n")
          }
          saveRDS(kegg.cache, file = kegg.cache.path)
        } else {
          node.id <- node.id-1
          break
        }
      }
      node.id <- node.id+1
    }
  }
  
  return(kegg.cache)
}

get.duplicated.nodes <- function(node.ids) {
  print('Retreiving duplicated nodes...')
  condensed.names <- sapply(unique(node.ids), function(node) {
    paste(names(node.ids)[node.ids == node], collapse = ";")
  })
  duplicated.condensed.names <- condensed.names[grep(";", condensed.names)]
  return(duplicated.condensed.names)
}

condense.graph <- function(igraph, duplicated.condensed.names) {
  print('Condensing graph...')
  verteces.to.delete <- c()
  for (duplicated.condensed.name in duplicated.condensed.names) {
    decollapses <- unlist(strsplit(duplicated.condensed.name, ";"))
    decollapses.nodes <- decollapses[decollapses %in% V(igraph)$name]
    for (decollapse.id in seq(decollapses.nodes)) {
      decollapse <- decollapses.nodes[decollapse.id]
      if (decollapse.id == 1) {
        first.collapse <- decollapse
      } else {
        for (neighbor in names(neighbors(igraph, decollapse, mode = 'in'))) {
          if (get.edge.ids(igraph, c(neighbor, first.collapse)) == 0) {
            igraph <- igraph::add_edges(igraph, c(neighbor, first.collapse))
          }
        }
        for (neighbor in names(neighbors(igraph, decollapse, mode = 'out'))) {
          if (get.edge.ids(igraph, c(first.collapse, neighbor)) == 0) {
            igraph <- igraph::add_edges(igraph, c(first.collapse, neighbor))
          }
        }
        igraph <- igraph::delete_vertices(igraph, decollapse)
      }
    }
  }
  E(igraph)$weight <- 1
  
  return(igraph)
}

from.ids.to.alias <- function(node.ids, kegg.cache) {
  print('Converting ids to alias...')
  node.alias <- sapply(names(node.ids), function(alias.ids) {
    alias <- node.ids[alias.ids]
    ids <- unlist(strsplit(alias, ';'))
    replaced.values <- sapply(ids, function(id) {
      if (id %in% names(kegg.cache)) {
        return(kegg.cache[[id]])
      } else {
        return(id)
      }
    })
    paste(replaced.values, collapse = ';')
  })
  names(node.alias) <- names(node.ids)
  return(node.alias)
}

rename.vertices <- function(graph, mapping) {
  print('Renaming vertices...')
  V(graph)$name <- mapping[V(graph)$name]
  return(graph)
}

save.adjacency.matrix <- function(graph, map) {
  print('Saving adjacency matrix...')
  adjacency.matrix <- as_adj(graph)
  
  if (! file.exists(file.path(adjacency.path, paste0('adjacency_matrix_', map, '_R.tsv')))) {
    dir.create(adjacency.path, recursive = TRUE)
    write.table(as.matrix(adjacency.matrix), file = file.path(adjacency.path, paste0('adjacency_matrix_', map, '_R.tsv')), sep = '\t', quote = FALSE, col.names = NA)
  }
}

cat(paste("\n>>> Log started at", Sys.time(), "<<<"), file = log.path, append = TRUE, sep = "\n")

species <- 'hsa'
choosen.map <- c('05022', '05200', '05010', '05014', '04010')
maps <- paste0(species, choosen.map)
for (map in maps) {
  print(paste('Working on map:', map))
  graphnel <- load.kegg.map(map)
  node.ids <- from.node.to.ids(graphnel)
  igraph <- get.cleaned.igraph(graphnel, node.ids)
  kegg.cache <- update.cache(node.ids, map)
  duplicated.condensed.names <- get.duplicated.nodes(node.ids)
  condensed.graph <- condense.graph(igraph, duplicated.condensed.names)
  node.alias <- from.ids.to.alias(node.ids, kegg.cache)
  renamed.graph <- rename.vertices(condensed.graph, node.alias)
  save.adjacency.matrix(renamed.graph, map)
}
