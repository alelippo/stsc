`stsc` 
==========

Matlab and R code used to implement and test source-to-sink communicability in the paper

		"Source-to-sink communicability: a new centrality measure for directed acyclic networks" 

by D.Bertaccini, L.Chiricosta, A.Filippo.

The repository consists in two folders: DATASET and CODE.

The DATASET folder contains:

Five .tsv files encoding the adjacency matrices of the five pathways analyzed in the paper.

The CODE folder contains: 

pathway_preprocessing.R 	
- R script that manages the download of the pathways from the Kyoto Encyclopedia of Genes and Genomes (KEGG) and their reorganization into directed graphs. The adjacency matrices of the graphs are then stored into .tsv files.

stsc.m
- MATLAB function that computes the source-to-sink communicability of the nodes of a directed acyclic graph with adjacency matrix A. 

pathway_analysis.m
- MATLAB script used to analyse the pathways in the DATASET folder. 

Last modified: 01/07/2025 by Alessandro Filippo 

>> Reference:
>>  D.Bertaccini(bertaccini@mat.uniroma2.it), L.Chiricosta(luigi.chiricosta@irccsme.it), A.Filippo (filippo@mat.uniroma2.it). "Source-to-sink communicability: a new centrality measure for directed acyclic networks", submitted.
>> 
