`stsc` 
==========

Matlab and R code used to implement and test source-to-sink communicability in the paper

		"Source-to-sink communicability: a new centrality measure for directed acyclic networks" 

by D.Bertaccini, L.Chiricosta, A.Filippo.

The repository consists in two folders: DATASET, CODE.

The DATASET folder contains:

adjacency matrices of five pathways described in the work

The CODE folder contains: 

stsc.m        - MATLAB function that computes the source-to-sink centrality of a directed acyclic graph with adjacency matrix A. 

pathway_preprocessing.R 	- R script that manages the download of the pathway and their reorganization into directed graphs.  

Last modified: 05/06/2025 by Alessandro Filippo 

>> Reference:
>>  D.Bertaccini(bertaccini@mat.uniroma2.it), L.Chiricosta(luigi.chiricosta@irccsme.it), A.Filippo (filippo@mat.uniroma2.it). "Source-to-sink communicability: a new centrality measure for directed acyclic networks", submitted.
>> 
