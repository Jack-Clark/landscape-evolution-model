#Landscape Evolution Modelling

This repository contains most of my work from my research internship at The School of Engineering and Computer Science, Durham University. 
It is only a fraction of the project's codebase, however due to potential licensing issues the full code must remain in a private repository,
and therefore I am only able to show my own code.

The code in this repository contains fast GPU algorithms for computing the flow accumulation from a single flow direction grid and a grid of runoff weights.
Before starting my internship, the fastest algorithm was the [SFD_NoPart_List algorithm](https://github.com/Jack-Clark/landscape-evolution-model/blob/master/parallel-SFD-List.cu), 
which is an optimised version of the correct flow algorithm, described in [this paper](http://community.dur.ac.uk/stephen.mcgough/CV/Papers/2012/Land_paper.pdf). During the internship, 
I managed to redesign the algorithm, making it 2x faster. The faster version is called the [Multiple Retries algorithm](https://github.com/Jack-Clark/landscape-evolution-model/blob/master/process_SFD_multiple_retries.cu).
