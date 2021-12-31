Inocula analysis file:

invivo_pairwise_inocula.txt

To be analysed with R script 

invivo_pairwise_inocula.R

File strain-strain interactions measured after 10 days:

invivo_pairwise.txt


To be analysed with R script 

invivo_pairwise.R 

-We analyse CFUs counts and plot them. 

-Then we check for contaminations. We make stacked bar plots to look qualitatively for contaminations and also we calculate for each sample the proportion of the non-focal strains (contaminations). We look at how this distributes. Decide for threshold to determine contamination = 0.1. 
List of samples that we lose:

20	4.
20	5.
78	3.
78	1.
78	2.
78	4.
8	3.
78	5.
72	5.
70	5.
69	1.
15	1.
70	3.
69	5.
70	2.
68	1.
70	1.
37	2.
30	3.

Monocolonisations 69, 70 and 78 (ESL0353, ESL0183 and ESL0351) were repeated and substituted with clean data.
And the other contaminated replicates in the list are deleted.

-Analyse number of counts and delete samples with < 1000 counts:

P72_5	60.
P72_1	122.
P66_5	155.
P76_5	165.
P64_3	265.
P62_1	300.

For some samples we don't have CFUs counts so I need to exclude them at this point from my dataset:

32_5.
44_5.
45_5.
74_3.
75_5.

New file= 
invivo_pairwise_up.txt

-Normalisation with CFUs

-We start to characterise interactions. The idea is to calculate each strains effect on the other. Once we have this we can use this information to define interactions. 

We calculate the effect of A on B taking the log2 of the normalised counts of A (by replicate) in presence of B divided by the mean of the normalised counts of A in mono-colonisation.

We combine the calculated effects to define interactions using cartesian plots. We only take one symmetry axis of the cartesian plot. 
We convert the cartesian coordinates into polar coordinates (r = distance from center, length of vector; theta = angle describing position of vector on plane) and we transform theta in order to have a range of angles between -90˚ and +90˚. 
> So the radius gives us an idea of the strength of the interactions between the two strains.

> And the angle theta gives us an idea of the type of interaction between the two strains (e.g. if included between -90˚and -45˚ the interaction is mutual inhibition.

-polar stats to indicate which interactions are significant.

-Once that the interactions are defined we want to correlate the interactions strength (radius) with the phylogenetic distances between the different strains.
Basically we calculated using ape (R package) the phylogenetic distance between the 12 strains based on an alignment including 40 core-genes. We correlated this with the radius. No correlation. We also correlated the radius with the Jaccard distances calculated based on the accessory genome (genes shared by 11 or less strains). No correlation. 

-At this point we started to analyse the strains-proportions within the co-colonisations.
I tried by calculating a number (calc50) that would say how much it is missing to get to an even distribution of the two strains (50-50). It looks like there is more distance from a 50-50 situation in pairs including strains from the same species. 



