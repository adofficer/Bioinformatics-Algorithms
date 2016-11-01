# Bioinformatics-Algorithms
Final Project Repo

Code and Data from Winter 2016 Bioinformatics Algorithms Final Project

Data used courtesy of the DeCoN Project: http://decon.fas.harvard.edu/

Special thanks to Loyal Goff for the master list of MGI ID's to query the DeCoN API


# Functionality

Uses a genetic algorithm to determine the optimum N RNA transcripts that can be used to reliably distinguish between three different cell type distributions within the mouse forebrain irrespective of developmental age: embryonic 15 weeks, 16 weeks, 18 weeks and pup 1 week.

In brief, several random subsets of N RNA transcripts are chosen and each subset is subjected to PCA followed by K-means clustering in the first two, or three, dimensions depending on user-defined fields.  The twelve data points, three cell types from each of the four different ages of mouse, are then binned into the 3 K-means clusters and these bins are subjected to the objective function which penalizes heterogeneous bins and bins that contain too many members.  The top scoring subsets are then subjected to genetic evolution via random crossovers and point mutations at customizable rates. The algorithm continues until a perfect result is found or until it converges at a local maximum and returns the ID's of the N RNA sets that scored the highest.

Includes code to scrape data from the DeCoN API, munge into a unified data structure and then evaluate the genetic algorithm and visualize results in either 2- or 3-D PCA plots using K-means clustering (no. clusters = 3) to differentiate the cell types.


# Notes and Observations

For N > 4 the energy landscape of the objective function has many global maxima (perfect recall) and numerous local maxima suggesting there are several subsets of RNA transcripts that will differentiate between these cell distributions

For N <= 4 the objective function begins to lose effectiveness due to the nature of performing PCA in low dimensions.  PCA into K-means is quite similar to a naive machine learning classifier, choosing alternative classification methods will overcome this.


# Usage

Raw data from the DeCoN dataset can be found in the all_decon_data.pkl, a pickled pandas DataFrame object.  Code for requerying the API can be found in final_project.py but this is not advised as it involves several tens of thousand queries to the DeCoN API and most queries return empty JSON files.  Contact Loyal Goff, loyalgoff@jhmi.edu, for alternative methods of accessing the raw data.

P-value filtering of the dataset is done to determine significantly regulated genes and to shrink the search space.  The current implementation uses a one-way ANOVA, also known as an f-test, to assign p-values, and then filters out those transcripts that do not pass a user-defined cutoff.  Other methods can also be used, most commonly a moderated-T test as used in the edgeR package, but this code has not been ported over to Python as of the date of writing.  Bonferroni corrections can be applied or one can filter by q-value, see the qvalue.py file for more information and one method of FDR correction for multiple hypothesis testing.

The genetic algorithm itself is run by the genet_alg function which takes as input a pandas DataFrame consisting of any number of rows that correspond to the twelve different samples as outlined in the DeCoN paper.  set_size defines the number of RNA transcripts the algorithm will be searching for and targ_max is the target score of the objective function, a perfect score is 100.  top_n is used to determine the number of subsets carried through each iteration and how many survive the pruning step.  Alternative arguments mut_rate and max_loop are to control the rate of mutation of the genetic algorithm, default 0.1, or 10% mutation rate, and the maximum number of iterations the code is allowed to proceed before completing and returning the results of the last iteration.

# References

Molyneaux, B. J.*, Goff, L. A.*, Brettler, A. C., Chen, H.-H., Brown, J. R., Hrvatin, S., et al. (2014). DeCoN: Genome-wide Analysis of In Vivo Transcriptional Dynamics during Pyramidal Neuron Fate Selection in Neocortex. Neuron (in press). doi:10.1016/j.neuron.2014.12.024. *Authors contributed equally
