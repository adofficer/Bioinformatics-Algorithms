# Bioinformatics-Algorithms
Final Project Repo

Code and Data from Winter 2016 Bioinformatics Algorithms Final Project

Data used courtesy of the DeCoN Project: http://decon.fas.harvard.edu/

Special thanks to Loyal Goff for the master list of MGI ID's to query the DeCoN API


# Functionality

Uses a genetic algorithm to determine the optimum N RNA transcripts that can be used to reliably distinguish between three different cell type distributions within the mouse forebrain irrespective of developmental age.

In brief, several random subsets of N RNA transcripts are chosen and each subset is subjected to PCA followed by K-means clustering in the first two, or three, dimensions depending on user-defined fields.  Data points from each of the four different ages of mouse are then binned into each of the 3 K-means clusters and these bins are subjected to the objective function which penalizes heterogeneous bins and bins that contain too many members.  The top scoring bins are then subjected to genetic evolution via random crossovers and point mutations at customizable rates. The algorithm continues until a perfect result is found or until it converges at a local maximum and returns the ID's of the N RNA sets that scored the highest.

Includes code to scrape data from the DeCoN API, munge into a unified data structure and then evaluate the genetic algorithm and visualize results in either 2- or 3-D PCA plots using K-means clustering (no. clusters = 3) to differentiate the cell types.


# Notes and Observations

For N > 4 the energy landscape of the objective function has many global maxima (perfect recall) and numerous local maxima suggesting there are several subsets of RNA transcripts that will differentiate between these cell distributions

For N <= 4 the objective function begins to lose effectiveness due to the nature of the dimensionality reduction of PCA and the algorithm does not efficiently sample the search space.  PCA into K-means is quite similar to a naive machine learning classifier, choosing alternative classification methods could overcome this shortcoming.
