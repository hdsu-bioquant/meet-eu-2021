# Meet-EU Session 2021-2022

## Important links

* [Main website](https://hdsu-bioquant.github.io/meet-eu-2021/)
* [Description of the project specifics (for the teachers)](https://docs.google.com/document/d/1UfugMVAMy_YdhosHKuZ0v9oKOideTGDI5eKdwWZIUK4/edit)
* [Discord link](https://discord.gg/kscA9FfUru)
* [Sheet to fill up for group definition and pairing](https://docs.google.com/spreadsheets/d/1Rn3PVAmTvKFZmIiuem7xGvYscvawnPb4Q-Dd7TjSo0E/edit?usp=sharing)
* [Slides of the opening session](http://www.lcqb.upmc.fr/meetu/dataforstudent/slides/opening2022.pdf)
* Videos of 2020 year:

    * Opening -- https://youtu.be/-x5gzjqB6ms  
    * Practical -- https://youtu.be/y5omyj-5X3o  
    * Half-way forum -- https://youtu.be/zvcf7v8Bdas



## Tasks and definitions 

### 1. TAD TEAM goal: create an algoritm for TAD detection with the best precision you can reach.
Implement a solution like Topdom, or Arrowhead (link to the paper on github, have a look at it!!)
Optimise it to only find high-quality TADs.

This challenge deals with optimisation and quality assessment. The goal is to find some idea(s) and develop a critical mind about the parameters in an algorithm. You will be evaluated in how you consider that one TAD has been correctly detected.

Some good ways to solve it :  

- Find only domains containing CTCF at their borders.  
- Implement 2 or 3 algorithms and compare their ouputs with the best parameters.  
- Design an evaluation strategy to validate the outputs of the algorithms.


### 2. CPT TEAM goal: Detect chromatin compartments the best way you can.

Defining two compartments is really easy and we give you most of the bricks to do this 'take in hand' task :  

- Load HiC matrix and epigenomic feature  
- Apply filtering, SCN, Observed under expected and then pearson correlation to HiC (see the slides!!)  
- Extract the first eigenvector on the pearson correlation matrix  
- Order it by the epigenomic feature at the same resolution  
    
DONE (all functions to do it are given.. almost)
    
    
Then, you have to improve it with one of this possible options (sorted by difficulty, but you have to choose one!):  

- Compare existing solutions to detect compartments, using inter-chromosomal contacts on intra-chromosomal contacts -- what is the best?  
- Increase the number of compartments (train hidden markov model for instance) and define for you what the best number of compartments is  
- Increase the resolution -- is the information still correct? 
    

The final product should give as output a set of TADs or a set of compartments. We strongly encourage you to visualise your results in 2D or 3D (see the provided examples). 

## Input and validation data 

We provide you with three types of input data, all available [here](http://www.lcqb.upmc.fr/meetu/)

- Hi-C matrices, taken from [this repository](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl). These are the basic input data that your algorithm should treat. They give observed contacts between chromatin regions. Each matrix corresponds to one or several chromosome(s) from a given organism, studied in a particular experimental condition.

- Annotations from the ChromHMM tool for GM12878 (most documented human cell type), which can be found [here](https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E116_15_coreMarks_stateno.bed.gz). These annotations provide information about epigenetic marks. These marks are indicators of the state of the chromatin and thus can help you determine the optimal number of **compartments** and detect them.

To validate your program, we can rely on gold-standard definitions of the TADs and the compartments. These can be found in the [data repository](http://www.lcqb.upmc.fr/meetu/).

## Hictoolbox 
Contains some useful functions that you might need in your code: how to load a matrix, filter a matrix in the easiest way, how to change the resolution, and how to normalise a matrix with SCN (Cournac & al). The functions are all written for **scipy** sparse matrix.

If you find some bug or have some suggestions of improvement, feel free to create an issue :)

This toolbox does not contain more than: https://gitlab.com/LeopoldC/shrec3d.

## Papers you can't avoid to read

### Introductory texts to understand the context and the basic notions  
We wrote two instructive (or made for you) articles, in French and in English. These are improved versions of the introductory lecture!  
- French versions: https://bioinfo-fr.net/author/mathurin. You must read the last article on TADs and compartments.  
- English version: https://docs.google.com/document/d/1Y8JeQGIkyTmOuSCtA2XfDUNiC4fDy8pKe67o8IEIKRA/edit?usp=sharing
and
 https://docs.google.com/document/d/1ZD12v3K7CvZP_bt-w7ZwHOXRJNMg-1S59gI1VfFoD_E/edit?usp=sharing

We also provide below a list of papers relevant to the project. They are sorted by task (or team type). Sometimes, we indicate the supplementary material of a paper rather than the paper itself.

- A really usefull link: https://gitlab.com/LeopoldC/shrec3d/-/blob/master/main.py. This is a working code example which loads a matrix and generates a 3D representation. Just look the begining of the **main** (50 first lines).

- This paper, from Rao et al., gives basic concepts for all teams: https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4. You should particularly look at the supplementary material. It contains methods for TAD detection AND compartment detection.


### MUST-READ for TAD teams

Basic review of state-of-the-art methods : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5389712/.

Arrowhead from Rao paper, a simple way to understand it : https://www.youtube.com/watch?v=yJSZ-IYbw_s

Lots of good papers synthetized here : https://github.com/3DGenomes/TADbit

For precision increase, pertinent strategies can be found in this paper : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5846465/

This paper contains interesting ideas : https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6088-0

A new interesting idea : https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03652-w

Sounds interesting but not essential : https://www.nature.com/articles/s41586-020-2151-x

This one might be of interest if you want to play with CTCF : https://twitter.com/colinlog/status/1295112780861800448

### MUST-READ for CPT teams

Like all teams, again, read the Rao paper : https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4

Some good ideas in plants, will work with human : https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1694-3

All you need to know about epigenetic and compartment (in drosophila, but they have pretty much the same marks as humans) : https://www.ncbi.nlm.nih.gov/pubmed/28826674


