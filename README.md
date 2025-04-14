# Introduction

Influenza A viruses are notorious and exemplary pathogens that cross host-species boundaries. By no means exclusively, yet the potential for IAV cross-species emergence is one of the most keenly felt biological threats in the world. The proclivity of these viruses to circulate in multiple vertebrate species, including many of the species reared agriculturally, poses great importance on understanding how and when host-jumps may occur. In addition, huge amounts of resources are spent annually to try and predict whether IAV will jump between species, and what those consequences may entail. To make these predictions, many different angles of influenza virus biology are analysed from evolutionary history to ecological contact networks to structural/molecular modelling of host-pathogen protein interactions to name but a few. However, exploring these features independently can obfuscate patterns and/or lead to pseudo-replication - is the relationship between protein structures already accounted for by viral phylogenetic histories?

Our cross-specialisation approach first uses protein compositional features to train machine learning models on classifying the original host the virus was obtained from. _The features are x and show y, included them because z_

Finally, we compared the results of both modelling and phylogenetic approaches to detect any overlap or, more interestingly, any inferences that were exclusive to only one method.

# Methodology

## Sequence Collection

Sequences were collected from NCBI Virus Database. This dataset was filtered to only include samples with: a) complete nucleotide sequences, missing no more than 10% of the average length for that segment b) metadata containing the virus subtype and original host, at least up to the Class level (Aves or Mammalia) through preferably detailing all the way down to the Family or Species level Nucleotide sequences were then assigned a random identifier (between AAA-000-AAA and ZZZ-999-ZZZ) to associate with metadata throughout analyses. 14,682 avian and 21,737 mammalian sequences were collated at this stage.

Samples from mammalian hosts were then retained only if the original host was within the Canidae, Suidae, Equidae, Hominidae or Phyllostomidae families, as viruses within these hosts were expected to display clear examples of adaptation. Viral sequences of mammalian-origin were then reduced further by excluding samples where the subtype does not circulate naturally within closed systems of that host (i.e. removing spillovers or stuttering chains of transmission, neither of which represent the kind of host adaptation we sought to detect). Included subtypes for each family: Canidae (H3N2, H3N8), Equidae (H7N7, H3N8) and finally Hominidae (H1N1, H1N2, H2N2 and H3N2).

Corresponding protein sequences were then obtained for each sample, most of which were available from the NCBI Virus Database. The few genomic sequences without matching proteins were translated locally using the 'ape' package in R [@paradis_2019aa].

In the future, we would've clustered per host and also possibly split the birds up more instead of assuming that **all** birds are **always** the source of infection.

## Clustering

Proteins were then aligned using MAFFT [@rozewicki_2019aa] into 20 resulting protein fasta files (two for each of the 10 major IAV proteins, one containing mammalian-origin viruses and the other avian-origin). These were then clustered using the linear algorithm 'easy-linclust' within the tool MMSeq2 [@steinegger_2018aa]. Three variables of the clustering algorithm were experimented with in optimising data downsampling: the percentage identity threshold required to delineate sequences from one another (`--min-seq-id 0.5, 0.75, 0.85, 0.9, 0.95, 0.99`), coverage of the sequences within each cluster (`-c 0.5, 0.7, 0.8`) and the way in which the coverage is calculated (`--cov-mode 0, 1`). Ultimately options were set to `--min-seq-id 0.95`, `-c 0.8` and `--cov-mode 1` in order to grant a representative number of samples without introducing redundancy. Ultimately, 6,697 proteins were nominated leading to sequences from 4551 individual viruses.

## Protein Features

Proteins from avian-origin and mammalian-origin viruses retained by MMSeq2 were then aligned into a single fasta file for each protein. Accessory proteins (PB1-F2 and PA-X) were then removed from the analyses as both lie within the coding regions of their respective genomic segments and thus the calculation of any features from sequence data would already include these regions. Similarly, those spliced from alternate reading frames (M2 and NEP) were excluded from any further analyses as they were too short and/or conserved to show patterns of adaptation

OmegaFeaturesCLI [@chen_2022aa] was run over the alignments to produce matrices of protein compositional and physio-chemical properties with the 'get_descriptor' command. Features, and the associated command, are as follows: 2mers ('DPC type 1'), Chou's pseudo-amino acid composition ('PseAAC'), conjoined triad ('CTriad'), and finally Composition (CTDc), & Transition (CTDt) & Distribution (CTDd) values.

Six datasets of protein features were derived from the iFeatureOmega processing: amino acid two-mers (DPC, 400 features), (CTDc, CTDd and CTDt, 39, 195 and 39 features respectively), paac (23 features) and triad (343 features). These were then annotated with the subtype of the origin virus and labelled with the original host.

## Machine Learning Models

Matrices of protein features were then used to create Random Forest classification models, to distinguish between the host from which viruses were sampled.

In order to create test data, subtypes that appeared most frequently (>20 times) were iteratively held out of model construction. In addition to these 23 subtypes, both bat-exclusive subtypes (H17N10 and H18N11) were held as test data. The remaining 70 subtypes comprised $\frac{308}{4551}$ (~7%) sequences. Six Random Forest models were created from each of the feature datasets, for each holdout. These models were then stacked with the use of the 'caretEnsemble' package [@deane-mayer_2024aa]. Models were weighted in order to accommodate for disproportionate sizes of host groups by using the inverse of that group's proportion.

These datasets were then used to train Random Forest models, leading to 2400 models in total: 8 proteins, 6 chemical and compositional features, 25 subtype holdouts and either two-class or multi-class models. Each models' training statistics were extracted for preliminary validation. Models for each protein and subtype were then compiled with the stacking algorithm by caret.

![Model assembly: One model was created for each parameter.](Figures%20&%20Presentables/Model%20Assembly.pdf)

## Phylogenetics

Using sequence data from each of the 8 proteins studied phylogenetic trees were estimated. Initial model parameters such as evolutionary and codon-partitioning models were estimated by ModelFinder within IQTree2 [@kalyaanamoorthy_2017aa,@minh_2020aa], and temporality was confirmed with TempEst [@rambaut_2016aa].

Maximum Likelihood trees were estimated with substitution models suggested by ModelFinder (Table X) and bootstrapped 1000 times for validation.

Phylogenetic models were first trialled on Maximum Likelihood trees with the use of BayesTraits v3.0 [@meade_2022aa]; in the R wrapper "btw" [@griffin_2018aa]. A discrete-trait analysis estimated the origin host of viruses and protein features were modelled as continuous traits. Trait analyses were first modelled using tips features and then ancestral trait reconstruction upon internal nodes, both using multi-state reverse-jump MCMC procedures. Initially, transition rates between discrete states (i.e. Host) were examined across each protein tree before the states of ancestral nodes were estimated.

Specific nodes were selected for ancestral trait reconstruction due to known historical host-switching events. Recognised host switch events include: human viruses that originated from avian and swine reassortments in 1918, 1957 and 1968; emergence of swine viruses into human pandemic strains (H1N1pdm09); avian IAV into horses (H7N7 and H3N8) and the jump of equine H3N8 viruses into canine populations.

## Sequence Reconstruction

RateAncestor: this variable can be set to 0 or 1. If RateAncestor = 1, the program is forced to do two additional analyses, which you can ignore if you do not need the results. First, under a model of variable rates across sites such as the gamma, RateAncestor = 1 forces the program to calculate rates for individual sites along the sequence (output that will be saved in the output rates file), using the empirical Bayes procedure (Yang and Wang 1995). Secondly, RateAncestor = 1 forces the program to perform the empirical Bayesian reconstruction of ancestral sequences (Yang et al. 1995b; Koshi and Goldstein 1996; [Yang 2006 pages 119-124(<http://abacus.gene.ucl.ac.uk/CME/>)]). Ancestral state reconstruction by parsimony (Fitch 1971; Hartigan 1973) is well known (i.e., it is implemented in the PAML program pamp). It can also be achieved using the likelihood/empirical Bayes approach. Often, the two approaches produce similar results, but the likelihood-based reconstruction has two advantages over parsimony: this method uses information from the branch lengths and the relative substitution rates between characters (nucleotides) and provides a measure of uncertainties in the form of posterior probabilities for reconstructed ancestral states.

The results are listed, by site, in the output file rst. You can also use the variable verbose to control how much information you want to be written in this output file. If verbose = 0, only the best nucleotide (the one with the highest posterior probability) at each node at each site is listed, while with verbose = 1 (try 2 if 1 does not work), the full posterior probability distribution from the marginal reconstruction is listed. If the model is homogenous (i.e., if nhomo = 0 or nhomo = 1 have been specified in the control file) and assumes one rate for all sites, both the joint and marginal ancestral reconstructions will be calculated. If the model assumes variable rates among sites like the gamma model, only the marginal reconstructions are calculated. More details about ancestral sequence reconstruction in the section below.

In order to test the models further, protein features of the internal nodes could be calculated. Of course, the majority of these internal nodes do not *have* sequences collected. To circumvent this, a further step of our analyses could be to reconstruct ancestral sequences. By estimating the sequence of an internal node, it can then be re-input into iFeaturesOmega. Then, if our hypothesis stands, the features of these inferred proteins ought to be an average of the two descendants.

CodeML is suggested as a tool for reconstructing ancestral protein sequences.`#0969DA`

# Results

## Dataset Curation

In total 36,419 sequences were found eligible for study (14,682 avian and 21,737 mammalian).

The clustering resulting from removing sequences with a 95% sequence identity left us with 6,697 protein sequences overall, representing 4,551 individual virus genomes. To assist with building the model architecture, any virus genome that had at least one representative protein was analysed fully - hence 4,551 samples for each of the 8 proteins analysed, despite the high sequence homology between many of them. The number of representative proteins called by the clustering algorithm is shown in Table X and full lists of the associated accession codes are available in Supplementary X.

> [!NOTE]
> Note that some viruses had multiple proteins

## Protein Models

Confusion matrices made between known origin host and those predicted by the stack models were used to assess model predictive strength and find evidence of possible host shifts.

Counter-intuitively, a lower model accuracy implied greater evidence of adaptation to the host species; the mislabelling of proteins is indicative of host-specific adaptations to non-native hosts. For example, should a human-origin virus be incorrectly labelled as swine-origin by the model, this supports the hypothesis that the protein retains some signal of adaptation to swine hosts. Overall, confusion matrices were evaluated based on accuracy and F1 statistics (shown in Table Y).

Main outliers (>1%) are pig-human-avian though maybe just caused by the relative proportion of these. The RF models do show good estimation of true-positives, especially for the most frequently appearing samples. Interestingly though, the biggest mistake that the model made was of misclassifying human-origin sequences as pig-origin. In fact, more of the human sequences were classified as swine (8.13%) than correctly identified (7.23%).

1. all models have high spec but sens varies wildly
2. again (maybe this is just repeating) high neg pred and varying pos pred
  - weirdly bat seems to perform better than canine and equine, maybe the sparsity of dog/horse plus their two different subtypes throw things off more than the bats which are always off in their own world
3. Dog is the only time that the models are basically guessing randomly

Overall though, the models perform pretty well given how complex this is and the fact that it is just based on 1039 protein features alone

### Protein by Protein

> [!NOTE]
> How do all these correlate to protein length? Is NS1 doing better purely because it is the shortest and therefore any difference in features is going to be proportionally more dramatic?

No obviously 'badly performing' protein

## Mismatches

When did the model predict things incorrectly? There was little to no pattern to the misclassifications.

## Protein Features

By running ML models and iteratively removing each protein feature, the model performance statistics can be used to infer which features cause performance to drop the greatest. Hence, taking out the features and re-running the model allows us assess the impact that feature had on model performance. The drop in AUC is used to assess the importance of each feature within each ML model.

## Phylogenetics

Alignments of each of the eight proteins were used to estimate trees using IQTree2. 

Could also stick some quick diversity measures in for each alignment, to show what we're working with?

**How correlated are protein features, length and diversity?**

### Protein Trees

**Is it worth adding bootstrap values? Or describing the tree (num of tips/nodes) just to clarify that some sequences were excluded from tree estimation?**

**Is there correlation between charge and branch length/placement? Would need a nice branch-branch measurement e.g. Faith's PD**

### Trait Analyses

Transition rates between hosts, based on entire trees when the ONLY feature input to the phylogenetic model is host. i.e. given the distribution of hosts throughout the tree, what is the likelihood that any given branch will transition from host<sub>1</sub> to host<sub>2</sub>?

Interesting how rarely avian is estimated as the donor, considering reality

Wonder if the pig-to-human values would drop if pdm<sub>2009</sub> were held as a separate group of sequences

*Nice horse-dog rates given how closely related they are*

Also interesting when matrices are asymmetrical, indicates a direction of transfer from the phylogenetic signal

### Ancestral Reconstruction

Describe subtree selection: 
- node number
- number of tips
- why were these subtrees selected
