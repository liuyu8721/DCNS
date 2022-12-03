# Early Detection of SARS-CoV-2 Variants through Dynamic Co-mutation Network Surveillance
A co-mutation network surveillance framework was devised to dynamically investigate the nucleotide covarying pattern of weekly sequences, which is expected to provide earlier detection of SARS-CoV-2 variants overcoming genomic data flood, aiding in the response to the COVID-19 pandemic.

We repeatedly executed weekly detection protocols for real-time tracking of circulating co-mutation network using our method. Two main steps are included in this process.

### Step 1 Weekly co-mutation community network
#### Step 1.1 The affinity model for identification of paired co-mutations  
We model a mutation's tendency to be present or absent in a genome where another mutation is already present. This tendency to co-occur is measured by a new index insensitive to prevalence, α, proposed by Mainali et al. (2022) and termed to be an affinity metric of co-occurrence (see the original paper). We develop the 'CooccurMatrix.R' to generate the co-occurrence matrix  between every mutation pair for each weekly data. Then the affinity metrics are calculated through the 'AffinityModel.R', where the algorithm, developed by Mainali et al. (2022) and embeded in the R package 'CooccurrenceAffinity', is invoked.

#### Step 1.2 Co-mutation network and co-mutation communities  
We pick up the co-mutation pairs with false discovery rate <0.001. Each pair of co-mutations will result in a connection or an edge leading to an adjacency matrix which defines the co-mutation network. Because the affinity model indiscriminately identifies homogeneous and heterogeneous co-mutation pairs, respectively abbreviated as HoCPs and HeCPs, we inherited the rate of the co-mutation (RCM) from Qin et al. (2021) to detect HoCPs. The HoCPs identified form an aggregated community structure with groups of strongly linked nodes. Then the Girvan-Newman partition algorithm (Newman et al. 2004) is used to discover these HoCP groups, named with co-mutation communities. All this work is executed by the 'CoNet.R'. 

#### Step 1.3 Weekly co-mutation community tree
The co-mutation communities exhibit hierarchical organization in weekly co-mutation network. This hierarchy can be captured by division of the viral genomes and their hierarchical containment according to the detected communities' presence or not. We built an arborescence, a directed rooted tree, to depict their concatenated containment between these divisions and then used its topological ordering to find the hierarchical relationship. These are done by 'ComutationCommunityTree.R'.

### Step 2 Dynamic creation of a co-mutation community dictionary tree
Step 2.1 Initial dictionary tree  

Step 2.2 Creation of weekly dictionary tree  

Step 2.2.1 Merging current week’s co-mutation communities into dictionary  

Step 2.2.2 Re-creation of current week’s co-mutation community tree  

Step 2.2.3 Union of last week’s dictionary tree and current week’s community tree    


### References
1. Mainali KP, Slud E, Singer MC, Fagan WF. A better index for analysis of co-occurrence and similarity. Sci Adv (2022) 8(4):eabj9204. doi: 10.1126/sciadv.abj9204.
2. Qin L, Ding X, Li Y, Chen Q, Meng J, Jiang T. Co-mutation modules capture the evolution and transmission patterns of SARS-CoV-2. Brief Bioinform (2021) 22(6):bbab222. doi: 10.1093/bib/bbab222.
3. Newman ME, Girvan M. Finding and evaluating community structure in networks. Phys Rev E Stat Nonlin Soft Matter Phys (2004) 69(2 Pt 2):026113. doi: 10.1103/PhysRevE.69.026113.

