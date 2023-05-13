# PhageComparativeGenomics
### dataset for use
The provided dataset comprises the viral genus Abidjanvirus and Phifelvirus, with 6 and 7 genomes of complete phages respectively  
File --> `examples.fasta`

## CDS prediction

## ANI and AAI

## Clustering analysis

Clustering of bacteriopahges can be done in two main ways:
1. Using whole genome similarity
2. Using predicted proteome similarity

### 1 - Genereting all-against-all blastn

```
blastn -query examples.fasta -subject examples.fasta -outfmt "6 qseqid sseqid pident qcovs evalue" > examples_blast.tab
```
The output of this analysis contains all the identity information for pahges in the fasta file, and usally multiple hits for each pair of pahges compared. 

### 2 - Recovering unique *hits* and preparing the data
From the output of the blastn analysis we will calculate the mean identity, coverage and e-value.
Also, if any hit doesnt reach the thresholds of 80% coverage and 1e-5 e-value, their identity values are zeroed out.  


```
import pandas as pd
import numpy as np

# Getting the data, grouping by query id and subject id followed by calcualtion of mean for each feature
data = pd.read_table("examples_blast.tab", names=["query", "subject", "ident", "cov", "evalue"]).groupby(["query", "subject"]).agg({"ident":"mean", "cov": "mean", "evalue":"mean"}).reset_index()

# Zeroing out any hit bellow the stablished threshold

data["ident"] = np.where((data["cov"] < 80.0) | (data["evalue"] > 1e-5), 0, data["ident"])

# keeping only the necessary information (query id, subject id and identity)

data = data[["query","subject","ident"]]

# Saving the dataset

data.to_csv("examples_blast.mean.tab", sep=" ", index=False, header=False)
```
### 3 - Clustering the data using MCL
This step will use the *examples_blast.mean.tab* generated before as input for the MCL commands.  
1. mcxload
`mcxload --stream-mirror -abc examples_blast.mean.tab -o examples.mci -write-tab examples.mci.tab`

Here, all the necessary entry files are generated from the *examples_blast.mean.tab*  

2. mci

`mcl examples.mci -I 2 -use-tab examples.mci.tab -o examples.cluster -te 8`

Here we generate the final whole genome based clusters. The output is a simple .txt file in which each line is a cluster, and contains all the pahges in the cluster separetd by tab (\t)

### 4 - Retrieving cluster information and plotting

From the mci output *examples.cluster* and input *examples_blast.mean.tab* we will generate a dataframe in python and from this generate a vizualization in python  

```
import pandas as pd
import networkx as nx
from pyvis.network import Network
import random

###################################
### Reading cluster information ###
###################################

tmp_cluster_list = []
count = 1

with open("examples.cluster", "r") as c:
    for row in c:
            splitted_row = row.strip().split("\t")
            tmp_cluster_list.append([f"cluster{count}", splitted_row])
            count+=1
clusters_df = pd.DataFrame.from_records(tmp_cluster_list, columns=["Cluster", "Phage in cluster"])
clusters_df = clusters_df.explode("Phages in cluster").reset_index(drop=True)

#########################
### Adding some color ###
#########################

## --> The colors generated are random <-- ##

clusters = clusters_df["Cluster"].unique().tolist()

number_of_colors = len(clusters)

color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

cluster_colors = dict(zip(clusters, color))
colors_df = pd.DataFrame.from_dict(cluster_colors, orient="index").reset_index()

colors_df.rename(columns={"index":"Cluster", 0:"color"}, inplace=True)

### Creating the final dataframe with colors ###

with_color_df = pd.merge(clusters_df, colors_df, on="Cluster", how="left")

## Saving data
with_color_df.to_csv("example_clusters_colors.csv", index=False)

## Finishing creating the dictionary to color our network ##

nodes_attr = with_color_df.set_index("Phage in cluster").to_dict(orient="index")

########################
### Plotting Network ###
########################

network_data = pd.read_csv("examples_blast.mean.tab", sep="\t", names=["source", "target", "ident"])

G = nx.from_pandas_edgelist(network_data,
                            source="source",
                            target="target",
                            edge_attr="ident"
                            )
G.remove_edges_from(nx.selfloop_edges(G))
nx.set_node_attributes(G,nodes_attr)

net = Network(notebook=True)
net.from_nx(G)

net.save_graph("examples_network.html")
```

With that done, we now have a interactive network in .html

![network_view](examples_network_names.png)

This rendered network was colored by clsuter, so each color represents a different cluster, the nodes represent all phages analyzed and edges representes how they are conected based on similarity.



