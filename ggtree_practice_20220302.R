## Practice R script
## Introduction to ggTree - March 2, 2022
## by Taylor Davedow


# Data used in this workshop was modified from
# Arteaga et al. 2020. https://doi.org/10.1099/mgen.0.000340

# ggtree book resource: https://yulab-smu.top/treedata-book/ 

### Getting set up
# create a folder called "intro_ggtree" and within that create sub directories
# for "data", "output", "script"  and "docs"
# set intro_ggtree as the working directory
# save script into "script"
# save slides into "docs"
# "output" is where we will send the final tree
# download the newick file and the two xlsx files and store them in "data"


##********************************************************
### 1. Install and load packages ----
##********************************************************
# You can use different methods to install packages
#   A) Install button in the Packages Window
#
#   B) Select from menu bar Tools > Install Packages > type "package name" > Install
#  
#   C) Remove the "#" and run the following 3 lines of code.
# 
# install.packages("readxl")
# install.packages("BiocManager")
# install.packages("treeio")
# install.packages("tidyverse")
# install.packages("phytools")

# ggtree package must be installed after installing and loading BiocManager

library(BiocManager)

BiocManager::install("ggtree")

library(readxl) # for reading in xl files
library(ggtree) # for building tree
library(treeio) # for read.newick function
library(phytools) # for midpoint.root (also has read.newick option)
library(tidyverse) 

# we also will be using ggplot2 which should automatically load in
# with ggtree

# install ggtree
install.packages("BiocManager")
BiocManager::install("ggtree")


##********************************************************
## 2. Load in files ----
##********************************************************

# tree file
tree <- read.newick("data/sample_tree.newick")

# metadata file 
metadata <- read_xlsx("data/metadata.xlsx")


# blast results
blast_raw <- read_xlsx("data/blast_results.xlsx")


##********************************************************
## 3. Create a basic tree ----
##********************************************************

# basic tree
ggtree(tree)

# or
 
tree %>% ggtree()

# check the help page for ggtree to find further usages, arguments and
# references to explore!
?ggtree

##********************************************************
## 4. Changing the tree layout ----
##********************************************************
##*
# LAYOUTS

ggtree(tree,
       layout = "rectangular") # did anything change? this is the default

ggtree(tree, 
       layout = "circular") 

ggtree(tree, 
       layout = "slanted") 

# be careful with branch length
# setting branch.length = "none" will create a cladogram

ggtree(tree,
       layout = "rectangular", branch.length = "none")

# other options to consider: rooting 
# try changing the position of the root node using root.position argument

# adding a midpoint.root (part of phytools package)
ggtree(midpoint.root(tree))


# adding/positioning the scale

ggtree(tree) + 
       geom_treescale(x = 0.2, y = 0,  # change the position of the scale
                      width = 0.01) # change the width of the scale from default (0.03)


# tree manipulation

# many different functions for tree manipulation
# example: scale, collapse, expand, flip or group clades
# we can also zoom in to a particular clade (viewClade function)

# first we have to check the node numbers so we can refer to them
# when using these functions

ggtree(tree) + 
  geom_text2(aes(subset=!isTip, 
                 label=node), 
             hjust = -.3)


# zoom in to a clade and show only clade of interest
ggtree(tree) %>% 
  viewClade(node = 37)

# or show both the original and zoomed in tree
ggtree(tree) + 
  viewClade(node = 37)

# scale
ggtree(tree) %>% 
  scaleClade(node = 41, 
             scale = 0.1)

# collapse
ggtree(tree) %>% 
  collapse(node = 41)


# flip
ggtree(tree) %>% 
  flip(39, 41) # select the nodes to flip



##********************************************************
## 5. Adding and customizing labels ----
##********************************************************

# basic tip labels
ggtree(midpoint.root(tree)) + 
  geom_treescale(x = 0, y = 0, # x and y position of the treescale
                 width = 0.01) + # width of scale
  geom_tiplab(size = 4) + # displaying tip labels 
  coord_cartesian(clip = 'off')+ # allows us to draw outside the plot
  theme(plot.margin = margin(1,2,1,1, "cm")) # add space around the plot

# if we want to check tree tip labels
# this gives us an overview of what the tip labels will look like
# looks like the tip labels are the SRA numbers

head(tree$tip.label) 

# since we are planning on linking the metadata and blast file with
# the tree, we need to make sure they have a variable that links to the tip
# tip label names. the file_name variable will create this link
# we can also use a vector to check if there are any file_name
# observations that are not in the tree
metadata$file_name[!tree$tip.label %in% metadata$file_name]

# character(0) means they all match match up

# what if we want to switch out the tip label for something more
# meaningful to us?
# the new operator %<+% is discussed in section 7.1 of the book

ggtree(midpoint.root(tree)) %<+% # operator used to attach annotation data to tree
  metadata + # our metadata 
  geom_treescale(x = 0, y = 0,
                 width = 0.01)+ 
  coord_cartesian(clip = 'off')+
  theme(plot.margin = margin(1,2,1,1, "cm")) +
  geom_tiplab(aes(label = strain_ID)) # change the tip label to strain_ID


# lets make it easier to read
# aligning tip labels
gg_simple <- ggtree(midpoint.root(tree)) %<+%
  metadata + 
  geom_treescale(x = 0, y = 0,
                 width = 0.01)+ 
  coord_cartesian(clip = 'off')+
  theme(plot.margin = margin(1,2,1,1, "cm")) +
  geom_tiplab(aes(label = strain_ID), 
            color = "blue", # changing font color
            size = 4,  # changing font size
            offset = 0.01,
            align = TRUE) # this creates a dotted leader line



# lets add another label next to strain_ID

gg_simple +
  geom_tiplab(aes(label = serogroup), # add in serogroup information
            color = "black", 
            offset = 0.05,        # horizontal adjustment so the tiplabs don't overlap
            size = 4, 
            align = TRUE, # this creates the dotted leader line
            linetype = NA) # remove the dotted line b/w strain_ID and serogroup

# you can keep adding layer by layer and adjusting the offset each time
# takes a bit of playing around with offset


##********************************************************
## 6. Merging Heatmap with a Tree ----
##******************************************************** 

# first we want to tidy up the blast_raw data and 
# change it into a data frame using as.data.frame()
blast_df <- blast_raw %>% filter(percent_identical >= 80) %>% 
  pivot_wider(names_from = gene_name,
              values_from = percent_identical)%>% 
  relocate(ctxA, .before = ctxB) %>% as.data.frame()


# Next: set row names
# row names before are set as 1,2,3, ... etc
rownames(blast_df)

# changing the row names to file_name
rownames(blast_df) <- blast_df$file_name

# row names after 
rownames(blast_df)

# finally: remove the file_name column
# since file_name is now the same as rownames
# it is redundant and keeping it will interfere with our heatmap
# only want the gene names left as the columns

blast_df <- select(blast_df, !file_name) # ! means not

view(blast_df)

blast_df %>% head(5)


# gheatmap first look
# everything is over lapping!
gheatmap(gg_simple, # tree
         blast_df)   # heatmap


# lets tweak the formatting a bit
gheatmap(gg_simple, 
         blast_df, 
         offset = 0.1, # offset distance of heatmap from tree
         width = 1.5) # width of heatmap

# adjust the column labels
gheatmap(gg_simple, 
         blast_df, 
         offset = 0.1,
         width = 1.5,
         font.size = 5,
         colnames_position = "top") # can be top or bottom 
         
# check out help page for additional arguments        
?gheatmap

# specifying color using gheatmap
gheatmap(gg_simple, 
         blast_df, 
         offset = 0.1,
         width = 1.5,
         font.size = 5,
         colnames_position = "top", 
         color = "black", # color of cell border
         legend_title = "% identity", # title of legend
         low = "lightblue", # color of lowest value
         high = "darkblue") # color of highest value


# specifying color using ggplot2 scale_fill_ to override
# gheatmap fill, and saving as tree_heat
tree_heat <- gheatmap(gg_simple, 
         blast_df, 
         offset = 0.1,
         width = 1.5,
         font.size = 5,
         colnames_position = "top", 
         color = "black") +
  scale_fill_gradient(name = "% identity", # title of legend
                      low = "lightblue", high = "darkblue",
                      na.value = "grey77")

# warning message: Scale for 'fill' is already present. Adding another 
# scale for 'fill', which will replace the existing scale

# find more about colors and scales here: https://ggplot2-book.org/scale-colour.html

# final look at tree
tree_heat


##********************************************************
## 7. Exporting the final tree ----
##******************************************************** 

ggsave("output/tree_heat.jpeg", tree_heat, dpi = 300)

# can also specify height and width if a large tree
# open in paint and scale down


##********************************************************
## END
##******************************************************** 