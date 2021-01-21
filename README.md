# Diverstiy-with-Depth

Hello Seb and Valeriano

There is a Rmarkdown and a R script. I did not know which format you prefer to work. 

From the initial database: PA_df 

the script does:  1) Generic richness profile
                  2) NMDS
                  3) beta.pair per depth
                  4) betadisper (HERE THRE'S A PROBLEM)
                  5) Mantel tests.
                  
you'll see that each analysis re-starts from the PA_df to make it easier for you. 

I am confused with the betadisper per different depth. 

I can make it work for all depths and sites together considering groups as different depths. I obtain (1) the average distance to median ("b-dissimilarity" per depth ?), (2) anova results and (3) permutest pair-wise differences between (depths)

Still working with all depths, I can either work with the mother matrix ## 1 ##  or straight from beta.pair.abund matrix ## 2 ## 

However, the PROBLEM is making a betadisper per depth and check differences between islands. I think it is not possible as I need replicates... I obtain only 0 

- Valeriano: you said using all quadrats per depth. However, it's impossible if we work with the index occupancy-frequency which is F = nb Quadrats with genus / Total nb of quadrats. 
In other words, we are pulling all quadrats so it is not possible to work from the quadrats. 

For ease, you can check the word generated from the Rmarkdown with the outputs of the commands without running them.
