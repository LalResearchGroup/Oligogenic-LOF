---
title: "1. Genomic distribution of LOF-intolerant genes"
output: html_document
---

# **Identification and quantification of oligogenic loss-of-function disorders**

All 18,197 human genes were retrieved from the UCSC Human Genome Browser in hg19 assembly. We defined three intervals (0.5 Mb, 1 Mb, 5 Mb). To showcase the genomic distribution of LOF-intolerant genes analysis we will focus here on thee 5 Mb interval only. We investigate the distribution of LOF-intolerant genes following two approaches (1) around each gene and (2) as a non-overlapping sliding window across the genome. Using bedtools (www. bedtools.readthedocs.io), we annotate LOF-intolerant genes to these artificial intervals and count for each gene all neighboring genes and all LOF-intolerant genes per interval. LOF-intolerant genes were defined by a probability of loss-of-function (LOF) intolerance (pLI)>0.9 based on variation observed in 60,706 exomes of the Exome Aggregation Consortium (ExAC). We annotated a total of 3,230 pLI genes.

Input: Tab-delimited table with all human genes


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# **Packages**

```{r packagees, message = FALSE, warning = FALSE, echo = FALSE}
library("devtools")
library("tidyverse")
library("scales")
library("here")
library("gridExtra")
library("rmarkdown")
```

# **Read-in files**

```{r read-in files, message = FALSE, warning = FALSE}
# read-in files
## all genes
genes <- read.delim(here("data", "genes.txt"))
# make genes.txt file available in output/data folder
write_delim(genes, here("output", "data", "genes.txt"),
            delim = "\t")
# rename columns
genes <- rename(genes, chromosome = chr, start = chromStart, stop = chromEnd)
# pli genes
hpli_genes <- read.delim(here("data", "fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"))
# filter for gene and pli column
hpli_genes <- hpli_genes %>% 
  select(chr, cds_start, cds_end, pLI, gene) %>% 
  filter(pLI >= 0.9)
# set chr
hpli_genes$chr <- sub("^", "chr", hpli_genes$chr)
# save as tab delim .txt file for input in bedtools
write_delim(hpli_genes, here("output", "data", "pLI-genes.txt"),
                       delim = "\t")
# rename column
hpli_genes <- rename(hpli_genes, gene_hpli = gene)

# Create interval - 5Mb - 2.5Mb around each gene
genes_2500 <- genes %>%
  mutate(start = start - 2500000,
         stop = stop + 2500000) %>% # 18,197 genes
  filter(start >= 0 & stop >= 0) # 17,644 after exclusion of all negative values for start and stop position because bedtools cant compute them
# Write file for input in bedtools
write_delim(genes_2500, here("output", "data", "genes_5mb.txt"),
            col_names = FALSE, delim = "\t")
```

# ***intersect group_by with bedtools***

1. bedtools intersect –a genes_2500.txt –b pLI-genes.txt –wao >genes_2500intersect_hpli.txt
2. bedtools groupby –i genes_2500intersect_hpli.txt –g 1,2,3,4 –c 9 –o collapse >genes_2500master_hpli.txt

# **Read-in bedtool files**

```{r read-in bedtool files}
# read bedtools output all genes
genes_2500master_hpli <- read.delim(here("output", "data", "genes_2500master_hpli.txt"), 
                               header = FALSE)
# rename columns
genes_2500master_hpli <- genes_2500master_hpli %>% 
  rename(chromosome = V1, start = V2, stop = V3, gene = V4) %>% 
  arrange(gene)
```

# **pLI gene count per distinguished interval**

```{r N deletions with at least 1 additional high pLI gene affected - 5Mb}
# select columns needed
g2500 <- genes_2500master_hpli %>% 
  mutate(count = ifelse(V5 == ".", 0, str_count(V5, ",") + 1))

# rename columns
g2500 <- rename(g2500, genes_hpli_5mb_interval = V5, count_hpligenes = count)

# add column for "is hpli gene" differentiation
g2500 <- g2500 %>% 
  mutate(is_hpli_gene = ifelse(gene %in% hpli_genes$gene_hpli, 1, 0))

write_csv(g2500, here("output", "data", "g2500_hpli_counts.csv"))
```

```{r make main table for plotting pli}
# select columns needed
g2500_df <- g2500 %>%
  select(gene, is_hpli_gene, genes_hpli_5mb_interval, count_hpligenes, start, stop) %>% 
  mutate(number_additional_hpli_genes_affected = count_hpligenes - is_hpli_gene)
# make gene column factor
g2500_df$gene <- as.factor(as.character(g2500_df$gene))
# filter for only high pLI genes
g2500_df_hpli <- g2500_df %>% 
  filter(is_hpli_gene > 0)
```

```{r 5mb bar plot hpli}
plot1 <- 
  ggplot(data = g2500_df_hpli, aes(x = reorder(start, stop), 
                                   y = count_hpligenes)) +
  geom_bar(stat = "identity",
           width = 1,
           fill = "black") +
  theme_minimal() +
  ggtitle("Number of LOF-intolerant genes within 5 Mb interval around each gene\n2.5 Mb on each side", 
          subtitle = "") +
  ylab("Number of LOF-intolerant genes within interval") +
  xlab("LOF-intolerant genes") +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12,
                                  face = "bold"),
        plot.tag = element_text(size = 12,
                                face = "bold")) +
  scale_y_continuous(limits = c(0, 40),
                     breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40))
plot1

ggsave(here("output", "plots", "plot1.jpeg"), plot1, dpi = 500, width = 350, height = 300, units = "mm")
```

```{r make main table for plotting non pli}
# read in all genes table
g2500_df_all_genes <- read_csv(here("output", "data", "master_2500kb_all_genes.csv"))
# select columns needed
g2500_df_all_genes <- g2500_df_all_genes %>%
  select(gene, is_hpli_gene, genes_all_2500_interval, count_allgenes, start, stop)
  # mutate(number_additional_hpli_genes_affected = count_hpligenes - is_hpli_gene)
# make gene column factor
g2500_df_all_genes$gene <- as.factor(as.character(g2500_df_all_genes$gene))
# filter for only non pLI genes
g2500_df_non_hpli <- g2500_df_all_genes %>% 
  filter(is_hpli_gene == 0)
```

```{r 5mb bar plot non hpli}
plot1.1 <- 
  ggplot(data = g2500_df_non_hpli, aes(x = reorder(start, stop), 
                                   y = count_allgenes)) +
  geom_bar(stat = "identity",
           width = 1,
           fill = "black") +
  theme_minimal() +
  ggtitle("Number of all genes within 5 Mb interval around each gene\n2.5 Mb on each side", 
          subtitle = "") +
  ylab("Number of all genes within interval") +
  xlab("All genes") +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12,
                                  face = "bold"),
        plot.tag = element_text(size = 12,
                                face = "bold")) +
  scale_y_continuous(limits = c(0, 150),
                     breaks = c(0, 50, 100, 150))
plot1.1

ggsave(here("output", "plots", "plot1.1.jpeg"), plot1.1, dpi = 500, width = 350, height = 300, units = "mm")
```

```{r 5mb bar plot both hpli col}
plot1.2 <- 
  ggplot(data = g2500_df_all_genes, aes(x = reorder(start, stop), 
                                   y = count_allgenes,
                                   fill = is_hpli_gene)) +
  geom_bar(stat = "identity",
           width = 1) +
  theme_minimal() +
  ggtitle("Number of all genes within 5 Mb interval around each gene\n2.5 Mb on each side", 
          subtitle = "") +
  ylab("Number of all genes within interval") +
  xlab("All genes") +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 12,
                                  face = "bold"),
        plot.tag = element_text(size = 12,
                                face = "bold")) +
  scale_y_continuous(limits = c(0, 150),
                     breaks = c(0, 50, 100, 150))
plot1.2

ggsave(here("output", "plots", "plot1.2.jpeg"), plot1.2, dpi = 500, width = 350, height = 300, units = "mm")
```

**Density of LOF-intolerant genes within 5Mb**

```{r intervals 5mb}
df1 <- data.frame(chr = c("chr1"),
                 start = seq(from = 0, to = (249250621 - 5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 249250621, by = 5000000)
                 )

df2 <- data.frame(chr = c("chr2"),
                 start = seq(from = 0, to = (243199373 - 5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 243199373, by = 5000000)
                 )

df3 <- data.frame(chr = c("chr3"),
                 start = seq(from = 0, to = (198022430 - 5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 198022430, by = 5000000)
                 )

df4 <- data.frame(chr = c("chr4"),
                 start = seq(from = 0, to = (191154276 - 5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 191154276, by = 5000000)
                 )

df5 <- data.frame(chr = c("chr5"),
                 start = seq(from = 0, to = (180915260 - 5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 180915260, by = 5000000)
                 )

df6 <- data.frame(chr = c("chr6"),
                 start = seq(from = 0, to = (171115067 - 5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 171115067, by = 5000000)
                 )

df7 <- data.frame(chr = c("chr7"),
                 start = seq(from = 0, to = (159138663 - 5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 159138663, by = 5000000)
                 )

df8 <- data.frame(chr = c("chr8"),
                 start = seq(from = 0, to = (146364022-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 146364022, by = 5000000)
                 )

df9 <- data.frame(chr = c("chr9"),
                 start = seq(from = 0, to = (141213431-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 141213431, by = 5000000)
                 )

df10 <- data.frame(chr = c("chr10"),
                 start = seq(from = 0, to = (135534747-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 135534747, by = 5000000)
                 )

df11 <- data.frame(chr = c("chr11"),
                 start = seq(from = 0, to = (135006516-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 135006516, by = 5000000)
                 )

df12 <- data.frame(chr = c("chr12"),
                 start = seq(from = 0, to = (133851895-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 133851895, by = 5000000)
                 )

df13 <- data.frame(chr = c("chr13"),
                 start = seq(from = 0, to = (115169878-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 115169878, by = 5000000)
                 )

df14 <- data.frame(chr = c("chr14"),
                 start = seq(from = 0, to = (107349540-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 107349540, by = 5000000)
                 )

df15 <- data.frame(chr = c("chr15"),
                 start = seq(from = 0, to = (102531392-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 102531392, by = 5000000)
                 )

df16 <- data.frame(chr = c("chr16"),
                 start = seq(from = 0, to = (90354753-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 90354753, by = 5000000)
                 )

df17 <- data.frame(chr = c("chr17"),
                 start = seq(from = 0, to = (81195210-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 81195210, by = 5000000)
                 )

df18 <- data.frame(chr = c("chr18"),
                 start = seq(from = 0, to = (78077248-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 78077248, by = 5000000)
                 )

df19 <- data.frame(chr = c("chr19"),
                 start = seq(from = 0, to = (59128983-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 59128983, by = 5000000)
                 )

df20 <- data.frame(chr = c("chr20"),
                 start = seq(from = 0, to = (63025520-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 63025520, by = 5000000)
                 )

df21 <- data.frame(chr = c("chr21"),
                 start = seq(from = 0, to = (48129895-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 48129895, by = 5000000)
                 )

df22 <- data.frame(chr = c("chr22"),
                 start = seq(from = 0, to = (51304566-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 51304566, by = 5000000)
                 )

dfX <- data.frame(chr = c("chrX"),
                 start = seq(from = 0, to = (155270560-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 155270560, by = 5000000)
                 )

dfY <- data.frame(chr = c("chrY"),
                 start = seq(from = 0, to = (59373566-5000000), by = 5000000),
                 stop = seq(from = 5000000, to = 59373566, by = 5000000)
                 )
# # bind 5mb
df.5 <- bind_rows(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19, df20, df21, df22, dfX, dfY)
# write .txt for bedtools input 5mb
write_delim(df.5, here("output", "data", "intervals5.txt"),
                   col_names = FALSE, delim = "\t")
```

**Perform *intersect* and *group_by* with bedtools**

1. bedtools intersect –a intervals5.txt –b pLI-genes.txt –wao >intervals_5intersect_hpli.txt
2. bedtools groupby –i intervals_5intersect_hpli.txt –g 1,2,3,4 –c 8 –o collapse >intervals_5master_hpli.txt

```{r read-in bedtool files for density analysis}
# read bedtools output hpli genes
df.interv5.hpli <- read.delim(here("output", "data", "intervals_5master_hpli.txt"), 
                               header = FALSE)

# rename columns
df.interv5.hpli <- df.interv5.hpli %>%
  select(V1, V2, V3, V5) %>% 
  dplyr::rename(chromosome = V1, start = V2, stop = V3, gene = V5) %>% 
  arrange(gene)
```

**Gene count per distinguished interval**

```{r Gene count per distinguished interval, hpli - 5mb}
# select columns needed
df.5.hpligenes <- df.interv5.hpli %>% 
  mutate(count = ifelse(gene == ".", 0, str_count(gene, ",") + 1))

# rename columns
df.5.hpligenes <- dplyr::rename(df.5.hpligenes, genes = gene, hpli_genes_count_5 = count)

write_csv(df.5.hpligenes, here("output", "data", "df.5_hpligenes_count.csv"))
```

**Make master df - 5Mb**

```{r df 5Mb}
d5 <- df.5.hpligenes
# make chr factor
d5$chromosome <- factor(d5$chromosome, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'))
# order df
d5 <- with(d5, d5[order(chromosome, start),])
# add ID
d5 <- tibble::rowid_to_column(d5, "ID")

# d5[is.na(d5)] <- 0
```

```{r density plot 5Mb}
plot2 <- 
  ggplot(data=d5, aes(x=start, y=1)) +
  facet_grid(chromosome ~ ., switch='y') +
  geom_tile(aes(fill=hpli_genes_count_5)) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y.left = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x=element_blank(),
        axis.text.x=element_text(),
        axis.line = element_blank(),
        strip.background = element_rect(color = "NA"),
        legend.title = element_text(face = "bold"),
        plot.tag = element_text(size = 18,
                                face = "bold"),
        plot.title = element_text(face = "bold"),
        title = element_text(face = "bold", 
                             hjust = 0.5),
        legend.justification = c(.9, .001),
        legend.position = c(1, .4),
        legend.text = element_text()) +
  labs(fill = "Number\nof genes") +
  scale_x_continuous(limits = c(0, 250000000),
                   breaks = seq(0, 250000000, 50000000),
                   labels = c("0Mb", "50Mb", "100Mb", "150Mb", "2000Mb", "250Mb")) +
  scale_fill_gradientn(colours = c("white", "#2ca25f"),
                       breaks = c(0, 15, 30)) +
  ggtitle("Density of LOF-intolerant genes within 5 Mb")
plot2

ggsave(here("output", "plots", "plot2.jpeg"), plot2, dpi = 500, width = 350, height = 300, units = "mm")
```
