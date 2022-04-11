library(tidyverse)
library(ComplexUpset)
library(svglite)
#source('../colors.R')

infile = 'all_biased_gene.tsv'
#svgfile = '../../data/anova/anova_upset.svg'

df = (
  read.delim(infile) %>%
  filter(cluster != "unannotated") %>%
  mutate(bias = ifelse(bias > 0, "female", ifelse(bias < 0, "male", "unbiased"))) %>%
  group_by(cluster, bias, tissue) %>% 
  count() %>%
  ungroup()
)
head(df)
dim(df)

print(df %>% filter(cluster == "hemocyte") %>% mutate(logn = log(n)))

dim(unique(df[, c('cluster', 'tissue')]))

df = df %>% spread('tissue', 'n', fill=0) %>% gather('tissue', 'n', -cluster, -bias)
df = df %>% spread('bias', 'n', fill=0) %>% gather('bias', 'n', -cluster, -tissue)

dim(df)
print(length(unique(df$tissue)))
print(length(unique(df$cluster)))


df = df %>% mutate(cluster = fct_reorder(cluster, n, sum, .desc=TRUE))  %>%
     mutate(logn = ifelse(n > 0, log(n), 0)) %>%
     mutate(y = n)

print(df %>% filter(cluster == "hemocyte"))

pdf("Rplots.pdf", width=20, height=20)

#ggplot(df, aes(x=cluster, fill=tissue)) +
#  geom_bar(data = subset(df, bias=="female"), aes(y=y), position="stack", stat="identity") +
#  geom_bar(data = subset(df, bias=="male"), aes(y=-y), position="stack", stat="identity") +
#  #scale_y_symmetric(labels = abs, trans=function(x) {ifelse(x < 0, -log2(-x), log2(x))})+
#  #scale_y_continuous(trans='log2') +
#  #coord_trans(y="log2") +
#  coord_flip()

library(ggpol)

ggplot(df, aes(x = cluster, y = n, fill = tissue)) +
  geom_bar(stat = "identity", position='stack') +
  facet_share(~bias, dir = "h", scales = "free", reverse_num = TRUE) +   # note: scales = "free"
  coord_flip() +
  theme_minimal() +
  labs(y = "Count", x = "num biased", title = " ")

#
#  write.table("sex_biased_gene_distributions.tsv", sep = "\t",  row.names = FALSE)


quit()

anova_res = (
    read.table(infile, row.names=1, header=T)
    %>% select(status)
    %>% separate(status, c('Karyotype', 'TRA', 'Interaction'), sep=1:3, remove=F)
    %>% mutate(across(2:4, as.logical))
    %>% arrange(status)
)
print(head(anova_res))


#svglite(svgfile, width=6, height=4)

upset(
    anova_res,
    # sets are named bottom to top
    rev(c('Karyotype', 'TRA', 'Interaction')),
    name = 'ANOVA groups',
    set_sizes = F,
    # make sure the sets are not converted to integer
    # otherwise the columns are set improperly
    encode_sets = F,
    # required to make the sets order intact
    sort_sets = F,
    sort_intersections = FALSE,
    intersections = list(
        'Outside of known sets', # empty intersection
        'Interaction',
        'TRA',
        c('Interaction', 'TRA'),
        'Karyotype',
        c('Karyotype', 'Interaction'),
        c('Karyotype', 'TRA'),
        c('Karyotype', 'Interaction', 'TRA')
    ),
    base_annotations = list(
        '# transcripts in ANOVA group' = intersection_size(
            counts = TRUE,
            mapping = aes(fill = status)
        )
        + scale_fill_manual(values = color_values)
        + theme(legend.position = "none")
    )
)

