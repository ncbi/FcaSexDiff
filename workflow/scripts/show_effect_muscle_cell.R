library(tidyverse)
library(ggrepel)
library(ggthemes)
#library(xlsx)
library(writexl)


process <- function(resol, do_append) {
with_muscle_file = paste0("exports/extras_without_sex_specific/scatter_plots/resol~", resol, "/scatter_count_bias_body_stringent_", resol, ".csv")

with_muscle = (
    read.csv(with_muscle_file)
    %>% mutate(type = "with_muscle_cells")
)

print(with_muscle)

without_muscle_file = paste0("exports/extras/scatter_plots/resol~", resol, "/scatter_count_bias_body_stringent_", resol, ".csv")



without_muscle = (
    read.csv(without_muscle_file)
    %>% mutate(type = "without_muscle_cells")
)

print(without_muscle)


df = (
  bind_rows(with_muscle, without_muscle)
  %>% mutate(xmin = pct_male - sd_male)
  %>% mutate(xmax = pct_male + sd_male)
  %>% mutate(ymin = pct_female - sd_female)
  %>% mutate(ymax = pct_female + sd_female)
  %>% mutate(label = ifelse(((pct_female > 2) | (count_bias_type == "Female")) | 
                            ((pct_male > 2) | (count_bias_type == "Male")),
                            as.character(cluster), ''))
  %>% arrange(cluster, type)
)


base <- (
    ggplot(df, aes(y=pct_female, x=pct_male, color=count_bias_type))
    + geom_abline(intercept = 0, slope = 2, linetype='dotted')
    + geom_abline(intercept = 0, slope = 0.5, linetype='dotted')
    + geom_abline(intercept = 0, slope = 1, linetype='solid', color="#008800")
    + geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.5, size=0.5, color="gray")
    + geom_errorbarh(aes(xmin=xmin, xmax=xmax), height=0.5, size=0.5, color="gray")
    + geom_text_repel(aes(label=label), min.segment.length = 0, max.overlaps=Inf)
    + geom_point(size=1)
    + facet_wrap(~type)
    + theme_few()
    + scale_color_manual(values=c("Female"="red", "Male"="blue", "Unbiased"="black"))
)

main <- (
    base
    + theme(legend.position = "bottom")
    + labs(
        y='% female cells (replicate average)',
        x='% male cells (replicate average)',
        title = resol
    )
)

print(main)

#write.xlsx(df, "muscle_effect.xlsx", sheetName=resol, 
#  col.names=TRUE, row.names=FALSE, append=do_append)

df %>% select(cluster, major_annotation, type, count_bias_type, count_bias_padj, pct_male, pct_female, sd_male, sd_female)

}

df_list = list()

resolutions = c(
  "L0.4",
  "L0.6",
  "L0.8",
  "L1.0",
  "L1.2",
  "L1.4",
  "L1.6",
  "L1.8",
  "L2.0",
  "L4.0",
  "L6.0",
  "L8.0",
  "L10.0"
)

pdf("muscle_effect.pdf",  height=11, width=20)

df_list[["annotation"]] <- process("annotation", FALSE)

for  (resol in resolutions) {
  df_list[[resol]] <- process(resol, TRUE)
}

write_xlsx(df_list, "muscle_effect.xlsx")
