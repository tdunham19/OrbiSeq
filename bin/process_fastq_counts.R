#!/usr/bin/env Rscript 

library (tidyverse)

# This script reads in a table of fastq counts and generates some plots and summaries
#
# Mark Stenglein Nov 9, 2023 

# ------------------------
# import fastq / bam read count info.   
# ------------------------
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  count_file_names = args
  output_directory="./"
} else {
  # if running via RStudio
  count_file_names = list.files(path = "../results/fastq_counts", pattern = "*count.txt$", full.names = T)
  output_directory="../results/"
}


read_count_file <- function(file_name) {
  # read the file 
  # wrapped with try here because empty files cause read.delim to fail/error
  df <- try(read.delim(file_name, header = F, stringsAsFactors = F))
  if (!inherits(df, 'try-error')) {

    colnames(df) <- c("sample_id", "count_type", "count")
    return(df)

  } else {
    # a vcf file with no variants return no dataframe
    print(paste0("no count data in file: ", file_name))
    return()
  }
}


# read in all the individual count files
# see: https://www.rdocumentation.org/packages/base/versions/3.3.1/topics/lapply?
df_list <- lapply(count_file_names, "read_count_file")

# this rbinds all the dfs together to make one big df
# see: https://stat.ethz.ch/pipermail/r-help/2007-April/130594.html
counts_df <- do.call("rbind", df_list)

# confirm that metadata exists for all datasets
dataset_names <- counts_df %>% 
  group_by(sample_id) %>% summarize() %>% pull(sample_id)

# put the counts in order
counts_df$count_type <- fct_relevel(counts_df$count_type, "initial")
counts_df$count_type <- fct_relevel(counts_df$count_type, "post_trimming", after=1)

write.table(counts_df, file=paste0(output_directory, "all_read_counts.txt"), sep="\t", row.names=F, col.names=T, quote=F)

count_types <- counts_df %>% group_by(count_type) %>% summarize() %>% pull(count_type)

if ("post_collapse" %in% count_types) {

   normalized_counts_df <- counts_df  %>% 
     pivot_wider(names_from = "count_type", values_from="count")  %>%
     mutate(post_collapse = post_collapse / initial,
            post_trimming = post_trimming / initial, 
            initial       = initial       / initial) %>%
     pivot_longer(cols=-sample_id, names_to = "count_type", values_to="count") 
} else {
normalized_counts_df <- counts_df  %>% 
  pivot_wider(names_from = "count_type", values_from="count")  %>%
  mutate(post_trimming = post_trimming / initial, 
         initial       = initial       / initial) %>%
  pivot_longer(cols=-sample_id, names_to = "count_type", values_to="count") 
}


# put the counts in order
normalized_counts_df$count_type <- fct_relevel(normalized_counts_df$count_type, "initial")
normalized_counts_df$count_type <- fct_relevel(normalized_counts_df$count_type, "post_trimming", after=1)

raw_counts_p <- ggplot(counts_df) +
 geom_line(aes(x=count_type, y=count, group=sample_id), size=0.25, linetype="dotted") +
 geom_boxplot(aes(x=count_type, y=count), fill="white", alpha=0.9) +
 theme_bw(base_size = 11) +
 scale_y_log10() + 
 ylab("Reads remaining") + 
 xlab("")

ggsave(paste0(output_directory, "preprocessing_read_counts.pdf"), raw_counts_p, units="in", width=7, height=9)

norm_counts_p <- ggplot(normalized_counts_df) +
 geom_line(aes(x=count_type, y=count, group=sample_id), size=0.25, linetype="dotted") +
 geom_boxplot(aes(x=count_type, y=count), fill="white", alpha=0.9) +
 theme_bw(base_size = 11) +
 ylab("Fraction reads remaining") + 
 xlab("")

ggsave(paste0(output_directory, "preprocessing_normalized_read_counts.pdf"), norm_counts_p, units="in", width=7, height=9)

