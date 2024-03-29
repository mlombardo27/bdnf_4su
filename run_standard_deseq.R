run_standard_deseq = function(folder_of_featurecounts, 
                               base_grep = "Ctrl",
                               contrast_grep = "TDPKD", 
                              grep_pattern = "",
                              suffix = ".Aligned.sorted.out.bam",
                              baseName = "control",
                              contrastName = 'TDPKD'
                              ){
    library(DESeq2)
    library(data.table)
    library(tidyverse)
    # First we are going to load in the functions that I've written as helper scripts
    create_feature_path = here::here('create_feature_count_table.R')
    make_deseq_path = here::here('make_deseq_dfs.R')
    #you'll want to adjust the file paths accordingly
    #source will bring the functions in these Rscripts into the current environment
    source(create_feature_path)
    source(make_deseq_path)
    # this function creates a table wtih the first column being Geneid, then each
    #column being the count in a sample, and the final column being the gene name
    #it takes as an input the folder path where all the feature_counts files are
    #the prefix is somethign which it appended by feature counts and the suffix
    #is the bam suffix. Typically the bam suffix if it was aligned through
    #our most recent version of the pipeline will be .Aligned.sorted.out.bam
    #this table will be the input to our next function
    timecourse_feature = create_feature_count_table(folder_of_featurecounts,suffix = suffix)
    

    #but we don't want that for now
    timecourse_counts = make_deseq_dfs(timecourse_feature,
                                       grep_pattern = grep_pattern,
                                       base_grep = base_grep, 
                                       contrast_grep = contrast_grep)$conv_df #note the $ we're only take one item of the list
    
    
    timecourse_meta = make_deseq_dfs(timecourse_feature, 
                                     grep_pattern = grep_pattern,
                                     base_grep = base_grep, 
                                     contrast_grep = contrast_grep)$coldata #note the $ we're only take one item of the list
    
    
    # When you do DESeq2, you should remove genes with very low counts,
    #it can mess up the way DESeq2 normalizes counts and give you funky p-valeus
    #here I wrote another helper function to do it
    timecourse_counts = filter_count_table(timecourse_counts)
    
    #This helper function will create another column using
    #the 'cond' column of the metadata DF and give it whatever name you find
    #more meaningful than 'base' and 'contrast'
    #it also turns it into a factor with the 'base' condition as the first level
    #the new column will be named 'comparison'

    
    timecourse_meta = rename_relevel_for_deseq(timecourse_meta,
                                               baseName = baseName,
                                               contrastName = contrastName)
    
    

    
    #now that we've done all the data sorting we actually get to the DESeq2 part
    #first we make an object using the counts, the meta_data, and the column we
    #want to compare. You should read the DESeq2 manual for
    #more complicated design explanations
    dds_timecourse = DESeqDataSetFromMatrix(timecourse_counts,
                                            colData = timecourse_meta, 
                                            design = ~ comparison)
    #that just created the object
    #to actually run the analysis we just call "DESeq"
    dds_timecourse = DESeq(dds_timecourse)
    
    #we can quickly view the results with the function results()
    results(dds_timecourse)
    #and then we can get a summary with summary on results
    results(dds_timecourse) %>% summary()
    #we can pull the results data frame out and use the feature_counts table to 
    #append the gene_names back on like this

    if("gene_name" %in% colnames(timecourse_feature)){
        res_timecourse = results(dds_timecourse) %>% 
            as.data.frame() %>% 
            rownames_to_column('Geneid') %>% 
            left_join(timecourse_feature %>% dplyr::select(Geneid,gene_name)) %>% 
            as.data.table()
    }else{
        res_timecourse = results(dds_timecourse) %>% 
            as.data.frame() %>% 
            rownames_to_column('Geneid') %>% 
            separate(Geneid,"ensgene") %>% 
            left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol)) %>% 
            as.data.table() %>% 
            dplyr::rename(gene_name = symbol)
    }

    
    
    return_list = list(dds_timecourse,res_timecourse)
    names(return_list) = c("deseq_obj","results_table")
    
    return(return_list)
    
}
