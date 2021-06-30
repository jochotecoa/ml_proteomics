forceLibrary <- function(list.of.packages, ...) {
  checkNewPackages <- function(list.of.packages) {
    new.packages.log = !(list.of.packages %in% installed.packages()[,"Package"])
    new.packages <- list.of.packages[new.packages.log]
    return(new.packages)
  }
  new.packages = checkNewPackages(list.of.packages)
  if (length(new.packages)) {
    print(paste('Trying to install the following packages:', paste(new.packages)))
    install.packages(new.packages, ...)
    new.packages = checkNewPackages(list.of.packages)
    if (length(new.packages)) {
      print(paste(paste(new.packages), 'were not installed through the easy way'))
      print("Let's try the hard way then")
      setRepositories(graphics = F, ind = 1:8)
      install.packages(new.packages, ...)
      new.packages = checkNewPackages(list.of.packages)
      if (length(new.packages)) {
        stop('forceLibrary was not able to install the following packages: ', 
             paste(new.packages))
      }
    }
  } 
  
  lapply(list.of.packages, library, character.only = T)
  
  invisible()
}

cleanProtIds = function(protein_table) {
  colprotids = grepl('\\|', protein_table[1, ])
  
  if (!any(colprotids)) {
    protein_table = protein_table %>% 
      rownames_to_column()
    colprotids = 'rowname'
  }
  # Filter out double IDs
  
  protein_table = protein_table[!grepl(protein_table[, colprotids], pattern = ':'), ]
  if (class(protein_table) != "data.frame") {
    protein_table = as.data.frame(protein_table)
  }
  names = strsplit(as.character(protein_table[, colprotids]), '\\|')
  names = as.character(lapply(names, '[', 2))
  protein_table$uniprot_gn = names
  return(protein_table)
}

naToZero = function(x) {
  x[is.na(x)] = 0
  return(x)
}

zeroToNa = function(x) {
  x[x == 0] = NA
  return(x)
}

pseudocount = function(x, addition = 1) {
  x = x + addition
  return(x)
}


transcrToGene = function(table, aggregate = F, prot_cod = F, ...) {
  forceLibrary(c('biomaRt', 'dplyr'))
  
  # sampl = table[nrow(table), ]
  enst_col = apply(table, 2, grepl, pattern = 'ENST') %>% 
    apply(2, any) %>%
    as.logical()
  if (sum(enst_col) == 0) {
    table[, 'rownames'] = rownames(table)
    enst.rown = grepl(pattern = 'ENST', x = table[, 'rownames'])
    if (!sum(enst.rown)) {stop(print(table[1, ]))}
    enst_col = grepl(pattern = 'rownames', x = colnames(table))
  }
  # If integer or else, it might try to aggregate it
  table[, enst_col] = as.character(table[, enst_col])
  version = grepl('\\.', table[, enst_col]) %>% sum()
  if (version) {version = T} else {version = F}
  if (version) {
    transcript_id = 'ensembl_transcript_id_version'
  } else {
    transcript_id = 'ensembl_transcript_id'
  }
  
  values = table[, enst_col]
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org') 
  if (prot_cod) {
    transcr_biotypes = getBM(attributes = c(transcript_id, 'transcript_biotype'), 
                             filters = transcript_id, values = values, 
                             mart = mart.human)
    isProtCod = transcr_biotypes$transcript_biotype == 'protein_coding'
    values = values[isProtCod]
    print(paste0(sum(!isProtCod), ' transcripts were not protein_coding'))
  }
  
  new_cols = getBM(attributes = c(transcript_id, 'ensembl_gene_id'), 
                   filters = transcript_id, values = values, mart = mart.human)
  
  table = merge.data.frame(x = table, y = new_cols, 
                           by.x = colnames(table)[enst_col], 
                           by.y = transcript_id, ...)
  if (nrow(table) < length(values)) {
    print(paste0(length(values) - nrow(table), 
                 ' transcripts had no gene matched'))
  }
  if (prot_cod & !aggregate) {
    table = tibble::rownames_to_column(table)
  }
  if (aggregate) {
    int_cols = grepl('integer', sapply(X = table[1, ], FUN = typeof))
    int_cols = int_cols + grepl('double', sapply(X = table[1, ], 
                                                 FUN = typeof))
    int_cols = as.logical(int_cols)
    table = aggregate(x = table[, int_cols], by = list(table$ensembl_gene_id), 
                      FUN = sum, na.rm = T)
    colnames(table)[1] = 'ensembl_gene_id'
  }
  return(table)
}

rmMirnas = function(x) {
  mirna.cols = grep(pattern = 'hsa', x = colnames(x))
  y = x[, -mirna.cols]
  return(y)
}

forceSetWd = function(x) {
  if (dir.exists(x)) {
    setwd(x)
  } else {
    dir.create(x)
    if (dir.exists(x)) {
      setwd(x)
    } else {
      warning(c('Warning: ', x, 
                ' could not be created as a dir due to permission issues'))
    }
  }
}

mergeFiles = function(files_patt =  '.', by_col = 'Name', 
                      row_names = F, progr_bar = T, path = path, all_true = F, 
                      recursive = F, ...) {
  if (progr_bar) {
    forceLibrary('pbmcapply')
  }
  forceLibrary('dplyr')
  files = list.files(pattern = files_patt, recursive = recursive, path = path, 
                     full.names = T, include.dirs = F)
  # files = files[!grepl('total', files)]
  # files = files[-1]
  print(paste('Number of files found:', length(files)))
  file = files[1]
  stopifnot(file.exists(file))
  voom_file = read.table(file, stringsAsFactors = F, ...)
  if (row_names) {
    voom_file = voom_file %>% tibble::rownames_to_column() %>% 
      dplyr::select(rowname, everything())
    by_col = 'rowname'
  }
  colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
  big_quant_voom = voom_file
  if (progr_bar) {
    pb = progressBar(max = length(files[-1]))
  } else {
    p = progress_estimated(length(files[-1]))
  }
  for (file in files[-1]) {
    voom_file = read.table(file, stringsAsFactors = F, ...)
    if (row_names) {
      voom_file = voom_file %>% tibble::rownames_to_column() %>% 
        dplyr::select(rowname, everything())
    }
    colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
    big_quant_voom = merge.data.frame(big_quant_voom, voom_file, by = by_col, all = all_true)
    if (progr_bar) {
      setTxtProgressBar(pb, grep(file, files[-1]))
    } else {
      p$tick()$print()
    }
  }
  if (progr_bar) {
    close(pb)
  }
  return(big_quant_voom)
} 

openMart2018 <- function(...) {
  forceLibrary('biomaRt')
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org', ...) 
}

filterProtCod = function(table) {
  forceLibrary('biomaRt')
  
  sampl = table[nrow(table), ]
  enst_col = grep(pattern = 'ENST', x = sampl)
  if (length(enst_col) > 1) {enst_col = enst_col[1]}
  if (length(enst_col) == 0) {
    table[, 'rownames'] = rownames(table)
    sampl = table[nrow(table), ]
    enst_col = grep(pattern = 'ENST', x = sampl)[1]
    if (is.na(enst_col)) {print(sampl)}
  }
  version = grepl('\\.', sampl[, enst_col])
  if (length(version) == 0) {version = F}
  if (version) {
    transcript_id = 'ensembl_transcript_id_version'
  } else {
    transcript_id = 'ensembl_transcript_id'
  }
  
  values = table[, enst_col]
  mart.human = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                       dataset = 'hsapiens_gene_ensembl',
                       host = 'http://apr2018.archive.ensembl.org') 
  
  transcr_biotypes = getBM(attributes = c(transcript_id, 'transcript_biotype'), 
                           filters = transcript_id, values = values, 
                           mart = mart.human)
  isProtCod = transcr_biotypes$transcript_biotype == 'protein_coding'
  values = values[isProtCod]
  print(paste0(sum(!isProtCod), ' transcripts were not protein_coding'))
  new_table = table[table[, enst_col] %in% values, ]
  return(new_table)
}
filterSamplesBySeqDepth = function(df) {
  seq_depth_ratio <- df %>% 
    colSums(na.rm = T) %>% 
    `/` (mean(.)) %>% 
    log2() 
  
  if (!all(seq_depth_ratio > -2)) {
    warning(sum(!(seq_depth_ratio > -2)), 
            ' sample(s) filtered out due to sequencing depth ', 
            names(seq_depth_ratio) [seq_depth_ratio <= -2], immediate. = T)
  }
  
  df = df[, seq_depth_ratio > -2]
}
apply_2D = function(df, FUN, col.x = NULL, col.y = NULL, complete_cases = NULL, 
                    y = NULL, ...) {
  library(dplyr)
  if (exists('result.df')) {
    rm(result.df)
  }
  if (is.null(complete_cases)) {
    complete_cases = length(col.x)
  }
  if (length(col.x) != length(col.y)) {
    stop('The number of columns is different for each group')
  }
  p = progress_estimated(nrow(df))
  for (row in rownames(df)) {
    if (!is.null(y)) {
      X = df[row, ] %>% as.numeric()
      row.y = grep(row, rownames(df))[1]
      Y = y[row.y, ] %>% as.numeric()
    } else {
      X = df[row, col.x] %>% as.numeric()
      Y = df[row, col.y] %>% as.numeric()
    }
    if (sum(complete.cases(X, Y)) >= complete_cases) {
      res = FUN(x = X, y = Y, ...)
      result.ul = res %>% unlist()
      if (exists('result.df')) {
        result.df = result.ul %>% unlist() %>% t() %>% as.data.frame() %>% 
          rbind(result.df, .)
        rownames(result.df)[nrow(result.df)] = row
      } else {
        result.df = result.ul %>% unlist() %>% t() %>% as.data.frame()
        rownames(result.df) = row
      }
    } 
    p$tick()$print()
  }
  return(result.df)
}

mergeFilesRds = function(files_patt = NULL, by_col = 'Name', 
                      row_names = F, progr_bar = T, path = ., ...) {
  if (progr_bar) {
    forceLibrary('pbmcapply')
  }
  forceLibrary('dplyr')
  files = list.files(path = path, pattern = files_patt, recursive = T, 
                     full.names = T)
  # files = files[!grepl('total', files)]
  # files = files[-1]
  print(paste('Number of files found:', length(files)))
  file = files[1]
  stopifnot(file.exists(file))
  voom_file = readRDS(file)
  if (row_names) {
    voom_file = voom_file %>% tibble::rownames_to_column() %>% 
      dplyr::select(rowname, everything())
    by_col = 'rowname'
  }
  colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
  big_quant_voom = voom_file
  if (progr_bar) {
    pb = progressBar(max = length(files[-1]))
  } else {
    p = progress_estimated(length(files[-1]))
  }
  for (file in files[-1]) {
    voom_file = readRDS(file)
    if (row_names) {
      voom_file = voom_file %>% tibble::rownames_to_column() %>% 
        dplyr::select(rowname, everything())
    }
    colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
    big_quant_voom = merge.data.frame(big_quant_voom, voom_file, by = by_col, ...)
    if (progr_bar) {
      setTxtProgressBar(pb, grep(file, files[-1]))
    } else {
      p$tick()$print()
    }
  }
  if (progr_bar) {
    close(pb)
  }
  return(big_quant_voom)
} 

iq = function(x, na.rm = F) {
	if (na.rm) {x = na.omit(x)}
	first_q = as.numeric(quantile(x, 0.25))
	third_q = as.numeric(quantile(x, 0.75))
	y = third_q - first_q
	return(y)
}

outlier = function(x, y=NULL, na.rm = F, onlyExtreme = F) {
	outliers = NULL
	if (na.rm) {x = na.omit(x)}
	q1 = quantile(x, 0.25)
	q3 = quantile(x, 0.75)
	
	if (!is.null(y)) {
		isOutlier = y < q1 - 1.5*iq(x) | y > q3 + 1.5*iq(x)
		if (onlyExtreme) {
			isOutlier = y < q1 - 3*iq(x) | y > q3 + 3*iq(x)
		}
		return(as.logical(isOutlier))
	}
	
	for (i in x) {
		isOutlier = i < q1 - 1.5*iq(x) | i > q3 + 1.5*iq(x)
		if (onlyExtreme) {
			isOutlier = i < q1 - 3*iq(x) | i > q3 + 3*iq(x)
		}
		if (isOutlier) {outliers = c(outliers, i)}
	}
	return(outliers)
}

shift_median = function(df, median_of_medians) {
	for (i in colnames(df)) {
		sample = df[, i]
		corr_factor = median_of_medians - median(sample, na.rm = T)
		df[, i] = df[, i] + corr_factor
	}
	return(df)
}

normalizeProteomics = function(df) {
  if (df %>% apply(2, is.character) %>% any()) {
    stop('at least one character column')
  }
  common_set = na.omit(df)
	medians = apply(common_set, 2, median)
	median_of_medians = median(medians)
	df = shift_median(df, median_of_medians)
	return(df)
}

rmDuplicatedColumns <- function(df) {
  df = df[!duplicated(as.list(df))]
  return(df)
}

cleanCircNames <- function(x) {
  x %>% 
    strsplit('\\|') %>%
    sapply('[[', 1) 
}

iq_ratio <- function(x, prba = F) {
  y = x #%>% 
  # dplyr::select(matches(paste0(col, collapse = '.*')))
  
  # colname = paste0('iq_ratio_', paste0(col, collapse = '_'), collapse = '_')
  # x[, colname] = NA
  ratio_vector = NULL
  if (prba) {  pb_514 <- progressBar(max = nrow(x))}
  for (rwo in rownames(x)) {
    expr_values = y[rwo,]
    expr_values = expr_values %>% unlist() %>% as.numeric()
    expr_values = naToZero(expr_values)
    quantiles = quantile(x = expr_values, probs = c(1/7, 0.5, 7/8))
    ratio = (quantiles[3]-quantiles[1])/quantiles[2]
    ratio = as.numeric(ratio)
    ratio_vector = c(ratio_vector, ratio)
    # x[rwo, colname] = xa
    if (prba) { 
      i = grep(rwo, rownames(x))
      setTxtProgressBar(pb_514, i)
      if (i == nrow(x)) {
        close(pb_514)
      }
    }
  }
  # x[is.na(x[, colname]), colname] = 1000
  ratio_vector[is.na(ratio_vector)] = 1000
  return(ratio_vector)
}

addVarsProt <- function(x, fnc_list, by_str) {
  fnc_str = fnc_list[1]
  fnc = get(fnc_str)
  transf_cols = colnames(x)[!grepl(by_str, colnames(x))]
  by_lst = x[, by_str] %>% list()
  names(by_lst) = by_str
  df = x %>% 
    .[, !grepl(by_str, colnames(x)), F] %>% 
    aggregate.data.frame(by = by_lst, FUN = fnc, na.rm = TRUE) 
  colnames(df)[!grepl(by_str, colnames(df))] = paste0(transf_cols, '_', fnc_str)
  if (length(fnc_list) > 1) {
    for (fnc_str in fnc_list[-1]) {
      fnc = get(fnc_str)
      transf_cols = colnames(x)[!grepl(by_str, colnames(x))]
      by_lst = x[, by_str] %>% list()
      names(by_lst) = by_str
      x_2 = x %>% 
        .[, !grepl(by_str, colnames(x)), F] %>% 
        aggregate.data.frame(by = by_lst, FUN = fnc, na.rm = TRUE) 
      colnames(x_2)[!grepl(by_str, colnames(x_2))] = paste0(transf_cols, '_', fnc_str)
      df = x_2 %>% 
        merge.data.frame(x = df, by = by_str)
    }
  }
  return(df)
}

getSubOptVars <- function(results) {
  rmses = results$RMSE
  vars = results$Variables
  
  opt_error = min(rmses) + diff(range(rmses)) /10
  opt_vars = results$Variables[results$RMSE < opt_error][1]
  return(opt_vars)
}

mergeFilesCsv = function(files_patt =  '.', by_col = 'Name', 
                      row_names = F, progr_bar = T, path = path, all_true = F, 
                      recursive = F, ...) {
  if (progr_bar) {
    forceLibrary('pbmcapply')
  }
  forceLibrary('dplyr')
  files = list.files(pattern = files_patt, recursive = recursive, path = path, 
                     full.names = T, include.dirs = F)
  # files = files[!grepl('total', files)]
  # files = files[-1]
  print(paste('Number of files found:', length(files)))
  file = files[1]
  stopifnot(file.exists(file))
  voom_file = read.csv(file, stringsAsFactors = F, ...)
  if (row_names) {
    voom_file = voom_file %>% tibble::rownames_to_column() %>% 
      dplyr::select(rowname, everything())
    by_col = 'rowname'
  }
  colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
  big_quant_voom = voom_file
  if (progr_bar) {
    pb = progressBar(max = length(files[-1]))
  } else {
    p = progress_estimated(length(files[-1]))
  }
  for (file in files[-1]) {
    voom_file = read.csv(file, stringsAsFactors = F, ...)
    if (row_names) {
      voom_file = voom_file %>% tibble::rownames_to_column() %>% 
        dplyr::select(rowname, everything())
    }
    colnames(voom_file)[-1] = paste(colnames(voom_file)[-1], file, sep = '_')
    big_quant_voom = merge.data.frame(big_quant_voom, voom_file, by = by_col, all = all_true)
    if (progr_bar) {
      setTxtProgressBar(pb, grep(file, files[-1]))
    } else {
      p$tick()$print()
    }
  }
  if (progr_bar) {
    close(pb)
  }
  return(big_quant_voom)
} 
