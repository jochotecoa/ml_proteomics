library(tibble); library(caret); library(ggplot2)

rf_rfe = '../output_rfe/split_by_sample/na_omit/rfProfile_repeatedcv_wout_seq_depth_mirna.rds' %>% readRDS
rfImp = varImp(rf_rfe)


rfImp[order(rfImp$Overall, decreasing = T)[1:10], , F] %>%
  rownames_to_column('name') %>% 
  arrange(Overall) %>%    # First sort by Overall. This sorts the dataframe but NOT the factor levels
  mutate(name=factor(name, levels=name)) %>%   # This trick updates the factor levels
  ggplot( aes(x=name, y=Overall)) +
  geom_segment( aes(xend=name, yend=0)) +
  geom_point( size=4, color="blue") +
  coord_flip() +
  theme_bw() +
  xlab("")
