## Facet Plot Function 
FacetPower <- function(transform.data, col.x, col.x.title, no.facet.flag = FALSE,plot.type = F,alpha = 0.05, iter.num = 1000) {
  plot.gg <- ggplot(transform.data, 
                    aes_string(x = col.x,
                               y = 'Power', 
                               col = 'Tests', 
                               linetype = 'Tests')) +
    geom_path(size = 0.95) + 
    geom_point(size = 1.2) + 
    #geom_errorbar(aes_string(ymax = 'Up' , ymin = 'Down' ,x = col.name.x , y = 'Power' , col = 'Tests'), width = 0.01) + 
    geom_hline(yintercept = 0.05, col = 'grey', lwd = 1.05, linetype = 'dotted') + 
    theme_light() + 
    # scale_colour_brewer(palette = "Dark2", labels = c("SD", "CQ", "Lopes", "SH.Ln.LB", "SH.Sn.LB","Simes")) +
    scale_linetype_manual("Tests",
                          breaks = c('SD', 'CQ', 'SH.Sn.LB', 'SH.Ln.LB', 'Simes'), 
                          values = c('solid', 'solid',  'longdash', 'longdash', 'solid'),
                          labels = c("SD", "CQ",  expression(SH[n/4]), expression(SH[n/2]),"Simes")) +
    guides(colour = guide_legend(override.aes = list(size = 0.85))) + theme_bw() + 
    # scale_colour_grey(start = 0.2, end = 1, labels = c("SD", "CQ", "SH.Ln.LB", "SH.Sn.LB","Simes")) + 
    scale_color_manual(values = c("#B15928", "#E31A1C", "#6A3D9A", "#CAB2D6",  "#33A02C"),
                       breaks = c('SD', 'CQ', 'SH.Sn.LB', 'SH.Ln.LB', 'Simes'), 
                       labels = c("SD", "CQ",  expression(SH[n/4]),expression(SH[n/2]),"Simes")) +
    theme(legend.position="bottom", strip.text = element_text(size = 12), text = element_text(size = 12),
          legend.key.height = unit(0.5, "in"), legend.key.width = unit(0.65, "in")) 
  if (!no.facet.flag) {
    return(plot.gg + facet_grid( Obs ~  False.Hypotheses, scales = 'free'))
  } else { 
    return(plot.gg + facet_grid(Obs ~ .))
  }
}


### Loading simulations power non-equal covariance matrices 
setwd(paste0(here::here(), '/Non-equal Covaraince/Power'))

dat.1 <- read.csv('nonequal.power.simulation.matrix.1.csv') ## AR(1)
dat.3 <- read.csv('nonequal.power.simulation.matrix.3.csv') ## Block Covariance 
dat.5 <- read.csv('nonequal.power.simulation.matrix.5.csv') ## Cai model 
dat.7 <- read.csv('nonequal.power.simulation.matrix.7.csv') ## Equicorrelated 



### Directory to save plots 
setwd(paste0(here::here(), '/Non-equal Covaraince/Images'))

# Plotting Setting 1 - AR(1) ----------------------------------------------

dat1.clean <- dat.1 %>% select(c(False.Hypotheses, Obs, Correlation, SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) %>%
  rename('SH.Ln.LB' = SH.Ln.LB, 'SH.Sn.LB' = SH.Sn.LB, 'CQ' = CQ, 'SD' = SD) 

## Transforming Data 
transform.data <- dat1.clean %>% melt(id.vars = colnames(dat1.clean)[1:3], variable.name = 'Tests',  value.name = 'Power') %>% 
  mutate(Is.SH = str_detect(Tests,'SH'))

ggsave('AR_NonEqual.pdf', FacetPower(transform.data, 'Correlation'), width = 7.5, height =8.5)



# Plotting Setting 3 - Block Covariance Matrix  ---------------------------

## Block Size 20 
dat.3.size20 <- dat.3 %>% 
  filter(Block_Size == 20) %>% 
  select(c(False.Hypotheses, Obs, Rho, SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) %>%
  rename('SH.Ln.LB' = SH.Ln.LB, 'SH.Sn.LB' = SH.Sn.LB, 'CQ' = CQ, 'SD' = SD, 'Correlation' = Rho) 

transform.data <- dat.3.size20 %>% melt(id.vars = colnames(dat.3.size20)[1:3], variable.name = 'Tests',  value.name = 'Power') %>% 
  mutate(Is.SH = str_detect(Tests,'SH'))

plot.3.size.20 <- FacetPower(transform.data, col.x = 'Correlation')

ggsave('Block20_NonEqual.pdf', plot.3.size.20, width = 7.5, height =8.5)

## Block Size 50 
dat.3.size50 <- dat.3 %>% 
  filter(Block_Size == 50) %>% 
  select(c(False.Hypotheses, Obs, Rho, SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) %>%
  rename('SH.Ln.LB' = SH.Ln.LB, 'SH.Sn.LB' = SH.Sn.LB, 'CQ' = CQ, 'SD' = SD, 'Correlation' = Rho) 

transform.data <- dat.3.size50 %>% melt(id.vars = colnames(dat.3.size50)[1:3], variable.name = 'Tests',  value.name = 'Power') %>% 
  mutate(Is.SH = str_detect(Tests,'SH'))

plot.3.size.50 <- FacetPower(transform.data, col.x = 'Correlation')

ggsave('Block50_NonEqual.pdf', plot.3.size.50, width = 7.5, height =8.5)

## Block Size 100 

dat.3.size100 <- dat.3 %>% 
  filter(Block_Size == 100) %>% 
  select(c(False.Hypotheses, Obs, Rho, SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) %>%
  rename('SH.Ln.LB' = SH.Ln.LB, 'SH.Sn.LB' = SH.Sn.LB, 'CQ' = CQ, 'SD' = SD, 'Correlation' = Rho) 
transform.data <- dat.3.size100 %>% melt(id.vars = colnames(dat.3.size50)[1:3], variable.name = 'Tests',  value.name = 'Power') %>% 
  mutate(Is.SH = str_detect(Tests,'SH'))
plot.3.size.100 <- FacetPower(transform.data, col.x = 'Correlation')

ggsave('Block100_NonEqual.pdf', plot.3.size.100, width = 7.5, height =8.5)



# Plotting Setting 5 - Model Cai 7 ----------------------------------------
dat.5.clean <- dat.5 %>% select(c(False.Hypotheses, Obs, SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) %>%
  rename('SH.Ln.LB' = SH.Ln.LB, 'SH.Sn.LB' = SH.Sn.LB, 'CQ' = CQ, 'SD' = SD) 

transform.data <- dat.5.clean %>% melt(id.vars = colnames(dat.5.clean)[1:2], variable.name = 'Tests',  value.name = 'Power') %>% 
  mutate(Is.SH = str_detect(Tests,'SH'))

plot.cai.model.7 <- FacetPower(transform.data, col.x = 'False.Hypotheses', no.facet.flag = TRUE)

ggsave('CaiModel7_NonEqual.pdf', plot.cai.model.7, width = 7.5, height =8.5)


# Plotting Setting 7 - Equicorrelated Covariance --------------------------

equicorr.dat <- dat.7 %>% 
  select(c(False.Hypotheses, Obs, Correlation, SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) %>%
  rename('SH.Ln.LB' = SH.Ln.LB, 'SH.Sn.LB' = SH.Sn.LB, 'CQ' = CQ, 'SD' = SD) 

transform.data <- equicorr.dat %>% melt(id.vars = colnames(equicorr.dat)[1:3], variable.name = 'Tests',  value.name = 'Power') %>% 
  mutate(Is.SH = str_detect(Tests,'SH'))

plot.equicorr <- FacetPower(transform.data, col.x = 'Correlation')

ggsave('Equicorr_NonEqual.pdf', plot.equicorr, width = 7.5, height = 8.5)



# Ala- Princeton  ---------------------------------------------------------
power.df <- bind_rows(dat.1, dat.3, dat.5, dat.7) 

## Split By Observations 
power.df %>% select(-c(Block_Size, Correlation))
power.test.obs20 <- power.df %>% 
  filter(Obs == 20) %>% 
  select(c(SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) 

power.test.obs50 <- power.df %>% 
  filter(Obs == 50) %>% 
  select(c(SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) 

power.test.obs100 <- power.df %>% 
  filter(Obs == 100) %>% 
  select(c(SD, CQ, SH.Ln.LB, SH.Sn.LB, Simes)) 

power.test.obs20  <- apply(apply(power.test.obs20, 1, function(x) x / max(x)), 1, min)
power.test.obs50  <- apply(apply(power.test.obs50, 1, function(x) x / max(x)), 1, min)
power.test.obs100 <- apply(apply(power.test.obs100, 1, function(x) x / max(x)), 1, min)

(prince.df <- round(rbind(power.test.obs20, power.test.obs50, power.test.obs100), 3))
write.csv(prince.df, 'nonequal.ala.princeton.csv')
