# fitting 

options(ggplot2.discrete.color = gen_colors(zurich, 12), ggplot2.discrete.fill = gen_colors(zurich, 12))

load('../../fit/skewlogis/separate.RData')

est_comb <- rbind(
    mn_backreport$est %>% mutate(sex='male', direction = 'other'),
    mn_selfreport$est %>% mutate(sex='male', direction = 'self'),
    wm_backreport$est %>% mutate(sex='female', direction = 'other'),
    wm_selfreport$est %>% mutate(sex='female', direction = 'self')
)

# parameters comparision
# 
save_to <- 'fig/param_compare_direction'
dir.create(save_to)
cd(save_to)

for (cc in unique(est_comb$ISO_A3)) {
    cat('\r', cc)
    plot <- est_comb %>% 
        rename(Scale=a_vec, Shape=b_vec, Skewness=g_vec) %>%
        filter(ISO_A3==cc) %>% 
        pivot_longer(-char(ISO_A3, sex, age, direction)) %>% 
        ggplot() + 
            geom_rect(aes(xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf, fill=sex), data=tibble(sex=char(male, female)), alpha=.1) +
            geom_point(aes(age, value, color=direction), alpha=.3) +
            facet_wrap(char(sex, name), scales='free') + 
            scale_shape_manual(values=c(1, 3)) +
            scale_fill_manual(values=c('transparent', 'grey70')) +
            guides(fill=FALSE) +
            theme(legend.position=c(.85, 1.17), legend.direction='horizontal', legend.title=element_text(size=8)) +
            labs(x='Age', color='Source of report', y='Coef.', title=cc)
    quartz_off(plot, cc, open=0)
}

cd('-')

# fitted density
save_to <- 'fig/fitted_density_by_country'
dir.create(save_to)
cd(save_to)

for (cc in unique(est$ISO_A3)) {
    cat('\r', cc)
    di <- est_comb %>% filter(ISO_A3==cc)
    plot_ages <- c(18, 20, 25, 30, 40, 50)
    pages <- seq(10, 60, .5)
    lapply(1:nrow(di), function(x) f_gllogisI(pages, di$a_vec[x], di$b_vec[x], di$g_vec[x])) %>% 
        do.call('rbind', .) %>%
        set_colnames(pages) %>% 
        cbind(dplyr::select(di, age, sex, direction)) %>%
        pivot_longer(-char(age, sex, direction), names_to='partner', values_to='dens') %>% 
        mutate(partner = as.double(partner)) %>% 
        filter(age %in% plot_ages) -> di
    
    fbi %>% 
        filter(age %in% plot_ages & ISO_A3 == cc) %>%   
        mutate(age = paste0('', age)) %>% 
        ggplot() + 
            geom_bar(aes(partner, y=..density.., fill=direction), binwidth=1, position='identity', stat='bin', alpha=.3, linetype='dashed') +
            geom_line(aes(partner, y=dens, color=direction), di, size=.7) +
            facet_grid(rows=vars(sex), cols=vars(age)) +
            coord_cartesian(ylim=c(0, .3)) +
            theme(legend.position=c(.85, 1.1), legend.title=element_text(size=7), legend.direction='horizontal') +
            scale_fill_manual(values=char(grey40, orange)) +
            scale_color_manual(values=char(grey10, darkorange)) +
            labs(title=cc, subtitle='Data and fitted', fill='Source of report', color='Source of report', x='Partner age') -> plot
    quartz_off(plot, cc, open=0)
}
pdftools::pdf_combine(list.files()[-2], 'all.pdf')
open_file('all.pdf')

cd('-')

