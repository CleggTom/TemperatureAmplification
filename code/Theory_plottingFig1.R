library(tidyverse)
library(cowplot)
library(magick)

#read in data
R = read_csv("./data/theory/R.csv", col_types = cols(), col_names = c("Competition","Neutral","Facilitation")) %>%
    mutate(T = seq(5,25,length.out = 24)) %>%
    pivot_longer(-T,names_to = "Interaction", values_to = "Respiration")

ER = read_csv("./data/theory/ER.csv", col_types = cols(), col_names = c("E")) %>%
    mutate(Interaction = factor(c("Competition","Neutral","Facilitation"),levels = c("Competition","Neutral","Facilitation")), E = abs(E))

cols <- c(c("Competition" = "black","Neutral" = "grey","Facilitation" = "red"))

options(repr.plot.width = 5, repr.plot.height = 5)

#plot E
E_x_ticks <- c(expression(italic(bar(a) < 0)),expression(italic(bar(a) == 0)),expression(italic(bar(a) > 0)))

E_plot <- ggplot(ER,aes(x = Interaction, y = E, color = Interaction))+
    geom_point(size = 3) +
    ylab(expression(paste(italic(E[eco])," (Ev)"))) +
    scale_color_manual(values = cols) +
    scale_x_discrete(breaks = c("Competition","Neutral","Facilitation"), labels = E_x_ticks) +
    theme_cowplot()+
    theme(legend.position = "none", text = element_text(size = 20))+
    ylim(0.75, 2.0)

#plot R
R_plot <- ggplot(R,aes(x=T, y = log(Respiration), group = Interaction, color = Interaction)) +
    geom_line()+
    scale_color_manual(values = cols) +
    theme_cowplot()+
    xlab("Temperature (°C)")+
    ylab(expression(paste("Ecosystem Respiraion ", bgroup("(",log(italic(R[eco])),")") )))+
    theme(text = element_text(size = 20),legend.position = c(0.6,0.3),legend.title = element_blank(), plot.margin = margin(0,0,0,0))+
    annotation_custom(ggplotGrob(E_plot),xmin = 5, xmax= 18, ymin = 2.5, ymax = 7)

#Plot temperature
Temp_diagram <- image_read_svg("./plots/Fig_1b.svg")
Temp_diagram <- ggdraw() + draw_image(Temp_diagram, interpolate = ) + theme(plot.margin = margin(0,0,0,0))
Temp_diagram

p_final <- plot_grid(Temp_diagram,R_plot, nrow = 2, labels = c("a","b"),label_size = 25, label_y = c(0.87,1.15), rel_heights = c(1,0.8))

ggsave("./plots/Fig1_theory.pdf",width = 5 , height = 10)