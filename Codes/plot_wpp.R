library(wpp2010)
library(ggplot2)
library(gridExtra)
data("tfr")
par(mfrow=c(1,5))
cex.main = 2
cex.lab =3.5
cex.axis = 2
df = as.data.frame(cbind(c(1953,1958,1963,1968,1973,1978,1983,1988,1993,1998,2003,2008),
                         as.numeric(sapply(c(tfr[tfr$country == "Germany",]),function(x) x)[-c(1,2)]),
                         as.numeric(sapply(c(tfr[tfr$country == "France",]),function(x) x)[-c(1,2)]),
                         as.numeric(sapply(c(tfr[tfr$country == "Switzerland",]),function(x) x)[-c(1,2)]),
                         as.numeric(sapply(c(tfr[tfr$country == "Luxembourg",]),function(x) x)[-c(1,2)]),
                         as.numeric(sapply(c(tfr[tfr$country == "Republic of Korea",]),function(x) x)[-c(1,2)])
                         ))
names(df) = c("year","Germany","France","Switzerland","Luxembourg","Republic_of_Korea")
df$year = names(tfr)[-c(1,2)]
my_theme <- theme_bw() +
  theme(strip.background = element_rect(fill = "white"), text = element_text(face="bold", size=12),
  )
theme_set(my_theme)

axis.text = element_text(size = 14)
axis.title = element_text(size = 16)
plot.title = element_text(size = 18, face = "bold")
             

plot_germany = ggplot(df,aes(x=year,y=Germany)) + geom_point() +
  xlab("5-year time period") +
  theme(
    axis.text = axis.text,
    axis.title = axis.title,
    plot.title = plot.title,
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_france = ggplot(df,aes(x=year,y=France)) + geom_point() +
  xlab("5-year time period") +
  theme(
    axis.text = axis.text,
    axis.title = axis.title,
    plot.title = plot.title,
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_switzerland = ggplot(df,aes(x=year,y=Switzerland)) + geom_point() +
  xlab("5-year time period") +
  theme(
    axis.text = axis.text,
    axis.title = axis.title,
    plot.title = plot.title,
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_luxembourg = ggplot(df,aes(x=year,y=Luxembourg)) + geom_point() +
  xlab("5-year time period") +
  theme(
    axis.text = axis.text,
    axis.title = axis.title,
    plot.title = plot.title,
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_republic_of_korea = ggplot(df,aes(x=year,y=Republic_of_Korea)) + geom_point() +
  ylab("Republic of Korea") + xlab("5-year time period") +
  theme(
    axis.text = axis.text,
    axis.title = axis.title,
    plot.title = plot.title,
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_combined = grid.arrange(plot_france, 
                             #plot_republic_of_korea,
                             plot_germany,
                             nrow=1)
ggsave("atelier/sim_final_n195_raw_tfr.pdf",
       plot_combined, 
       device="pdf", width=12, height=6)
