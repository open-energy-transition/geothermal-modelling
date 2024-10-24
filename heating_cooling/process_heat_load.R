rm(list=ls())

data_dir <- "."
data_fl <- "residential_load_by_states.csv"

library(ggplot2)
library(lubridate)
library(usmap)
library(tidyverse)

load_df <- read.csv(file.path(data_dir, data_fl))
load_df["time"] = ymd_hms(unlist(load_df["timestamp"]))


time_snap <- "2018-02-01 07:15:00"
load_snap_df <- load_df %>%
    dplyr::filter(timestamp %in% time_snap)


# test on data for selected states --------------------------------------------
test_df <- load_df %>%
    dplyr::filter(`in.state` %in% c("AL", "WI"))

ts_plot <- ggplot(
    test_df,
    aes(x=time,
        y=overall_resload,
        color=`in.state`,
        by_group=`in.state`)
    )+
    geom_line(alpha=0.75)
ggsave("test_ts_al_wi.pdf", width=7, height=3.5)

# test on all the available data ----------------------------------------------
ts_plot <- ggplot(
    load_df,
    aes(x=time,
        y=overall_resload,
        color=`in.state`,
        by_group=`in.state`)
    )+
    geom_line(alpha=0.3)
ggsave("test_ts_all.pdf", width=7, height=3.5)


# check cross-correlations ----------------------------------------------------
res_wide <- load_df %>% 
    pivot_wider(
        names_from = `in.state`,
        values_from = overall_resload
    )
cor_mx <- cor(res_wide[, -1])

pdf("cor_heatmap.pdf", width=12, height=7)
    levelplot(
        cor_mx, col.regions = heat.colors,
        # scales=list(log="e", x=list(cex=.3), y=list(cex=.3)), xlab=list(cex=.05)
        scales=list(x=list(cex=0.85), y=list(cex=.75), xlab=list(cex=.5)),
        aspect="fill"
    )
dev.off()

which(cor_mx==min(cor_mx), arr.ind=TRUE)
# >   row col
# > FL   9   3
# > AZ   3   9

cor_mx[, 3]
# >        AL        AR        AZ        CA        CO        CT        DC        DE        FL        GA        IA 
# > 0.4525742 0.5453370 1.0000000 0.8688099 0.7769535 0.5426383 0.4654390 0.4942377 0.3240430 0.4834505 0.5782203 
# >        ID        IL        IN        KS        KY        LA        MA        MD        ME        MI        MN 
# > 0.7546693 0.5337831 0.5159754 0.6380850 0.4827804 0.4426038 0.5211375 0.5026393 0.5583261 0.5258189 0.5970221 
# >        MO        MS        MT        NC        ND        NE        NH        NJ        NM        NV        NY 
# > 0.5622910 0.4612566 0.6718847 0.4773994 0.6160334 0.6411915 0.5597905 0.5046856 0.8725898 0.8985851 0.5155146 
# >        OH        OK        OR        PA        RI        SC        SD        TN        TX        UT        VA 
# > 0.4952146 0.6256118 0.7077678 0.5083118 0.5063606 0.4664931 0.6370892 0.4836999 0.5628487 0.8195506 0.4957709 
# >        VT        WA        WI 
# > 0.5461889 0.6919418 0.5665383 

# test on data for selected states
demo_cor <- load_df %>%
    dplyr::filter(`in.state` %in% c("AZ", "FL", "NV"))

ts_plot <- ggplot(
    demo_cor,
    aes(x=time,
        y=overall_resload,
        color=`in.state`,
        by_group=`in.state`)
    )+
    geom_line(alpha=0.75)
ggsave("demo_cor_ts_az_fd_nv.pdf", width=7, height=3.5)



# ts_plot <- ggplot(
#     load_df,
#     aes(x=timestamp, y=overall_resload, color=in.state)
#     )+
#     geom_line()   
