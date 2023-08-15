library(tidyverse)

# set dir relative to repo
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

# read data for registered users
ga3 <- data.table::fread("data/20230808_CaeNDR_User_GA3Report.csv")

# read data for site visits ()
ga4 <- data.table::fread("data/20230808_CaeNDR_User_Report_GA4.csv") %>%
  dplyr::mutate(tot.n.users = cumsum(`New users`),
                tot.users = cumsum(Users),
                date = lubridate::mdy(date),
                month = lubridate::month(date)) %>%
  dplyr::filter(date >= '2023-04-01' & date < '2023-08-01') %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(nu.in.month = sum(`New users`),
                u.in.month = sum(Users)) %>%
  dplyr::ungroup()

ga4_summarized <- ga4 %>%
  dplyr::distinct(month, .keep_all = T) %>%
  dplyr::mutate(avg.nu.per.month = mean(nu.in.month),
                avg.u.per.month = mean(u.in.month)) %>%
  dplyr::select(-1:-3)

# plot the new user trends
p1 <- ggplot(ga4) +
  geom_line(aes(x = date, y = tot.n.users, color = "new")) +
  geom_line(aes(x = date, y = tot.users, color = "all")) +
  theme_classic() +
  labs(y = "Cumulative number of users", x = "") +
  scale_color_manual(name='User type',
                     breaks=c('new', 'all'),
                     values=c('new'='red', 'all'='black')) +
  theme(legend.position = c(0.8, 0.2))
p1

p2 <- ggplot(ga4 %>% dplyr::distinct(month, .keep_all = T)) +
  geom_bar(aes(y = nu.in.month, x = date), stat = "identity") +
  geom_label(aes(label = nu.in.month, y = nu.in.month/2, x = date), vjust = 0) +
  theme_classic() +
  geom_hline(yintercept = mean(ga4 %>% dplyr::distinct(month, .keep_all = T) %>% dplyr::pull(nu.in.month)), linetype = 2)+
  geom_label(aes(label = mean(ga4 %>% dplyr::distinct(month, .keep_all = T) %>% dplyr::pull(nu.in.month)),
                 x = lubridate::date('2023-06-01'),
                y = mean(ga4 %>% dplyr::distinct(month, .keep_all = T) %>% dplyr::pull(nu.in.month)))) +
  labs(y = "Number of new users", x = "")
p2
  
# put them together
p3 <- cowplot::plot_grid(p1, p2, ncol = 2, align = "vh", labels = c("A", "B"))

cowplot::ggsave2(p3, filename = "plots/CaeNDR_engagment.png", width = 7.5, height = 5)  
