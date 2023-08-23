library(UpSetR)

# read the csv file into a data frame
Day1_ADvsWT_data <- read.csv("./input/angelina_code/day 1 AD vs WT upset plot data.csv")
colnames(Day1_ADvsWT_data) <- c("gene","Tau vs WT", "AB vs WT", "Tau:AB vs WT")

# create the upset plot
pdf(file = "./output/figure4_a.pdf",width=5.38,height=5.13)
upset(Day1_ADvsWT_data,
      mainbar.y.label = "Intersections",
      sets.x.label = "Insoluble proteins",
      #sets.bar.color = c("grey60","#6699CC","#EE99AA","#EECC66"),
      text.scale = 2,
      point.size = 3,
      order.by = "freq")
dev.off()
