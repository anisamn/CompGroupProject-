

san_data <- read.csv("cellswithMYBPHL_SAN.csv", row.names = 1, header = TRUE)
lpf_data <- read.csv("cellswithMYBPHL_LPF.csv", row.names = 1, header = TRUE)
rpf_data <- read.csv("cellswithMYBPHL_RPF.csv", row.names = 1, header = TRUE)
avn_data <- read.csv("cellswithMYBPHL_AVN.csv", row.names = 1, header = TRUE)


san_mybphl <- san_data["Mybphl",]
san_mybphl = t(san_mybphl)
hist(san_mybphl)

lpf_mybphl <- lpf_data["Mybphl",]
lpf_mybphl = t(lpf_mybphl)
hist(lpf_mybphl)

rpf_mybphl <- rpf_data["Mybphl",]
rpf_mybphl = t(rpf_mybphl)
hist(rpf_mybphl)

avn_mybphl <- avn_data["Mybphl",]
avn_mybphl = t(avn_mybphl)
hist(avn_mybphl)
