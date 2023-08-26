# Copyright and License Information -------------------------------------------
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Contacts and bug reports ----------------------------------------------------
# 
# Chao Xiao
# xiaochao369@hust.edu.cn


# Description -----------------------------------------------------------------
# 
# The panel II model is a neural network (neuralnet) model, including 10 
# sperm-specific AR-CpG markers, namely cg06304190, cg06979108, cg12837463, 
# cg12277678, cg13872326, cg20602007, cg25187042, cg21843517, cg04123357, 
# and cg24812634. This model can be used to estimate individual age from 
# semen DNA or better from sperm DNA.
#
# Input: a genotyping table or a sizing table in '.txt' format exported
# using a GeneMapper ID/IDX Software. Ensure the exported data are
# organized as shown below and contain at least the columns named "Sample.Name",
# "Panel", "Marker", "Allele.1", "Allele.2", "Height.1", and "Height.2".
#
# Sample.Name	Panel	Marker	Allele.1	Allele.2	Height.1	Height.2	...
# SN001	Sperm_AR-CpG_Panel_II	BC	C	T	NA	1311	...
# SN001	Sperm_AR-CpG_Panel_II	cg06304190	C	T	1006	783	...
# SN001	Sperm_AR-CpG_Panel_II	cg06979108	C	T	1743	646	...
# SN001	Sperm_AR-CpG_Panel_II	cg12837463	C	T	1546	1703	...
# SN001	Sperm_AR-CpG_Panel_II	cg20828122	C	T	128	1619	...
# SN001	Sperm_AR-CpG_Panel_II	cg12277678	C	T	1003	495	...
# SN001	Sperm_AR-CpG_Panel_II	cg13872326	C	T	1178	398	...
# SN001	Sperm_AR-CpG_Panel_II	cg20602007	C	T	1237	422	...
# SN001	Sperm_AR-CpG_Panel_II	cg25187042	C	T	1606	590	...
# SN001	Sperm_AR-CpG_Panel_II	cg21843517	C	T	1459	649	...
# SN001	Sperm_AR-CpG_Panel_II	cg04123357	C	T	2452	442	...
# SN001	Sperm_AR-CpG_Panel_II	cg24812634	C	T	112	2789	...
#
# Output: the donor's age estimated by the panel II model.


# Encode statement ------------------------------------------------------------
# -*- coding: UTF-8 -*-


## Load packages --------------------------------------------------------------
if (!requireNamespace("caret", quietly = TRUE))
  install.packages("caret", dependencies = c("Depends", "Suggests"))
library(caret)  # for model calling
if (!requireNamespace("dplyr", quietly = TRUE))                     
  install.packages("dplyr")
library(dplyr)  # for data wrangling
if (!requireNamespace("tidyr", quietly = TRUE))                     
  install.packages("tidyr")
library(tidyr)  # for data wrangling


## Load panel-II model ---------------------------------------------------------
# https://github.com/XiaoChao369/Xiao-et-al-2023-age-semen-models
setwd(dirModel)
panel_II_model <- readRDS("panel_II_neuralnet_model.rds")

## Prepare data ---------------------------------------------------------------
# Read the genotyping table or the sizing table
setwd(dirData)
geno_data <- read.table("SN001_panel_II.txt", header = TRUE)

# Sample.Name	Panel	Marker	Allele.1	Allele.2	Height.1	Height.2	...
# SN001	Sperm_AR-CpG_Panel_II	BC	C	T	NA	1311	...
# SN001	Sperm_AR-CpG_Panel_II	cg06304190	C	T	1006	783	...
# SN001	Sperm_AR-CpG_Panel_II	cg06979108	C	T	1743	646	...
# SN001	Sperm_AR-CpG_Panel_II	cg12837463	C	T	1546	1703	...
# SN001	Sperm_AR-CpG_Panel_II	cg20828122	C	T	128	1619	...
# SN001	Sperm_AR-CpG_Panel_II	cg12277678	C	T	1003	495	...
# SN001	Sperm_AR-CpG_Panel_II	cg13872326	C	T	1178	398	...
# SN001	Sperm_AR-CpG_Panel_II	cg20602007	C	T	1237	422	...
# SN001	Sperm_AR-CpG_Panel_II	cg25187042	C	T	1606	590	...
# SN001	Sperm_AR-CpG_Panel_II	cg21843517	C	T	1459	649	...
# SN001	Sperm_AR-CpG_Panel_II	cg04123357	C	T	2452	442	...
# SN001	Sperm_AR-CpG_Panel_II	cg24812634	C	T	112	2789	...

# Subset data
panel <- "Sperm_AR-CpG_Panel_II"
markers <- c("cg06304190", "cg06979108", "cg12837463", "cg12277678",
             "cg13872326", "cg20602007", "cg25187042", "cg21843517", 
             "cg04123357", "cg24812634")
geno_data <- geno_data %>%
  dplyr::filter(Panel == panel) %>%
  dplyr::filter(Marker %in% markers)

# Check data
stopifnot(exprs = {
  isTRUE(sum(unique(geno_data$Marker) %in% markers) != 0)
  isTRUE(sum(markers %in% unique(geno_data$Marker)) != 0)
  isTRUE((nrow(geno_data) %% length(markers)) %in% 0)
  all.equal(unique(geno_data$Panel), panel)
})

# Calculate methylation level for each marker
geno_data$Height.1 <- as.numeric(geno_data$Height.1)
geno_data$Height.2 <- as.numeric(geno_data$Height.2)
meth_data <- geno_data %>%
  dplyr::group_by(Sample.Name, Panel, Marker) %>%
  dplyr::summarise(Meth.Level = Height.1/(Height.1 + Height.2))

# # A tibble: 10 ?? 4
# # Groups:   Sample.Name, Panel [1]
#    Sample.Name Panel                 Marker     Meth.Level
#    <chr>       <chr>                 <chr>           <dbl>
#  1 SN001       Sperm_AR-CpG_Panel_II cg04123357     0.847 
#  2 SN001       Sperm_AR-CpG_Panel_II cg06304190     0.562 
#  3 SN001       Sperm_AR-CpG_Panel_II cg06979108     0.730 
#  4 SN001       Sperm_AR-CpG_Panel_II cg12277678     0.670 
#  5 SN001       Sperm_AR-CpG_Panel_II cg12837463     0.476 
#  6 SN001       Sperm_AR-CpG_Panel_II cg13872326     0.747 
#  7 SN001       Sperm_AR-CpG_Panel_II cg20602007     0.746 
#  8 SN001       Sperm_AR-CpG_Panel_II cg21843517     0.692 
#  9 SN001       Sperm_AR-CpG_Panel_II cg24812634     0.0386
# 10 SN001       Sperm_AR-CpG_Panel_II cg25187042     0.731

# Save methylation data into a .csv file
setwd(dirResult)
write.csv(meth_data, "meth_data_panel_II.csv", quote = FALSE, row.names = FALSE)

# Convert long data to wide data
meth_data <- ungroup(meth_data) %>%
  dplyr::select(-Panel) %>%
  tidyr::spread(key = Marker, value = Meth.Level)

# # A tibble: 1 ?? 12
#   Sample.Name cg04123357 cg06304190 cg06979108 cg12277678 cg12837463  ...
#   <chr>            <dbl>      <dbl>      <dbl>      <dbl>      <dbl>  ...
# 1 SN001            0.847      0.562      0.730      0.670      0.476  ...

## Estimate age using the panel II model ---------------------------------------
# Estimate age
estimated_age <- data.frame(
  Sample.Name = meth_data$Sample.Name,
  Estimated.Age = predict(panel_II_model, newdata = meth_data[,-1])
)
# Save estimated ages into a .csv file
setwd(dirResult)
write.csv(estimated_age, "estimated_age_panel_II.csv", quote = FALSE, 
          row.names = FALSE)
# Print estimated ages to the R console
print(estimated_age)

#   Sample.Name Estimated.Age
# 1       SN001      48.98424





