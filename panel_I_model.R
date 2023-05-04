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
# The panel-I model is a support vector machine with polynomial kernel (svmPoly)
# model, including 11 AR-CpG markers, namely cg01789162, cg11262154, cg19998819,
# cg27231587, cg18037145, cg19983027, cg27111970, cg03634854, cg03030301,
# cg04119405, and cg25715498. This model can be used to estimate individual age
# from semen DNA or better from sperm DNA.
#
# Input: a genotyping table or a sizing table in '.txt' format exported
# using a GeneMapper ID/IDX Software. Ensure the exported data are
# organized as shown below and contain at least the columns named "Sample.Name",
# "Panel", "Marker", "Allele.1", "Allele.2", "Height.1", and "Height.2".
#
# Sample.Name	Panel	Marker	Allele.1	Allele.2	Height.1	Height.2  ...
# SN001	Sperm_AR-CpG_Panel_I	cg01789162	G	A	213	1524  ...
# SN001	Sperm_AR-CpG_Panel_I	cg11262154	G	A	2291	284 ...
# SN001	Sperm_AR-CpG_Panel_I	cg19998819	G	A	1406	271 ...
# SN001	Sperm_AR-CpG_Panel_I	cg27231587	G	A	541	1715  ...
# SN001	Sperm_AR-CpG_Panel_I	cg18037145	G	A	1765	754 ...
# SN001	Sperm_AR-CpG_Panel_I	cg19983027	G	A	131	1721  ...
# SN001	Sperm_AR-CpG_Panel_I	cg27111970	G	A	396	2367  ...
# SN001	Sperm_AR-CpG_Panel_I	cg03634854	G	A	801	1874  ...
# SN001	Sperm_AR-CpG_Panel_I	cg03030301	G	A	701	1146  ...
# SN001	Sperm_AR-CpG_Panel_I	cg04119405	G	A	1670	1476  ...
# SN001	Sperm_AR-CpG_Panel_I	cg25715498	G	A	3351	383 ...
#
# Output: the donor's age estimated by the panel-I model.


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


## Load panel-I model ---------------------------------------------------------
# https://github.com/XiaoChao369/Xiao-et-al-2023-age-semen-models
setwd(dirModel)
panel_I_model <- readRDS("panel_I_svmPoly_model.rds")

## Prepare data ---------------------------------------------------------------
# Read the genotyping table or the sizing table
setwd(dirData)
geno_data <- read.table("SN001_panel_I.txt", header = TRUE)

# Sample.Name	Panel	Marker	Allele.1	Allele.2	Height.1	Height.2  ...
# SN001	Sperm_AR-CpG_Panel_I	cg01789162	G	A	213	1524  ...
# SN001	Sperm_AR-CpG_Panel_I	cg11262154	G	A	2291	284 ...
# SN001	Sperm_AR-CpG_Panel_I	cg19998819	G	A	1406	271 ...
# SN001	Sperm_AR-CpG_Panel_I	cg27231587	G	A	541	1715  ...
# SN001	Sperm_AR-CpG_Panel_I	cg18037145	G	A	1765	754 ...
# SN001	Sperm_AR-CpG_Panel_I	cg19983027	G	A	131	1721  ...
# SN001	Sperm_AR-CpG_Panel_I	cg27111970	G	A	396	2367  ...
# SN001	Sperm_AR-CpG_Panel_I	cg03634854	G	A	801	1874  ...
# SN001	Sperm_AR-CpG_Panel_I	cg03030301	G	A	701	1146  ...
# SN001	Sperm_AR-CpG_Panel_I	cg04119405	G	A	1670	1476  ...
# SN001	Sperm_AR-CpG_Panel_I	cg25715498	G	A	3351	383 ...

# Subset data
panel <- "Sperm_AR-CpG_Panel_I"
markers <- c("cg01789162", "cg11262154", "cg19998819", "cg27231587", 
             "cg18037145", "cg19983027", "cg27111970", "cg03634854", 
             "cg03030301", "cg25715498", "cg04119405")
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

# # A tibble: 11 ¡Á 4
# # Groups:   Sample.Name, Panel [1]
#    Sample.Name Panel                Marker     Meth.Level
#    <chr>       <chr>                <chr>           <dbl>
#  1 SN001       Sperm_AR-CpG_Panle_I cg01789162     0.123 
#  2 SN001       Sperm_AR-CpG_Panle_I cg03030301     0.380 
#  3 SN001       Sperm_AR-CpG_Panle_I cg03634854     0.299 
#  4 SN001       Sperm_AR-CpG_Panle_I cg04119405     0.531 
#  5 SN001       Sperm_AR-CpG_Panle_I cg11262154     0.890 
#  6 SN001       Sperm_AR-CpG_Panle_I cg18037145     0.701 
#  7 SN001       Sperm_AR-CpG_Panle_I cg19983027     0.0707
#  8 SN001       Sperm_AR-CpG_Panle_I cg19998819     0.838 
#  9 SN001       Sperm_AR-CpG_Panle_I cg25715498     0.897 
# 10 SN001       Sperm_AR-CpG_Panle_I cg27111970     0.143 
# 11 SN001       Sperm_AR-CpG_Panle_I cg27231587     0.240

# Save methylation data into a .csv file
setwd(dirResult)
write.csv(meth_data, "meth_data_panel_I.csv", quote = FALSE, row.names = FALSE)

# Convert long data to wide data
meth_data <- ungroup(meth_data) %>%
  dplyr::select(-Panel) %>%
  tidyr::spread(key = Marker, value = Meth.Level)

# # A tibble: 1 ¡Á 12
#   Sample.Name cg01789162 cg03030301 cg03634854 cg04119405 cg11262154  ...
#   <chr>            <dbl>      <dbl>      <dbl>      <dbl>      <dbl>  ...
# 1 SN001            0.123      0.380      0.299      0.531      0.890  ...

## Estimate age using the panel-I model ---------------------------------------
# Estimate age
estimated_age <- data.frame(
  Sample.Name = meth_data$Sample.Name,
  Estimated.Age = predict(panel_I_model, newdata = meth_data[,-1])
)
# Save estimated ages into a .csv file
setwd(dirResult)
write.csv(estimated_age, "estimated_age_panel_I.csv", quote = FALSE, 
          row.names = FALSE)
# Print estimated ages to the R console
print(estimated_age)

#   Sample.Name Estimated.Age
# 1       SN001      47.62784





