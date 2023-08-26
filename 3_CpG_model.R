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
# The 3-CpG model is a support vector machine with polynomial kernel (svmPoly)
# model, including 3 AR-CpG markers, namely cg21843517, cg12837463, and cg19998819. 
# This model can be used to estimate individual age from semen DNA 
# or better from sperm DNA.
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
# Output: the donor's age estimated by the 3-CpG model.


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


## Load 3-CpG model ---------------------------------------------------------
# https://github.com/XiaoChao369/Xiao-et-al-2023-age-semen-models
setwd(dirModel)
panel_3_CpG_model <- readRDS("3_CpG_svmPoly_model.rds")

## Prepare data ---------------------------------------------------------------
# Read the genotyping table or the sizing table
setwd(dirData)
geno_data <- read.table("SN001_panels_I_&_II.txt", header = TRUE)

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
markers <- c("cg21843517", "cg12837463", "cg19998819")
geno_data <- geno_data %>%
  dplyr::filter(Marker %in% markers)

# Check data
stopifnot(exprs = {
  isTRUE(sum(unique(geno_data$Marker) %in% markers) != 0)
  isTRUE(sum(markers %in% unique(geno_data$Marker)) != 0)
  isTRUE((nrow(geno_data) %% length(markers)) %in% 0)
})

# Calculate methylation level for each marker
geno_data$Height.1 <- as.numeric(geno_data$Height.1)
geno_data$Height.2 <- as.numeric(geno_data$Height.2)
meth_data <- geno_data %>%
  dplyr::group_by(Sample.Name, Panel, Marker) %>%
  dplyr::summarise(Meth.Level = Height.1/(Height.1 + Height.2))
meth_data$Marker <- factor(meth_data$Marker, levels = markers)
meth_data <- meth_data %>%
  dplyr::arrange(Marker)

# # A tibble: 3 × 4
# # Groups:   Sample.Name, Panel [2]
# Sample.Name Panel                 Marker     Meth.Level
# <chr>       <chr>                 <fct>           <dbl>
#   1 SN001       Sperm_AR-CpG_Panel_II cg21843517      0.692
# 2 SN001       Sperm_AR-CpG_Panel_II cg12837463      0.476
# 3 SN001       Sperm_AR-CpG_Panel_I  cg19998819      0.838

# Save methylation data into a .csv file
setwd(dirResult)
write.csv(meth_data, "meth_data_panel_3_CpG.csv", quote = FALSE, row.names = FALSE)

# Convert long data to wide data
meth_data <- ungroup(meth_data) %>%
  dplyr::select(-Panel) %>%
  tidyr::spread(key = Marker, value = Meth.Level)

# # A tibble: 1 × 4
# Sample.Name cg21843517 cg12837463 cg19998819
# <chr>            <dbl>      <dbl>      <dbl>
#   1 SN001            0.692      0.476      0.838

## Estimate age using the 3-CpG model ---------------------------------------
# Estimate age
estimated_age <- data.frame(
  Sample.Name = meth_data$Sample.Name,
  Estimated.Age = predict(panel_3_CpG_model, newdata = meth_data[,-1])
)
# Save estimated ages into a .csv file
setwd(dirResult)
write.csv(estimated_age, "estimated_age_panel_3_CpG.csv", quote = FALSE, 
          row.names = FALSE)
# Print estimated ages to the R console
print(estimated_age)

# Sample.Name Estimated.Age
# 1       SN001      47.86754





