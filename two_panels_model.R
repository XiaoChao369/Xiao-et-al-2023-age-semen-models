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
# The two-panels model is a Bayesian regularized neural network (brnn) model,
# including 21 sperm-specific AR-CpG markers, namely 
# cg01789162, cg11262154, cg19998819, cg27231587, cg18037145, cg19983027, 
# cg27111970, cg03634854, cg03030301, cg04119405, cg25715498, cg06304190, 
# cg06979108, cg12837463, cg12277678, cg13872326, cg20602007, cg25187042,
# cg21843517, cg04123357, and cg24812634. This model can be used to estimate 
# individual age from semen DNA or better from sperm DNA.
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
two_panels_model <- readRDS("two_panels_brnn_model.rds")

## Prepare data ---------------------------------------------------------------
# Read the genotyping table or the sizing table
setwd(dirData)
geno_data <- read.table("SN001_two_panels.txt", header = TRUE)

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
panel <- c("Sperm_AR-CpG_Panel_I", "Sperm_AR-CpG_Panel_II")
markers <- c("cg01789162", "cg11262154", "cg19998819", "cg27231587", 
             "cg18037145", "cg19983027", "cg27111970", "cg03634854", 
             "cg03030301", "cg25715498", "cg04119405", "cg06304190", 
             "cg06979108", "cg12837463", "cg12277678", "cg13872326", 
             "cg20602007", "cg25187042", "cg21843517", "cg04123357", 
             "cg24812634")
geno_data <- geno_data %>%
  dplyr::filter(Panel %in% panel) %>%
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

# # A tibble: 21 ¡Á 4
# # Groups:   Sample.Name, Panel [2]
#    Sample.Name Panel                Marker     Meth.Level
#    <chr>       <chr>                <chr>           <dbl>
#  1 SN001       Sperm_AR-CpG_Panel_I cg01789162     0.123 
#  2 SN001       Sperm_AR-CpG_Panel_I cg03030301     0.380 
#  3 SN001       Sperm_AR-CpG_Panel_I cg03634854     0.299 
#  4 SN001       Sperm_AR-CpG_Panel_I cg04119405     0.531 
#  5 SN001       Sperm_AR-CpG_Panel_I cg11262154     0.890 
#  6 SN001       Sperm_AR-CpG_Panel_I cg18037145     0.701 
#  7 SN001       Sperm_AR-CpG_Panel_I cg19983027     0.0707
#  8 SN001       Sperm_AR-CpG_Panel_I cg19998819     0.838 
#  9 SN001       Sperm_AR-CpG_Panel_I cg25715498     0.897 
# 10 SN001       Sperm_AR-CpG_Panel_I cg27111970     0.143 
# # ¡­ with 11 more rows

# Save methylation data into a .csv file
setwd(dirResult)
write.csv(meth_data, "meth_data_two_panels.csv", quote = FALSE, row.names = FALSE)

# Convert long data to wide data
meth_data <- ungroup(meth_data) %>%
  dplyr::select(-Panel) %>%
  tidyr::spread(key = Marker, value = Meth.Level)

# # A tibble: 1 ¡Á 22
#   Sample.Name cg01789162 cg03030301 cg03634854 cg04119405 cg04123357  ...
#   <chr>            <dbl>      <dbl>      <dbl>      <dbl>      <dbl>  ...
# 1 SN001            0.123      0.380      0.299      0.531      0.847  ...

## Estimate age using the panel-I model ---------------------------------------
# Estimate age
estimated_age <- data.frame(
  Sample.Name = meth_data$Sample.Name,
  Estimated.Age = predict(two_panels_model, newdata = meth_data[,-1])
)
# Save estimated ages into a .csv file
setwd(dirResult)
write.csv(estimated_age, "estimated_age_two_panels.csv", quote = FALSE, 
          row.names = FALSE)
# Print estimated ages to the R console
print(estimated_age)

#   Sample.Name Estimated.Age
# 1       SN001       48.6893





