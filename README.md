# Xiao-et-al-2023-age-semen-models

Three models for age estimation from semen using sperm-specific age-related
CpG markers

Accurate age estimation from semen has the potential to significantly narrow
down the pool of unidentified suspects in sexual assault investigations. In a
study, we trained three models using methylation data of 21 sperm-specific 
age-related CpG (AR-CpG) markers in 253 sperm DNA samples (22 years to 67.19 years).
Given that these methylation data are generated using two methylation SNaPshot
assays, designated as panel I and panel II, we named these three models as
panel-I model, panel-II model, and two-panels model, respectively.

The following steps are outlined to demonstrate the functionality of these
three models. These models are to be used with methylation data generated
using methylation SNaPshot assays.

## Download files
Before starting, download files required to run the three models:
```
# https://github.com/XiaoChao369/Xiao-et-al-2023-age-semen-models
```
## Set environment
```
options(stringsAsFactors = FALSE)
rm(list = ls())
```
## Set directory
For convenience, set the directory containing the data or saving the results:
```
dirModel <- "J:/Xiao-et-al-2023-age-semen-models"  #  model files
dirData <- paste(dirModel, "/example_data", sep = "")  #  genotyping tables
dirResult <- paste(dirModel, "/example_results", sep = "")  # for saving results
```

## panel-I model
The panel-I model is a support vector machine with polynomial kernel (svmPoly)
model, including 11 AR-CpG markers, namely cg01789162, cg11262154, cg19998819,
cg27231587, cg18037145, cg19983027, cg27111970, cg03634854, cg03030301,
cg04119405, and cg25715498. This model can be used to estimate individual age
from semen DNA or better from sperm DNA.

Prepare data:
```
# There are several ways to prepare data used for age estimation. One way is
# to use a GeneMapper ID/IDX Software to export the genotyping table or 
# the sizing table as a '.txt' file. Another way is to maker a '.txt' file
# organized as shown below and contain at least the columns named 
# "Sample.Name", "Panel", "Marker", "Allele.1", "Allele.2", "Height.1", 
# and "Height.2".

# Sample.Name	Panel	Marker	Allele.1	Allele.2	Height.1	Height.2  ...
# SN001	Sperm_AR-CpG_Panle_I	cg01789162	G	A	213	1524  ...
# SN001	Sperm_AR-CpG_Panle_I	cg11262154	G	A	2291	284 ...
# SN001	Sperm_AR-CpG_Panle_I	cg19998819	G	A	1406	271 ...
# SN001	Sperm_AR-CpG_Panle_I	cg27231587	G	A	541	1715  ...
# SN001	Sperm_AR-CpG_Panle_I	cg18037145	G	A	1765	754 ...
# SN001	Sperm_AR-CpG_Panle_I	cg19983027	G	A	131	1721  ...
# SN001	Sperm_AR-CpG_Panle_I	cg27111970	G	A	396	2367  ...
# SN001	Sperm_AR-CpG_Panle_I	cg03634854	G	A	801	1874  ...
# SN001	Sperm_AR-CpG_Panle_I	cg03030301	G	A	701	1146  ...
# SN001	Sperm_AR-CpG_Panle_I	cg04119405	G	A	1670	1476  ...
# SN001	Sperm_AR-CpG_Panle_I	cg25715498	G	A	3351	383 ...
```
Run the panel-I model on your data using the "panel_I_model.R" script provided.
The only thing you have to do is to copy the prepared data file into
the 'dirData' folder.
```
setwd(dirModel)
source('panel_I_model.R')
```
The output 'estimated_age' is in years and is automatically printed in the R
console but is also saved as a '.csv' file, which can be found in the 
'dirResult' folder.
```
print(estimated_age)

# Sample.Name Estimated.Age
# 1       SN001      47.62784
```

## panel-II model
The panel-II model is a support vector machine with radial basis function
kernel (svmRadial) model, including 10 sperm-specific AR-CpG markers, namely
cg06304190, cg06979108, cg12837463, cg12277678, cg13872326,
cg20602007, cg25187042, cg21843517, cg04123357, and cg24812634. This model can
be used to estimate individual age from semen DNA or better from sperm DNA.

Prepare data:
```
# There are several ways to prepare data used for age estimation. One way is
# to use a GeneMapper ID/IDX Software to export the genotyping table or 
# the sizing table as a '.txt' file. Another way is to maker a '.txt' file
# organized as shown below and contain at least the columns named 
# "Sample.Name", "Panel", "Marker", "Allele.1", "Allele.2", "Height.1", 
# and "Height.2".

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
```

Run the panel-II model on your data using the "panel_II_model.R" script
provided. The only thing you have to do is to copy the prepared data file into
the 'dirData' folder.
```
setwd(dirModel)
source('panel_II_model.R')
```
The output 'estimated_age' is in years and is automatically printed in the R
console but is also saved as a '.csv' file, which can be found in the 
'dirResult' folder.
```
print(estimated_age)

#   Sample.Name Estimated.Age
# 1       SN001      48.98424
```

## two-panels model
The two-panels model is a Bayesian regularized neural network (brnn) model,
including 21 sperm-specific AR-CpG markers, namely 
cg01789162, cg11262154, cg19998819, cg27231587, cg18037145, cg19983027, 
cg27111970, cg03634854, cg03030301, cg04119405, cg25715498, cg06304190, 
cg06979108, cg12837463, cg12277678, cg13872326, cg20602007, cg25187042,
cg21843517, cg04123357, and cg24812634. This model can be used to estimate 
individual age from semen DNA or better from sperm DNA.

Prepare data:
```
# There are several ways to prepare data used for age estimation. One way is
# to use a GeneMapper ID/IDX Software to export the genotyping table or 
# the sizing table as a '.txt' file. Another way is to maker a '.txt' file
# organized as shown below and contain at least the columns named 
# "Sample.Name", "Panel", "Marker", "Allele.1", "Allele.2", "Height.1", 
# and "Height.2".

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
```

Run the two-panels model on your data using the "two-panels_model.R" script 
provided. The only thing you have to do is to copy the prepared data file into
the 'dirData' folder.
```
setwd(dirModel)
source('two_panels_model.R')
```
The output 'estimated_age' is in years and is automatically printed in the R
console but is also saved as a '.csv' file, which can be found in the
'dirResult' folder.
```
print(estimated_age)

#   Sample.Name Estimated.Age
# 1       SN001       48.6893
```


Contacts and bug reports
========================

Chao Xiao
xiaochao369@hust.edu.cn


Copyright and License Information
=================================

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
