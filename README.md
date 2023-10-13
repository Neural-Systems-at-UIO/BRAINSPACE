# BRAINSPACE

This repository contains scripts and files relating to the BRAINSPACE project for the analysis of data with/from the QCAlign software and Nutil software.

More information on the [BRAINSPACE project](https://www.humanbrainproject.eu/en/collaborate-hbp/partnering-projects/brainspace/). 

More information on the [QUINT workflow](https://quint-workflow.readthedocs.io/en/latest/).

# Overview of files

The `Nutil_Validation` folder contains two test datasets that were used to validate Nutil Quantifier v0.7.0 for the parameters selected for the BRAINSPACE project. The test datasets contained objects of known size and spatial location and were run through Nutil with the parameters specified in the .NUT file. For both datasets, the Nutil output matched the ground truth for the dataset. 

The `Intermediate_hierarchy.txt` lists the custom regions used for QUINT analysis in the BRAINSPACE project. Each custom region is comprised of regions from the Allen Mouse Brain Atlas Common Coordinate Framework v3 (CCFv3). The assigned custom region name is listed in row 1, with atlas IDs assigned to this region listed below. 

The `NutilFileExample.nut` lists the Nutil parameters used for the study. 

The `SupplementaryTables` folder contains supplementary tables for the BRAINSPACE article. 

The `Scripts` folder contains R scripts relating to analysis of Nutil output and QCAlign output. 

* `ANOVA_Wilcoxon_Statistics.R` 

* `CompileNutilOutput_IntermediateRegions.R` 

* `CompileQCAlignOutput.R` 

* `MultilevelCorrelations_SampleMatch.R`



 

# How to cite


**BRAINSPACE Project** 

Gurdon B, Yates SC, Csucs G, Groeneboom NE, Hadad N, Telpoukhovskaia M, Ouellette A, Ouellette T, O’Connell K, Singh S, Murdy M, Merchant E, Bjerke I, Kleven H, Schlegel U, Puchades MA, Leergaard TB, Bjaalie JG, and Kaczorowski CC. Detecting the effect of genetic diversity on brain composition in an Alzheimer’s disease mouse model. BioRxiv Preprint. https://doi.org/10.1101/2023.02.27.530226 

**QUINT**

Yates SC, Groeneboom NE, Coello C, Lichtenthaler SF, Kuhn PH, Demuth HU,Hartlage-Rübsamen M, Roßner S, Leergaard T, Kreshuk A, Puchades MA, Bjaalie JG. QUINT: Workflow for quantification and spatial analysis of features in histological images from rodent brain. *Front Neuroinform.* 2019 Dec 3;13:75. https://doi.org/10.3389/fninf.2019.00075.

**QuickNII (RRID:SCR_016854)**
   
Puchades MA, Csucs G, Lederberger D, Leergaard TB and Bjaalie JG. Spatial registration of serial microscopic brain images to three-dimensional reference atlases with the QuickNII tool. PLosONE, 2019, 14(5): e0216796. https://doi.org/10.1371/journal.pone.0216796

**VisuAlign (RRID:SCR_017978)**

Gurdon B, Yates SC, et al. Detecting the effect of genetic diversity on brain-wide cellular and pathological changes in a novel Alzheimer’s disease mouse model. Manuscript in preparation.

**ilastik**

Berg S., Kutra D., Kroeger T., Straehle C.N., Kausler B.X., Haubold C., et al. (2019) ilastik:interactive machine learning for (bio) image analysis. Nat Methods. 16, 1226–1232. https://doi.org/10.1038/s41592-019-0582-9

**Nutil (RRID: SCR_017183)**
   
Groeneboom NE, Yates SC, Puchades MA and Bjaalie JG. Nutil: A Pre- and Post-processing Toolbox for Histological Rodent Brain Section Images. Front. Neuroinform. 2020,14:37. https://doi.org/10.3389/fninf.2020.00037

**QCAlign (RRID:SCR_023088)**

Gurdon B, Yates SC, et al. Detecting the effect of genetic diversity on brain-wide cellular and pathological changes in a novel Alzheimer’s disease mouse model. Manuscript in preparation. 

# Acknowledgements

The BRAINSPACE project received support from the EBRAINS infrastructure with funding from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under the Framework Partnership Agreement No. 650003 (HBP FPA).

# Contact us

For advice on using the tools in the QUINT workflow, contact: support@ebrains.eu



