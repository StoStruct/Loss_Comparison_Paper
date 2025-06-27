Two MATLAB scripts that create box plots comparing Expected Annual Loss (EAL) data for steel buildings.
Scripts

Figure_4_NormalizedEAL.m
Plots: Box plots of EAL by number of stories (1, 5, 9, 14, 19) comparing Assembly, Component, and SLF methodologies.

Data Files:

Building_Info_URGENT.xlsx - Building parameters
HAZUS_EAL_NormEAL_MAXEDP.xlsx - HAZUS Assembly Loss Data
FEMAP58_EAL_NormEAL_Latest.xlsx - FEMA P-58 Loss Data
SLF_EAL_NormEAL.xlsx - SLF Loss Data

Columns: Currently uses Column 11 (Total losses). Change to Column 7 (Structural), 8 (IDR-NS), or 9 (PFA-NS) for other loss types.


Figure_8_NormalizedEAL_Comparison_EDP_Statistic.m
Plots: Box plots comparing 6 different EDP statistics: HAZUS Max, HAZUS Mean, HAZUS Weighted Averages, Component, and SLF.

Data Files:

HAZUS_EAL_NormEAL_MAXEDP.xlsx - HAZUS Max
HAZUS_NormEAL_Mean.xlsx - HAZUS Mean
HAZUS_NormEAL_Weighted_Testing4.xlsx - HAZUS WA Max-Mean
HAZUS_NormEAL_Weighted_Testing_25-75Percentile.xlsx - HAZUS WA 25-75
FEMAP58_EAL_NormEAL_Latest.xlsx - Component
SLF_EAL_NormEAL.xlsx - SLF

Columns: Uses Column 11 (Total losses) only.