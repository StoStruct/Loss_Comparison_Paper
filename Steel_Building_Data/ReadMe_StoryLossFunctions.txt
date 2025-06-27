# Story Loss Functions (SLFs) 

## What are Story Loss Functions (SLFs)?

**Story Loss Functions (SLFs)** represent the aggregated loss of all damageable components present in a specific story with respect to Engineering Demand Parameter (EDP) levels. SLFs provide a computationally efficient alternative to detailed component-based FEMA P-58 assessments while maintaining reasonable accuracy.

Instead of analyzing hundreds of individual components, SLFs create simplified relationships between:
- **Input:** Inter-story drift ratio (IDR) or peak floor acceleration (PFA)  
- **Output:** Expected monetary loss per story (USD)

**Author:** Shiva Baddipalli  
**Institution:** Utah State University, USA  
**Email:** shivalinga.baddipalli@usu.edu

## Reference
Shahnazaryan, D., O'Reilly, G.J., and Monteiro, R. (2021). "Story loss functions for seismic design and assessment: Development of tools and application." *Earthquake Spectra*, 37(4), 2813-2839.

## How SLFs are Computed in This Study

### Step 1: Component Inventory Development
Python codes require **inventory data** containing:
- **Component types** (A1-A6, B1-B4, C1-C17)
- **Quantities** per story based on building geometry
- **Fragility functions** (damage probabilities)
- **Cost functions** (repair/replacement costs)

### Step 2: Building Geometry-Based Cases
SLFs are generated for different building configurations:

#### **Structural Components (9 Cases)**
Based on building geometry combinations:
- **Number of LFRS Bays:** 1, 3, 5
- **Bay Width:** 20 ft, 30 ft, 40 ft
- **Cases:** My_Inventory_CASE1.csv through My_Inventory_CASE9.csv

#### **Non-Structural Drift-Sensitive (3 Cases)**
Based on bay width variations:
- **Bay Width:** 20 ft, 30 ft, 40 ft
- **Cases:** My_Inventory_CASE1.csv through My_Inventory_CASE3.csv

#### **Non-Structural Acceleration-Sensitive**
- **Floor Level Components:** 3 cases (bay width variations)
- **Building Level Components:** 15 cases (stories × bay width combinations)

### Step 3: SLF Generation Process
For each case, Python codes (`storeyloss` library):

1. **Load inventory data** from CSV files
2. **Group components** by EDP sensitivity:
   - IDR-Structural: Column base plates, splices, main components
   - IDR-Non-Structural: Drift-sensitive components
   - PFA-Non-Structural: Acceleration-sensitive components
3. **Run Monte Carlo simulations** (20,000 realizations)
4. **Generate loss functions** relating EDP to expected loss
5. **Export results** as JSON and combined CSV files

## Folder Structure and Files

```
StoryLossFunctions_20000/
├── MyInventory_Structural/              # Structural component SLFs
│   ├── MyInventory_Structural_4comp/          # Main structural components
│   │   ├── My_Inventory_CASE1.csv             # INPUT: Inventory for case 1
│   │   ├── My_Inventory_CASE2.csv             # INPUT: Inventory for case 2
│   │   ├── ...                                # INPUT: Cases 3-9
│   │   ├── output_1.json                      # Raw SLF data case 1
│   │   ├── output_2.json                      # Raw SLF data case 2
│   │   ├── ...                                # Raw data cases 3-9
│   │   └── CombinedOutput_SLF_Struct_4comps.csv  # OUTPUT: Final SLFs
│   ├── MyInventory_Structural_Columnbase/     # Column base plates
│   └── MyInventory_Structural_Splices/        # Column splices
├── MyInventory_NSD/                     # Non-structural drift-sensitive
└── MyInventory_NSA/                     # Non-structural acceleration-sensitive
    ├── MyInventory_NSA_FloorLevelComponents/
    └── MyInventory_NSA_BuildingLevelComponents/
```

## Input Files (Inventory Data)
**Format:** `My_Inventory_CASE{X}.csv`  
**Content:** Component inventory (Type, fragility and cost functions) for specific building configuration

| Column | Description |
|--------|-------------|
| ITEM_ID | Component identifier (A1, A2, etc.) |
| EDP | Engineering demand parameter type |
| Component | Component name/description |
| Group | Component group number |
| Quantity | Number of components per story |
| Damage_States | Number of damage states |
| DS1-DS5 | Median demand for each damage state |
| Repair costs | Cost parameters for each damage state |

## Output Files (Generated SLFs)

### Individual JSON Files
- **Format:** `output_{case}.json`
- **Content:** Raw SLF data with statistical measures for each case

### Combined CSV Files  
- **Format:** `CombinedOutput_SLF_*.csv`
- **Content:** Consolidated SLFs for all cases

| Column | Description |
|--------|-------------|
| EDP Range | Engineering demand parameter values |
| SLF Case 1 | Story loss function for building configuration 1 |
| SLF Case 2 | Story loss function for building configuration 2 |
| ... | Additional cases |

## Component Categories

### Structural Components (IDR-Sensitive)
- **A1:** Column base plates
- **A2:** Column splices  
- **A3:** Shear tab connections
- **A4:** RBS one-sided connections
- **A5:** RBS two-sided connections
- **A6:** Slab components

### Non-Structural Components
- **B1-B4:** Drift-sensitive (curtain walls, partitions)
- **C1-C12:** Floor-level acceleration-sensitive
- **C13-C17:** Building-level acceleration-sensitive

## Usage in Seismic Loss Assessment

These SLFs are used in the main MATLAB loss assessment code:

1. **SLF Selection:** Based on building geometry (number of bays, bay width)
2. **Demand Prediction:** Calculate EDP levels from ground motion analysis
3. **Loss Interpolation:** Use appropriate SLF to estimate story loss
4. **Total Loss:** Sum losses across all stories and component types

## Analysis Parameters

- **Monte Carlo Realizations:** 20,000 per case
- **EDP Ranges:**
  - IDR: 0% to 300% (0.001 increment) 
  - PFA: 0g to 500g (0.05g increment)
- **Tool:** `storeyloss` Python package

For more information, check https://github.com/davitshahnazaryan3/SLFGenerator---


MATLAB Visualization Codes
Each folder contains MATLAB plotting codes that generate SLF visualizations:
Available MATLAB Files

PLOT_SLFs.m - Main plotting script for SLF visualization
PLOT_SLFs_fourcomp.m - Specialized plotting for four-component analysis
PLOT_SLFs_splices.m - Plotting script for splice components
PLOT_SLFs_NSD.m - Plotting for non-structural drift-sensitive components
PLOT_SLFs_BldgLvl.m - Plotting for building-level NSA components

How to Use MATLAB Plotting Codes

Navigate to the appropriate folder (e.g., MyInventory_Structural/)
Load the combined output CSV file (e.g., CombinedOutput_SLF_Struct_4comps.csv)
Run the corresponding MATLAB script to generate plots
View SLF curves for all cases on a single graph

**Contact:** shivalinga.baddipalli@usu.edu