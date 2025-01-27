**README File**

**Initial Setup**

Before running any scripts, start by executing the "0.1_Setup.R" script. Within this script, set the working directory (RepoFolder) to your desired working directory.

The working directory must include the data folder, which contains all the data required to run the scripts.

The figures and output folders will be automatically created by the setup script within the specified repository folder (RepoFolder). Additionally, a temporary folder named data_tmp will be generated inside the data folder to store intermediate data produced during the analysis.

**Script Workflow**

All subsequent scripts must be run in sequential order, as some depend on results from previous scripts. These results will be stored in the output or data_tmp folders. The scripts are organized according to the sections in the STAR Methods:

	1. Null models of beak morphospace volume and shape
	2. Archetypal analysis of beak and body shape
	3. Correspondence between archetypes and physical tasks
	4. Sensitivity analyses
	5. Directionality in trait evolution
	6. Mapping speciation rates and extinction risk across morphospace

**Datasets Overview**

The analysis uses the following datasets:

	"Data1_FeedingTechniques.csv" : Includes the main foraging niche (32 categories) and the percentage of use across 10 feeding techniques (n = 9,993 species).
	"Data2_BeakTaskEnrichment.csv" : Contains percentages for each of the three beak physical tasks: Crush, Engulf, and Reach (n = 9,993 species).
	"Data3_BodyTaskEnrichment.csv" : Contains percentages for each of the three body physical tasks: Fly, Walk, and Swim (n = 9,993 species).
	"Data4_Task_sensitivity/" : A folder containing six datasets for sensitivity analyses of physical task scoring. Each dataset represents upweighted task enrichment data for a specific beak or body physical task (n = 9,993 species).	
 	"Data5_ExtinctionRisk_Confounds.csv" : Includes extinction risk categories from the IUCN Red List (2022), along with species trait data such as body mass, insularity, generation length, and habitat breadth (n = 9,987 species).
	"Data6_Diversification_DR.csv" : Provides the diversification rate metric (DR) derived from the BirdTree phylogeny (n = 9,993 species). 
 	"phylo/" : A folder containing the phylogenetic tree files used in different analyses.
