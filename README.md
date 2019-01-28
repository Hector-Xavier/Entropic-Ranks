# Entropic Ranks readme

Depends:
RankProd (reliable compatibility with RankProd versions up to 2.44.0), entropy, factoextra

# Entropic Ranks (docker) readme

## entropic_ranks

Description: Performs an Entropic Ranks analysis on a data set, returning a list containing downregulated and upregulated features. May be used supervised, returning the full feature list and printing the suggested cutoff points for later manual trimming, or unsupervised, returning only the information-rich feature list. In the unsupervised mode, the lists of information-rich features may be exported as tab-delimited .txt files automatically.

## Usage:
> entropic_ranks(data_under_analysis,population_vector,data_origin=NULL,granularity=1,supervised=FALSE,process_log=FALSE,export_plots=FALSE,create_output_files=FALSE,is_logged=TRUE,logbase=2,huge_feature_list=FALSE)

## Arguments:
**data_under_analysis** - Tab-delimited .txt table with rows representing features, columns representing samples and cells containing the values to be compared. Rownames and column names must be unique. (see included test data)

**population_vector** - Tab-delimited .txt table with a single column, one row per sample and 0 or 1 as the table values. Header and rownames must be included. Denotes the two sample subpopulations to be compared. (see included test data)

**data_origin** - Tab-delimited .txt table with a single column, one row per sample. Header and rownames must be included. The data must be labels differentiating the origin of each sample. If NULL, it defaults to assuming that data are of the same origin. To be used only if the data are from different experiments or publications. (default: null)

**granularity** - The sliding window step, corresponding to the granularity of the partitioning process (feature-by-feature, or partitioning by 5-feature steps). (default: 1)

**supervised** - If TRUE, the full list of differentially behaving features is returned and the tables of suggested cutoff points are printed. If FALSE, only the list of information-rich features is returned. (default: FALSE)

**process_log** - If TRUE, statistics of the entropic_ranks execution will be printed and plots of the entropy distributions and clustering qualities will be generated. (default: FALSE)

**export_plots** - If TRUE, png plots of the entropy distributions and clustering qualities will be exported as files in a folder system created in the current working directory. (default: TRUE)

**create_output_files** - If TRUE, the feature lists of information-rich features will be automatically exported in the working directory as tab-delimited .txt files. Ignored if supervised is set to TRUE. (default: TRUE)

**is_logged** - Set to TRUE if the values are log-transformed and you want to export the Fold Change instead of the Log Fold Change in .txt files. Ignored if supervised is set to TRUE. (default: TRUE)

**logbase** - The base of the log transformation. Ignored if supervised is set to TRUE or if create_output_files is set to FALSE. (default: 2)

**huge_feature_list** - Only set to TRUE if the entropic_ranks fails to run due to huge feature lists returned by RankProd (e.g. more than 25000-30000 features) and you can be reasonably sure that less than 1000 are differentially expressed and information-rich. If TRUE, entropic_analysis will only investigate the first 20000 features and isolare information-rich features from among them. The issue may consistently appear due to RAM shortage when analyzing methylation data. (default: FALSE)


## isolate_significant_elements

Description: Calls entropic_analysis repeatedly on an ordered vector to generate a set of possible cutoff points over a range of different window sized and bins. Identifies the most consistent cutoff point. May be used as an unsupervised procedure, returning the cutoff point, or as a supervised procedure, returning the set of possible cutoff points for the researcher to further investigate.

## Usage:
> isolate_significant_elements (ordered_vector,granularity=1,supervised=FALSE,process_log=FALSE,export_plots=FALSE,path=NULL)

## Arguments:
**ordered_vector** - An ordered vector (produced by RPadvance) to be analyzed. It corresponds to the rank product score distribution produced by a two-class case experiment.

**granularity** - The sliding window step, corresponding to the granularity of the partitioning process (feature-by-feature, or partitioning by 5-feature steps).

**supervised** - If TRUE, an integer vector of possible cutoff points is returned. If FALSE, only the most consistently appearing cutoff point is returned.

**process_log** - If TRUE, statistics of the isolate_significant_elements execution will be printed and plots of the entropy distributions and clustering qualities will be generated.

**export_plots** - If TRUE, png plots of the entropy distributions and clustering qualities will be exported as files in the directory provided in the path variable.

**path** - Of type character, specifies the path for plot export. If NULL, it defaults to working directory/Entropic Ranks plots. Ignored if export_plots is set to FALSE.


## entropic_analysis
Description: The function performs an entropic analuysis of an ordered vector produced by RPadvance to identify the end point of its information-rich area. Entropy scores are calculated with the use of a sliding window.

## Usage:
> entropic_analysis (ordered_vector,step_up=1,window_size,bins,verbose=FALSE,export_plots=FALSE,path=NULL)

## Arguments:
**ordered_vector** - An ordered vector (produced by RPadvance) to be analyzed. It corresponds to the rank product score distribution produced by a two-class case experiment.

**step_up** - The sliding window step, corresponding to the granularity of the partitioning process (feature-by-feature, or partitioning by 5-feature steps).

**window_size** - The size of the sliding window used to calculate entropy scores. Must be less than 80% of the total number of features. Recommended values: 50-200.

**bins** - The number of bins used for the discretization part of entropy calculation. Recommended values: 5 to 25.

**verbose** - If TRUE, statistics of the entropic_analysis execution will be printed and plots of the entropy distribution and clustering quality will be generated.

**export_plots** - If TRUE, png plots of the entropy distribution and clustering quality will be exported as files in the directory provided in the path variable.

**path** - Of type character, specifies the path for plot export. If NULL, it defaults to working directory/Entropic Ranks plots. Ignored if export_plots is set to FALSE.
