# DREAMIT
https://doi.org/10.5281/zenodo.13175583

Here we present DREAMIT (Dynamic Regulation of Expression Across Modules in Inferred Trajectories). DREAMIT aims to analyze dynamic regulatory patterns along trajectory branches, implicating transcription factors (TFs) involved in cell state transitions within scRNAseq datasets. DREAMIT uses pseudotime ordering within a robust subrange of a trajectory branch to group individual cells into bins. It aggregates the cell-based expression data into a set of robust pseudobulk measurements containing gene expression averaged within bins of neighboring cells. It then smooths trends after searching for an optimal fitting spline across the bins. DREAMIT rejects further analyzing branches that produce highly variable smoothing estimates (covariation in spline fitting parameters found to be greater than 1.0 across 80% subsampling) as these branches may represent sparse or noisy parts of the data that could produce unreliable TF inferences. 

Using the transformed smoothed data, it calculates the association between a TF and all of its predicted targets according to the TRRUST database assessed using multiple metrics (e.g., Pearson correlation, Mutual Information, Dynamic Time Warping distance, etc). Finally, DREAMIT uses a Relational Set Enrichment Analysis (RSEA) test to evaluate the significance of the TF-to-target associations and identify a core set of targets compared to a background model, which consists of arbitrarily selected targets.

![Figure 1](https://github.com/nathanmaulding/DREAMIT/blob/main/Images/Figure1.png)
Associating transcription factors (TFs) to trajectory branches via identification of TF-to-target coexpression along pseudotime. A. Expression across the cell transition in trajectory branches is used by DREAMIT to infer a dynamic view of TF-to-target gene regulation. B. Expression levels of a hypothetical gene in individual cells (y-axis) illustrating the division into arbitrary “start” and “end” states along the pseudotime of a theoretical differentiation process from stem cells to differentiated erythrocytes (x-axis). C. Differences in the expression level between “start” and “end” states may not exist which may cause Differential Expression (“DE”) approaches to miss other patterns in the data (e.g. concordant fluctuations in the middle of pseudotime). D. DREAMIT models the entirety of the expression on the branch and assesses TF-to-target relationships that look for a consistent relationship between the expression levels (y-axis) of a TF (green line) and its target genes (blue lines) along pseudotime (x-axis). E. Alignment plot showing one TF (y-axis) aligned to a “typical” target from a target set (x-axis) illustrating how allowing for a lag or delay (red line) can help a metric pick up on an association between a TF and its targets over a subset of pseudotime (blue line) in which all the targets have the same lag in expression relative to the TF.

<p align="center">
    <img src="https://github.com/nathanmaulding/DREAMIT/blob/main/Images/Figure4.png" alt="Example Image" width="65%">
</p>
Detecting TF-to-target relations using pseudotime focusing, spline smoothing, target focusing, and significance assessment via random target selection. A. Cells with outlier pseudotime assignments (red dots) compared to the other cells (green dots) are shown. B. Pseudotime focusing removes outliers from the analysis and retains cells within 1.5 times the interquartile range of all pseudotime values of the branch (see Methods). C. Cells are grouped into bins containing at least 10 cells per bin. The average expression of a gene is calculated from all the cells in a bin (blue line) and this bin-averaged expression is used for all subsequent analysis. D. Spline smoothing incorporates information from cells in neighboring bins to further smooth out the expression changes in pseudotime (green curve). E. DREAMIT quantifies TF-target relationships through pairwise tests of the spline smoothed expression of the TF (green line) and its target genes (blue lines). F. Illustration of a “rolling” metric incorporating pseudotime lag. A significant lagged correlation will be detected when several targets share the same delay. G. Target focusing employs Relational Set Enrichment Analysis (RSEA, see Methods) to identify a “core” set of targets with high association to the target (blue lines) while excluding the targets with weak or poor association (red curve). H. 75% of the targets with the highest concordance to the factor are retained. I. Significance is assessed by comparing TF-to-target metric scores to a random background in which random targets are chosen to be of the same size as the TFs original regulon (yellow lists) with a Kolmogorov-Smirnov test. J. The core set of targets (blue nodes) found by RSEA are used in the statistical analysis.


# Tutorial
For a toy example of running DREAMIT follow these steps:

1. clone the repository
```bash
git clone https://github.com/nathanmaulding/DREAMIT.git
cd DREAMIT
```

2. Unzip the data
```bash
gunzip Data/paul_toy_data.json.gz
```
The toy example json file contains a trajectory inferred by PAGA for a hematopoetic lineage.

3. Make sure you have python3 and pip3 installed (https://www.python.org/downloads/). Download dependencies
```bash
pip3 install numpy scipy matplotlib scikit-learn dtw statsmodels
```
4. Run DREAMIT (change the thread count based on your machine)
```bash
python3 DREAMIT.py -json_dictionary Data/paul_toy_data.json -tf_dict Data/trrust_rawdata.mouse.tsv -threads 20 -outdir test_toy_run
```
Note: Runtime is expected to be around 1 hour when using ~20 threads. Output may have minor variance due to spline fit, subsampling, and other factors.

# Parameters

```bash
python3 DREAMIT.py --help

options:
  -h, --help            show this help message and exit
  -json_dictionary JSON_DICTIONARY
                        input a JSON file with branching information
  -tf_dict TF_DICT      input must be a factor to all targets file with scores (could be edited)
  -threads THREADS      input the number of threads to use. Default = 1
  -number_of_bins NUMBER_OF_BINS
                        Input option for the number of bins to use (default = 20)
  -core_cut CORE_CUT    Input the percentile cutoff for the core set of targets (e.g. 25 cuts the bottom 25 percentile of targets in association with TF)
  -gam GAM              Choose "bins" or "nobins". "bins" is better tested performance.
  -bin_by_var BIN_BY_VAR
                        Option to bin by any string. Not developed, leave as default="psuedotime"
  -even_distribution EVEN_DISTRIBUTION
                        Option for creating bins with an equal distribution of cells in each bin. Not developed, (default = False)
  -outdir OUTDIR

```

# More info
TF-target data can be found here https://www.grnpedia.org/trrust/ or used directly from the Data/ folder. The format of any TF-target data should be a tsv file with the Transcription Factor listed first, its target second, the relationship third, and the weight fourth. For example:

```bash
AATF    BAX     Repression      22909821
AATF    CDKN1A  Unknown 17157788
AATF    KLK3    Unknown 23146908
AATF    MYC     Activation      20549547
```

The inferred trajectory solution should be provided through a JSON file.Once the user has a trajectory inferred to their liking they should format the data into a JSON file similar to that found in this toy example within the Data/ folder where each branch is its own key. The JSON file should be formatted to look like this example where the cell attributes can contain additional information, but must contain the "cell_id", "pseudotime", and the genes to be included in the analysis:

```bash
{"branches":{"Branch_1":{"cell_attribute1": [a,b,c,...,n], "cell_attribute2": [a,b,c,...,n], "cell_attributen": [a,b,c,...,n]}, "Branch_n":{"cell_attribute1": [a,b,c,...,n], "cell_attribute2": [a,b,c,...,n], "cell_attributen": [a,b,c,...,n]}}}
```

# Rights and Permissions
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
