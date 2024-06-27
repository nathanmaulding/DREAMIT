# DREAMIT
Here we present DREAMIT (Dynamic Regulation of Expression Across Modules in Inferred Trajectories). DREAMIT aims to analyze dynamic regulatory patterns along trajectory branches, implicating transcription factors (TFs) involved in cell state transitions within scRNAseq datasets. DREAMIT uses pseudotime ordering within a robust subrange of a trajectory branch to group individual cells into bins. It aggregates the cell-based expression data into a set of robust pseudobulk measurements containing gene expression averaged within bins of neighboring cells. It then smooths trends after searching for an optimal fitting spline across the bins. DREAMIT rejects further analyzing branches that produce highly variable smoothing estimates (covariation in spline fitting parameters found to be greater than 1.0 across 80% subsampling) as these branches may represent sparse or noisy parts of the data that could produce unreliable TF inferences. 

Using the transformed smoothed data, it calculates the association between a TF and all of its predicted targets according to the TRRUST database assessed using multiple metrics (e.g., Pearson correlation, Mutual Information, Dynamic Time Warping distance, etc). Finally, DREAMIT uses a Relational Set Enrichment Analysis (RSEA) test to evaluate the significance of the TF-to-target associations and identify a core set of targets compared to a background model, which consists of arbitrarily selected targets.

# Tutorial
For toy example of running DREAMIT follow these steps:

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
pip3 install numpy scipy matplotlib scikit-learn dtw-python statsmodels
```
4. Run DREAMIT (change the thread count based on your machine)
```bash
python3 DREAMIT.py -json_dictionary Data/paul_toy_data.json -tf_dict Data/trrust_rawdata.mouse.tsv -threads 20 -outdir test_toy_run
```

# Parameters

```python3
self.parser.add_argument('-json_dictionary',
                                 help="input a JSON file with branching information (created by NM_slingshot_v2.R)")
self.parser.add_argument('-tf_dict',
                         help="input must be a factor to all targets file with scores (could be edited)")
self.parser.add_argument('-threads', type=int, default=1,
                         help='input the number of threads to use. Default = 1')
self.parser.add_argument('-number_of_bins', type=int, default=None,
                         help='Input option for the number of bins to use (default = 20)')
self.parser.add_argument('-core_cut', type=int, default=25,
                         help='Input the percentile cutoff for the core set of targets (e.g. 25 cuts the bottom 25 percentile of targets in association with TF)')
self.parser.add_argument('-gam', type=str, default='bins',
                         help='Choose "bins" or "nobins". "bins" is better tested performance.')
self.parser.add_argument('-bin_by_var', type=str, default="psuedotime",
                         help='Option to bin by any string. Not developed, leave as default="psuedotime"')
self.parser.add_argument('-even_distribution', type=bool, default=False,
                         help='Option for creating bins with an equal distribution of cells in each bin. Not developed, (default = False)')
self.parser.add_argument('-outdir', type=str, default='./dreamit_analysis')
```


For toy example of running DREAMIT download Data/paul_toy_data.json.gz and Data/trrust_rawdata.mouse.tsv (which can also be found here https://www.grnpedia.org/trrust/)

The inferred trajectory solution is contained in the JSON file. Once the user has a trajectory inferred to their liking they should format the data into a JSON file similar to that found in this toy example where each branch is its own key. The toy example contains a trajectory inferred by PAGA for a hematopoetic lineage.
From the terminal, unzip the toy data file with:

gunzip path/to/paul_toy_data.json.gz

Also download the DREAMIT.py file along with its dependencies.

Run DREAMIT on the toy example with the following command with the appropriate paths to the script, data, and desired output folder location and name:

python3 DREAMIT.py -json_dictionary path/to/paul_toy_data.json -tf_dict path/to/trrust_rawdata.mouse.tsv -threads 20 -outdir path/to/testtoy

Note: Runtime is expected to be around 1 hour when using ~20 threads. Output may have minor variance due to spline fit, subsampling, and other factors.
