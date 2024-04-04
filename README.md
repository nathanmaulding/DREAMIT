For toy example of running DREAMIT download data/paul_toy_data.json.gz and data/trrust_rawdata.mouse.tsv (which can also be found here https://www.grnpedia.org/trrust/)

From terminal unzip the toy data file with:

gunzip path/to/paul_toy_data.json.gz

Also download the DREAMIT.py file along with its dependencies.

Run DREAMIT on the toy example with the following command with the appropriate paths to the script, data, and desired output folder location and name:

python3 DREAMIT.py -json_dictionary path/to/paul_toy_data.json -tf_dict path/to/trrust_rawdata.mouse.tsv -threads 20 -outdir path/to/testtoy

Note: Runtime is expected to be around 1 hour when using ~20 threads. Output may have minor variance due to spline fit, subsampling, and other factors.
