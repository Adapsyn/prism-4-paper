## prism-4-paper

Code used to conduct the analyses in Skinnider et al., "Comprehensive prediction of secondary metabolite structure and biological activity from microbial genome sequences."

Data is available from Zenodo:
- [10.5281/zenodo.3985982](http://dx.doi.org/10.5281/zenodo.3985982): 'Gold standard’ database of 1,281 BGCs linked to known secondary metabolites
- [10.5281/zenodo.3985977](http://dx.doi.org/10.5281/zenodo.3985977): PRISM output for all 10,121 genomes analyzed

Please also see (https://prism.adapsyn.com/)[https://prism.adapsyn.com/], and the tutorials there, for further details.

### Details

The code has been tested with python (version 3.7.6), and requires the following python libraries:

```
collections
gzip
itertools
joblib
json
math
numpy
os
pandas
pickle
random
rdkit
re
sklearn
sys
tqdm
```

To download a copy of the code, execute the following command in a terminal:

```
git clone https://github.com/Adapsyn/prism-4-paper
```

This should take no more than a few seconds.

Details regarding how to run each python script are provided within the comments and documentation of each file. Run time may range from seconds to several hours, depending on the script in question and the details of the user's hardware, but will allow the interested user to reproduce the analyses in the manuscript.

As an example, to extract all the biosynthetic gene clusters from a set of antiSMASH 5 results and write them to a comma-delimited file, enter the paths to the antiSMASH results on your local filesystem into the `io` tuple, then run:

```python extract-antismash-clusters.py```
