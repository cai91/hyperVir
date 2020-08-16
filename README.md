hyperVir
========

Tool to generate amino acid hypervariation plots

<b>Requirements:</b>
* Python (tested v3.6.7)
* Numpy (tested v1.17.0)
* Scipy (tested v1.3.1)
* MAFFT (tested v7.453)
* Matplotlib (tested v3.0.2)


<b>Usage:</b> 
```
hyperVir.py <input_file>
```

input_file: FASTA file containing a list of related proteins <br />

<b>Output:</b>
```
<input_file>.out: Hypervariation signal across multiple sequence alignment
<input_file>.png: Hypervariation plot
