GB-Score 
--
A scoring function based on Gradient Boosting Trees algorithm for predicting ligand-protein binding affinity. PDBbind 2019v collection of the general and refined sets minus core set is used for training GB-Score and core set is used as an independent test set.

![](https://github.com/miladrayka/GB_Score/blob/main/Capture300.JPG)

Contact
---
Milad Rayka, milad.rayka@yahoo.com

Citation
--
[GB-Score: Minimally Designed Machine Learning Scoring Function Based on Distance-weighted Interatomic Contact Features](https://onlinelibrary.wiley.com/doi/10.1002/minf.202200135)

Installation
--
Below packages should be installed for using GT-Score. Dependecies:

- python = 3.8.8

- numpy = 1.21.2

- pandas = 1.2.4

- seaborn = 0.11.1

- joblib = 1.0.1

- matplotlib = 3.3.4

- biopandas = 0.2.8

- scipy = 1.7.1

- sklearn = 0.24.1

- progressbar2 = 3.53.1

For installing first make a virtual environment and activate it.

On windows:

>    python py -m venv env
>    .\env\Scripts\activate

On macOS and Linux:

>    python3 -m venv env
>    source env/bin/activate

Which env is the location to create the virtual environment. Now you can install packages:

>    pip install *package_name*

Usage
--
1- Preparing ligand and protein file

a- Ligand and protein structure should be saved in .mol2 and .pdb format files respectively.
b- Each ligand and protein files for a specific complex must be placed in a same folder.

for example:

>     ./1a1e/1a1e_ligand.mol2
>     ./1a1e/1a1e_protein.pdb
>     ./1a4k/1a4k_ligand.mol2
>     ./1a4k/1a4k_protein.pdb

2- Generating features

generate_features.py is used for generating features of GB-Score and return a csv file:

>    python generate_feature.py -h  
>    python generate_feature.py -d file_directory  

3 - Repeat and extend current report

Jupyter notebook (analysis.ipynb) is provided to reproduce all results and figures in the current report.  

Generated features and saved model
--

All generated features (.csv), best saved model (.joblib), and PDBIDs (.txt) of all train an test sets can be found in files, saved_model, and pdbids folders, respectively.
