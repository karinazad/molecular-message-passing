# :pill: :envelope_with_arrow: Molecular Message Passing
Self-attention based message passing graph neural network for analysis of molecular graphs


### Task
structure-bioactivity/structure–property relationships
(QSAR/QSPR) of compounds

### Background


### Dataset
Datasets of molecular lipophilicity (4200 molecules, CHEMBL) and aqueous solubility (OCHEM, 1311)

#### Data Processing
From the dataset, remove duplicates and items that are not recognized by RDKit.
SMILES representations were converted to directed graphs (note: using Deepchem and Chemprop's MPN encoder). Representation of directed graphs included lists of:

1. Node to Edge 
2. Edge to Node 
3. Edge to Reverse Node
4. Node to Neighboring Nodes 


Tenfold stratified cross-validation (80%, 10%, 10%).
