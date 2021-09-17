
# Description

## Script for simulation in SLiM-3.6

**Divergent selection:** ```steppingStone.secondaryContact.divergentSelection.txt```

**Recurrent selection:** ```steppingStone.secondaryContact.recurrentSweeps.txt```

## Post-processing of .tree files produced from SLiM

**Impose neutral mutations on samples and generate VCFs:**

  ```treeseq2vcf.py```
  
**Run local ancestry estimation software ELAI on output VCFs:**

  ```master.sh```, which sub-launches:
  
  ```vcf2bimbam.py```
      
  ```bimbam2elai.sh```
      
**Compute entropy on local ancestry:**

  ```get_entropy_by_segments.jl```
  
**Compute the test statistic (Pearson's correlation coefficient) and test for significance:**

  ```get_correlations.jl```

