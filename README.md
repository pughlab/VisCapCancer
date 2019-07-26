# VisCapCancer
VisCapCancer is a tool to infer somatic copy number alterations in tumours from targeted sequencing data. VisCap calculates the fraction of overall sequence coverage assigned to genomic intervals and computes log2 ratios for tumours with respect to a panel of normals. Candidate somatic CNVs are called when log2 ratios exceed user-defined thresholds.

```
Rscript VisCapCancer.R VisCapCancer.cfg   /path/to/normals /path/to/tumours
```
