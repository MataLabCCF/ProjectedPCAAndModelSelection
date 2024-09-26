# ProjectedPCAAndModelSelection

Script with the goal to automatize the Projected PCA and the model selection using stepwise regression. We export all covar (all covar + all PCAs) and de stepwise regression 

## Parameters

```
PCA and regression

options:
  -h, --help            show this help message and exit

Data arguments:
  -A AUTOSOMAL, --autosomal AUTOSOMAL
                        Genotyped file name for autosomal chromosome
  -t TABLECOVAR, --tableCovar TABLECOVAR
                        File with covariatives to be added to the model
  --threads THREADS     Number of processors to be used (default = 1)

Reference data arguments:
  -R AUTOSOMALREF, --AutosomalRef AUTOSOMALREF
                        Data with parental reference data to run autosomal PCA (optional)

Output arguments:
  -n NAME, --name NAME  Analysis name
  -f FOLDER, --folder FOLDER
                        Folder to output files

PCA arguments:
  --model MODEL [MODEL ...]
                        Regression model. If you do not provide a model the script willuse the selectModel.R to build the model

Programs:
  --plink2 PLINK2       Path of PLINK 2 (default = plink2)
  --plink1 PLINK1       Path of PLINK 1 (default = plink)
  --gcta GCTA           Path of gcta
  --selectModel SELECTMODEL
                        Path of selectModel script
```

## Example of command line
```
python covarProjectedPCA.py \
   -A /home/peixott/beegfs/Analysis/DataClean/CleanData/LARGE_NewPipeline/FinalData/LARGE_Phase2_QCed_Autosomal \
   -t covar.txt -R /home/peixott/beegfs/Analysis/DataClean/CleanData/Shriner/OneThousand_All \
   -n LARGE_TryNew -f ./NewPCs_OutModel \
   --gcta /home/peixott/beegfs/Programs/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 \
   --selectModel ./selectModel.R --plink1 plink --plink2 /home/peixott/beegfs/Programs/plink2
```
