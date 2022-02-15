# nf-ninjamap2

## Usage

### Test Run

```{bash}
aws batch submit-job \
    --profile maf \
    --job-name nf-nm2-0201-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-ninjamap2,\
"--manifest","s3://genomics-workflow-core/Results/ninjamap2/00_TEST/00_manifest/test.csv",\
"--project","00_TEST"
```

### Selected Strain Dropouts

```{bash}
aws batch submit-job \
    --profile maf \
    --job-name nf-nm2-0201-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=s3://nextflow-pipelines/nf-ninjamap2,\
"--manifest","s3://genomics-workflow-core/Results/ninjamap2/20220201_Strain_Dropout_Verification/00_manifest/dropouts_in_question.manifest.csv",\
"--project","20220201_Strain_Dropout_Verification"
```
