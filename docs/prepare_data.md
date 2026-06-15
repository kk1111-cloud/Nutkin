# Preparing the Test Data from Raw GEO Files

The test dataset `test_data.h5ad` was generated from raw scRNA-seq data originally obtained
from GEO (accession [GSM7872694](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7872694)).

If you wish to reproduce `test_data.h5ad` from scratch rather than using the file provided
in the `data/` directory, follow the steps below.

## Raw Input Files

Download the following files from GEO:

| File | Size |
|---|---|
| `GSM7872694_BARtab_cell-barcode-anno.tsv.gz` | 103.6 KB |
| `GSM7872694_barcodes.tsv.gz` | 68.5 KB |
| `GSM7872694_features.tsv.gz` | 245.0 KB |
| `GSM7872694_matrix.mtx.gz` | 145.7 MB |

## Processing Steps

The following steps were used to create `test_data.h5ad`:

1. **Load the raw count matrix** in MTX format and transpose it so that rows correspond
   to cells and columns to genes.

2. **Assign gene names and cell barcodes** using the accompanying feature and barcode files.

3. **Parse lineage barcode annotations** from `GSM7872694_BARtab_cell-barcode-anno.tsv.gz`
   to compute:
   - Number of barcodes detected per cell
   - Maximum UMI count per barcode

4. **Filter minor barcodes**: remove barcodes supported by fewer than half of the maximum
   UMI count within each cell.

5. **Merge** the processed count matrix and filtered BARtab annotations into a Scanpy
   AnnData object. The resulting object contains the count matrix, gene metadata
   (`gene_ids` and gene symbols), and cell-level BARtab annotations aligned to `adata.obs`.

6. **Subset to four clones** to reduce dataset size for the test example:
   - `mCHERRY_Barcode_1614`
   - `mCHERRY_Barcode_5774`
   - `mCHERRY_Barcode_1755`
   - `mCHERRY_Barcode_65372`

7. **Save** the filtered AnnData object as `test_data.h5ad`.

## Notes

- The full unsubsetted dataset is considerably larger. The four-clone subset used here
  contains 1,932 cells.
- The `clone_id` column in `adata.obs` identifies which clone each cell belongs to and
  is used as the `group_col` throughout this manual.
