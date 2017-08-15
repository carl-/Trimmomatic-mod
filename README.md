# Trimmomatic-mod
A modified version of trimmomatic; 

original version could be found at http://www.usadellab.org/cms/index.php?page=trimmomatic

## Modifications

- [x] 1. add support to search adapter files in provided folders when file not found in cwd; [ Completed ]

- [x] 2. future plan: add ncbi ngs sdk to support processing ability for sra file; [ Completed ]

## HOWTO

#### Adapter
adapter files come with trimmomatic could be used without specifying fullpath, the version of trimmomatic are able to search the adapter folder where the jar file of trimmomatic locates


see branch 'auto_adapters'

#### NCBI SRA file support
 NCBI NGS SDK (ngs-sdk and ngs-java) and NCBI VDB (ncbi-vdb) should be installed to use this feature

see branch 'ngs-java'
