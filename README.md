The TSS meta table is available for [download here](https://s3.msi.umn.edu/zhoup-cage/91_share/01.rtss.xlsx), which has the following columns:
- `tss_id`: TSS ID
- `tpm`: pooled mean TPM from all 34 samples used to call TSS clusters
- `shapeIQR`: interquartile width of the TSS tag count vector
- `shapeEntropy`: Shannon Entropy (base log2) of the tag count vector, `shape_index (SI) = 2 - shapeEntropy`
- `support`: number of expressed samples (TPM>=1)
- `tid`, `gid`: AGPv5 transcript and gene ID
- `peakType`: genomic feature of the TSS location
- `npos`: number of unique tag positions in the TSS cluster
- `chrom`, `start`, `end`, `srd`: genomic range of the TSS cluster
- `peakPos`: coordinate of the peak tag position ("dominant TSS")
- `vtype`: B73-A632 variant information of this region (actually, this region plus downstream 20bp), can be one of:
  - `identical`: the A632 sequence of this region is identical to that of B73
  - `1snp`: only one SNP difference between the A632 sequence and B73 sequence for this region
  - `2snp`: only two SNP difference
  - `other`: all other scenarios (more than 2 SNPs, any number of InDels or SVs, poorly covered A632 region, etc.)
- `ibd`: whether this region locate in one of the 42 IBD regions between B73 and A632
- `prop_tpm`, `dominant`: relative proportion of TPM contribution to all TSS clusters identified within a multi-TSS gene, called "dominant" if this proportion exceeds 15%
- `tpm.shoot.B`, `tpm.shoot.A`, ..., `tpm.root.H`: pooled mean TPM in 6 conditions: shoot/root B73/A632/Hybrid
- `ss-root.B-shoot.B`, `ss-root.A-shoot.A`, ..., `ss-shoot.H-shoot.C`: shifting scores calculated for 10 pairwise comparisons
  - this value is only available when both conditions of the pair has TPM>=1
  - "shoot.C" means 1:1 mixed TPM vector from "shoot.B73" and "shoot.A632" (i.e., mid-parent expectation of F1)


Bigwig and BAM files of the maize B73-A632 CAGE project can be downloaded:
[https://s3.msi.umn.edu/zhoup-cage/CAGE_tracks_AGPv5.tar.gz](https://s3.msi.umn.edu/zhoup-cage/CAGE_tracks_AGPv5.tar.gz)
- (~90GB; can be downloaded using `wget` and uncompressed using `tar -xzf`)

Directory structure is as follows:

```
CAGE_tracks_AGPv5/
├── 00.meta.tsv
├── 01_bam
│   ├── husk.B73.01.bam
│   ├── husk.B73.01.bam.csi
│   ├── husk.B73.02.bam
│   ├── husk.B73.02.bam.csi
│   ├── ...
├── 01_bigwig
│   ├── husk.B73.01.minus.bw
│   ├── husk.B73.01.plus.bw
│   ├── husk.B73.02.minus.bw
│   ├── husk.B73.02.plus.bw
│   ├── ...
├── 05_bam_merged
│   ├── root.A632.bam
│   ├── root.A632.bam.bai
│   ├── root.B73.bam
│   ├── root.B73.bam.bai
│   ├── root.F1.bam
│   ├── root.F1.bam.bai
│   ├── shoot.A632.bam
│   ├── shoot.A632.bam.bai
│   ├── shoot.B73.bam
│   ├── shoot.B73.bam.bai
│   ├── shoot.F1.bam
│   └── shoot.F1.bam.bai
└── 05_bigwig_merged
    ├── root.A632.minus.bw
    ├── root.A632.plus.bw
    ├── root.B73.minus.bw
    ├── root.B73.plus.bw
    ├── root.F1.minus.bw
    ├── root.F1.plus.bw
    ├── shoot.A632.minus.bw
    ├── shoot.A632.plus.bw
    ├── shoot.B73.minus.bw
    ├── shoot.B73.plus.bw
    ├── shoot.F1.minus.bw
    └── shoot.F1.plus.bw
```

- `00.meta.tsv`: sample meta table 
- `01_bam`: sample-wise BAM file (e.g., `shoot.B73.01.bam`)
- `01_bigwig`: sample-wise TSS tag coverage (bigwig) file for plus and minus strands (e.g., `shoot.B73.01.plus.bw` and `shoot.B73.01.minus.bw`)
- `05_bam_merged`: replicate-merged BAM file (e.g., `shoot.B73.bam`)
- `05_bigwig_merged`: replicate-merged bigwig file (e.g., `shoot.B73.plus.bw` and ``shoot.B73.minus.bw``)

Note: the custom reference fasta and gff files used here are a bit different from the standard AGPv5 reference (chromosome names being 'chr01','chr02', ..., 'chr10' rather than '1', '2', ..., '10') and can be downloaded here:
- [https://s3.msi.umn.edu/zhoup-genome/Zmays_B73v5/10.fasta](https://s3.msi.umn.edu/zhoup-genome/Zmays_B73v5/10.fasta)
- [https://s3.msi.umn.edu/zhoup-genome/Zmays_B73v5/50_annotation/10.gff](https://s3.msi.umn.edu/zhoup-genome/Zmays_B73v5/10.fasta)

