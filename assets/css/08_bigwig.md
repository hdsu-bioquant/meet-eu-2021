# Generating signal files

So far, we have worked with the following files

* `fastq` files, which represent raw reads;
* `bam` files, which represent aligned reads;
* `bed` or  `narrowPeak` files, which represent discrete  genomic regions containing signal.

In a final step, we will generate *signal files* in the `bigwig` format, which represent continuous signals along the genome. Another format to do this is the `bedgraph` format. For both formats, the files can then be loaded into a genomic browser in order to explore the signal and the peak regions.

For ChIP-seq, we have the **IP** file representing the signalm, and the **control** file (generally input), which represents the level of background noise. When generating the signal file, we can (1) either generate separate bigwig files for IP and control, or (2) generate a single bigwig file, in which we substract the background noise from the IP.

## Generating separate bigwig

We will use a tool form the `bedtools` suite, namely `genomecov`, to build the bigwig files. We need to proceed in 2 steps: (1) create a bedgraph file, (2) convert bedgraph to bigwig

```
cd

## create a directory
mkdir -p myanalysis/bigwig

# create and sort bedgraph file for the CTCF IP
bedtools genomecov -bg -ibam /vol/volume/HCT116/analysis/CTCF/Bowtie2/CTCF_Rep1_ENCFF001HLV_trimmed_aligned_filt_sort_nodup.bam | sort -k1,1 -k2,2n > myanalysis/bigwig/CTCF_Rep1.bg

# convert bedgraph to bigwig
bedGraphToBigWig myanalysis/bigwig/CTCF_Rep1.bg /vol/volume/HCT116/data/hg38.genome myanalysis/bigwig/CTCF_Rep1.bw

# delete the bedgraph file
rm myanalysis/bigwigCTCF_Rep1.bg
```

### Exercice

> create the bigwig file for the control `/vol/volume/HCT116/analysis/CTCF/Bowtie2/CTCF_Control_ENCFF001HME_trimmed_aligned_filt_sort_nodup.bam`; MAKE SURE TO GIVE A DIFFERENT NAME TO THE OUTPUT FILES!!

## Generating a merged bigwig file

We can also create a single bigwig file by substracting the Control from the IP; scaling Control and IP is highly non-trivial, a good method is the **SES** method implemented in DeepTools

```
bamCompare -b1 /vol/volume/HCT116/analysis/CTCF/Bowtie2/CTCF_Rep1_ENCFF001HLV_trimmed_aligned_filt_sort_nodup.bam \
 -b2 /vol/volume/HCT116/analysis/CTCF/Bowtie2/CTCF_Control_ENCFF001HME_trimmed_aligned_filt_sort_nodup.bam \
 --scaleFactorsMethod SES -p 4 \
 --operation subtract \
  -o myanalysis/bigwig/CTCF_ses_subtract.bw
```

## Loading the files into the IGV app       

We can now use the [IGV web app](https://igv.org/app/) to visualize the different datasets into a genomic browser; here, we will load
* the individual bigwig files for IP and control
* the compbined bigwig file
* the peak files in narrowPeak format

Procedure

1. go to to the [IGV web app](https://igv.org/app/)
2. in the top left menu *Genome*, select the `hg 38` genome version
3. using Cyberduck, download the following files onto your local computer
   * `CTCF_ses_subtract.bw`
   * `CTCF_Rep1.bw`
   * `CTCF_peak.narrowPeak` (should be in the `myanalysis/MACS2/CTCF` folder!)
4. in the IGV app, go to `Tracks > local files...` and load the bigwig and narrowPeak file
5. you can now select a specific chromosome in the drop-down menu, or type in a gene name to see this gene locus.

### Question

> do the called peaks correspond to regions of higher signal in the bigwig tracks?
> choose two peaks in the bigwig file which seem to have different height; click on the corresponding peak in the narrowPeak track (blue rectangle) to see the P-values and Q-values; do you see a difference?
