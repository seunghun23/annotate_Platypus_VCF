# Annotate_Platypus_VCF
Given a raw vcf file generated using Platypus variant caller (https://github.com/andyrimmer/Platypus), 
annotate each variant with 
1) Depth of sequence coverage 
2) Number of reads supporting variants 
3) Percentage of reads supporting the variant versus those supporting reference reads
4) Gene name
5) Variant Type 
6) Variant Consequence 
7) Minor Allele Frequency 

For items number 4-7, the tool uses REST API for Variant Effector Predictor (https://rest.ensembl.org/#VEP)



Installation
-------------
Click the green Clone or Download button above, copy the presented web URL, and then paste it into your terminal, preceded by the command 'git clone'.
  `git clone https://github.com/seunghunh23/annotate_Platypus_VCF.git`
  
  
  
  
  
  
Parameters
-------------

* `input_VCF` : Name of the input vcf file. The input file should be located in the sample folder as the annotate.py file
*  `Genome_build` : Genome build of the reference used for variant calling. Users can choose between `GRCh38` or `GRCh37` The default is GRCh37
*  `output_dir` : Output directory for the annotated TSV file. The default is the current folder
*  `number_thread` : Number of threads to use for multi-threading. The default is 4 cores
*  `output_filename` : Filename for the annotated TSV file. 
*  `print_raw_cols` : If set to "yes", the output will include a table with all raw VCF columns in addition to the default shortened annotated output. The default is "no".


Usage
--------------
`python annotate.py --input_VCF test_vcf_data.txt --output_dir results --number_thread 8 --output_filename test --print_raw_cols yes`

Output
--------------
<b>The default annotated output includes only the variants that passed the filtering based on the information in the FILTER column of the
  raw vcf.</b>
  To obtain the full list of variants along with the original columns from the input VCF, use the option `--print_raw_cols yes`

Example default annotated variant output row:



|VEP_ID	|SYMBOL|	VARIANT_CLASS|	Consequence|	MAF|	VAF	|RAF	|Read_Depth	|Supporting_Reads	|ALT_Type	|sample	
|---|---|---|---|---|---|---|---|---|---|---|
 |1_1895177_C/T|	C1orf222|	SNV|	intron_variant|	0.1665|	99.074074|	0.925926|	108|	107|	Biallelic|	1/1:-300.0,-29.5,0.0:2:99:108:107|	
 
 
 * `VEP_ID` : Variant ID in the format of CHROM_POS_REF/ALT
 * `SYMBOL` : Gene name 
 * `VARIANT_CLASS` : Type of the variant. (SNV, indel, substitution)
 * `Consequence` : Effect of the variant. (missense_variant, frameshift_variant, start_lost...etc).    
       Refer to https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
 * `MAF` : Minor Allele Frequency of a variant. This value is provided only when the alt allele of the variant is the minor allele for the position.
 * `VAF` : Variant Allele Frequency, calculated as the number of reads supporting the variant divided by the total number of sequencing reads at the position.
 * `RAF` : Reference Allele Frequency, calculated as the number of reads supporting the ref allele divided by the total number of sequencing reads (also 100 - VAF).
 * `Read_Depth` : Total number of reads at the site of the variant (Depth of sequence coverage).
 * `Supporting_Reads` : Total number of reads supporting the variant
 * `ALT_Type` : Indicates whether a variant is Multi-Allelic or Biallelic. Multi-Allelic variants are those with more than one possible ALT allele in the original VCF. They are decomposed into individual variants per alt allele. For example, if a variant has value of 1/2 for GT column in the vcf, this variant is divided into separate entities each with GT of 0/1. It is recommended that users review the raw sequencing files (BAM or CRAM) for multi-allelic variants.
 * `sample` : The sample column value extacted from the original VCF. The format is `GT:GL:GOF:GQ:NR:NV`


Note on VEP annotation
------------------
By default, VEP annotates every genomic feature (e.g. transcripts) with which a variant colocalizes. To obtain annotations for only the most relevant transcript per variant, the '--pick' option was utilized. For a comprehensive explanation of the selection criteria, please visit https://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick.
