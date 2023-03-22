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

For item number 4-7, the tool uses REST API for Variant Effector Predictor (https://rest.ensembl.org/#VEP)



Installation
-------------
Click the green Clone or Download button above, copy the presented web URL, and then paste it into your terminal, preceded by the command 'git clone'.
  `git clone https://github.com/seunghunh23/annotate_Platypus_VCF.git`
  
  
  
  
  
  
Parameters
-------------

* `input_VCF` : Name of the input vcf file. The input file needs to be inside the sample folder as the annotate.py file
*  `Genome_build` : Genome build used when calling the variants. User can choose between `GRCh38` or `GRCh37` Default is GRCh37
*  `output_dir` : Output directory to write the annotated tsv file. Default is the current folder
*  `number_thread` : Number of thread to use for multi-threading. Default is 4 cores
*  `output_filename` : Filename for the annotated tsv file. 
*  `print_raw_cols` : Outputs a table with all raw vcf columns on top of the default shortened annotated output if "yes". Default is "no" 


Usage
--------------
`python annotate.py --input_VCF test_vcf_data.txt --output_dir results --number_thread 8 --output_filename test`

Output
--------------
Example default annotated variant output row:



|VEP_ID	|SYMBOL|	VARIANT_CLASS|	Consequence|	MAF|	VAF	|RAF	|Read_Depth	|Supporting_Reads	|ALT_Type	|sample	
|---|---|---|---|---|---|---|---|---|---|---|
 |1_1895177_C/T|	C1orf222|	SNV|	intron_variant|	0.1665|	99.074074|	0.925926|	108|	107|	Biallelic|	1/1:-300.0,-29.5,0.0:2:99:108:107|	
 
 
 * `VEP_ID` : Variant ID in the format of CHROM_POS_REF/ALT
 * `SYMBOL` : Gene name 
 * `VARIANT_CLASS` : Type of the variant. (SNV, indel, substitution)
 * `Consequence` : Effect of the variant. (missense_variant, frameshift_variant, start_lost...etc).    
       Refer to https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
 * `MAF` : Minor Allele Frequency of a variant. This value is given only when the alt allele of the variant is actually a minor allele for the position
 * `VAF` : Variant Allele Frequency. This is number of reads supporting the variant divided by the total number of sequencing reads at the position
 * `RAF` : Reference Allele Frequency. Number of reads supporting the ref allele divided by the total number of sequencing reads (Also `100 - VAF`)
 * `Read_Depth` : Total number of reads at the site of variant (Depth of sequence coverage).
 * `Supporting_Reads` : Total number of reads supporting the variant
 * `ALT_Type` : Indicate whether a variant is `Multi-Allelic` or `Biallelic`. Variants with Multi-Allelic flag are the ones originally had one row in the vcf with more than one possible ALT allele. They are decomposed into individual variant per alt allele. We recommend users to review the raw sequencing file (BAM or CRAMS) for the multi-allelic variants
 * `sample` : The same sample column value extacted from the original VCF. 
