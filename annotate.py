#====================================#
# Annotate_Platypus_VCF.py           #
# Version 1.0                        #
# Seunghun Han, Harvard University   #
# seunghunhan@g.harvard.edu          #
# Mar 22, 2023                       #
#====================================#


#Given raw vcf file generated using the Platypus variant caller (https://github.com/andyrimmer/Platypus), 
#This script annotate each variant with the following information:
# 1) Depth of sequence coverage 2) Number of reads supporting variants 3) Percentage of reads supporting the variant versus those supporting reference reads
# 4) Gene name  5) Variant Type 6) Variant Consequence 7) Minor Allele Frequency

#This script uses REST API for Variant Effector Predictor (https://rest.ensembl.org/#VEP)


##Libraries

import concurrent.futures
import requests, sys
from tqdm import tqdm
import pandas as pd
import numpy as np
import argparse
import json
import time
import os

# Python 2/3 adaptability
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError


#parse input parameters
parser = argparse.ArgumentParser(description='Annotate Platypus VCF')
parser.add_argument('--input_VCF', metavar='input_VCF', type=str, help= 'input variant call file', required= True)
parser.add_argument('--Genome_build', metavar='Genome_build', default='GRCh37', type=str, help = 'Genome build, default = GRCh37')
parser.add_argument('--output_dir', metavar='output_dir', type=str, default='./', help = 'Output Folder')
parser.add_argument('--number_thread', metavar='number_thread',type=int, default = 4 , help = "Number of threads to be used for multi-tasking")
parser.add_argument("--output_filename", metavar="output_filename",type=str, default = "VEP_annotated_variants", help = 'File name for the output')
parser.add_argument("--print_raw_cols", metavar="print_raw_cols",type=str, default ='no', help = "Decide whether to output the raw vcf annotated with information")
args = parser.parse_args()

input_VCF = args.input_VCF
output_dir = args.output_dir
Genome_build = args.Genome_build
number_thread = args.number_thread
output_filename = args.output_filename
print_raw_cols = args.print_raw_cols


## Functions

#Function to extract specific field from INFO column
def extract_INFO(row,ID):
    INFO=row['INFO']
    Split_INFO=INFO.split(";")
    Extracted_Field=[x for x in Split_INFO if ID in x][0]
    Extracted_Field=int(Extracted_Field.split("=")[1])
    return Extracted_Field


#Function for Decomposing Multiallelic Variants
def decompose_multiallelic(df_original,df_to_decompose):
    df_decomposed = pd.DataFrame(columns = df_original.columns)  

    for i, row in df_to_decompose.iterrows():
        alts = row['ALT'].split(',')
        alt_count=len(alts)
        
        #Get indices of relevant columns
        format_split = row['FORMAT'].split(":")
        NR_index=format_split.index("NR")
        NV_index=format_split.index("NV")
        GL_index=format_split.index("GL")
        GT_index=format_split.index("GT")
        
        #Decomposing multi-allelic per allele
        for alt in alts:
            new_row = row.copy()
            new_row['ALT']= alt
            
            #Decomposing INFO field 
            info_fields = row['INFO'].split(';')
            info_dict = {field.split('=')[0]: field.split('=')[1] for field in info_fields}
            for key in info_dict:
                if ',' in info_dict[key]:
                    info_dict[key] = info_dict[key].split(',')[alts.index(alt)]
            new_row['INFO'] = ';'.join([f"{key}={info_dict[key]}" for key in info_dict])
            
            #Decomposing sample field 
            samples_fields = row['sample'].split(':')
            samples_dict = {samples_fields.index(field):field for field in samples_fields}
            for key in samples_dict:
                if (GL_index!=key) & (',' in samples_dict[key]): #GL has values with "," regardless of alt count, which should not be split
                    samples_dict[key] = samples_dict[key].split(',')[alts.index(alt)]
            
            #Further decompose GT to 0/0, 0/1, or 1/1
            GT_value=samples_dict[GT_index]
            if alt_count>1:    #Only decompose GT of multiallelic variants
                alt_index=alts.index(alt)+1
                GT_fields=samples_dict[GT_index].split("/")
                GT_fields=[int(x) for x in GT_fields]
                gt_dict= {0:"0/0",1:"0/1",2:"1/1"}
                alt_allele_count=GT_fields.count(alt_index)
                samples_dict[GT_index]=gt_dict[alt_allele_count]
            new_row['sample'] = ':'.join([f"{samples_dict[key]}" for key in samples_dict])        
          
            df_decomposed=pd.concat([df_decomposed,pd.DataFrame(new_row).T])
           
            
    return df_decomposed

#Function for VEP annotation
class EnsemblVepAPI(object):
    if Genome_build=='GRCh37':
        server_url='https://grch37.rest.ensembl.org/vep/human/region'
    else:
        server_url='https://rest.ensembl.org/vep/human/region'   #For GRCh38 

    def __init__(self, server=server_url, reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
    
    def run_rest_API(self, endpoint, hdrs=None,params=None):
        if hdrs is None:
            hdrs = {}
        
        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'
        
        if params:
            endpoint += "?" + urlencode(params)
            
        data = None
        
        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        try:       
            
            request = Request(self.server + endpoint, headers=hdrs)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1
            
        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.run_rest_API(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
                
        return data
    
    def get_VEP_annotation(self, chroms, start, end, alt, ref):
        if len(ref)>1:
            end=end+len(ref)-1
        annots = self.run_rest_API(
            endpoint="/{0}:{1}-{2}/{3}".format(chroms,start,end,alt),
            params={"pick" : "1", "variant_class": "1", "allele_number":1}
        )
        
        VEP_id=annots[0]['id']
        if "transcript_consequences" in annots[0]:
            gene_name=annots[0]['transcript_consequences'][0]['gene_symbol']
            variant_consequence=annots[0]['transcript_consequences'][0]['consequence_terms'][0]           
            variant_class=annots[0]['variant_class']
        else:
            gene_name="Not Found"
            variant_consequence=annots[0]['most_severe_consequence']
            variant_class=annots[0]['variant_class']
        
        MAF="Not Available"
        
        #Assigning MAF only when the alt allele is the actual minor allele for the location
        if "colocated_variants" in annots[0]:
            if 'minor_allele_freq' in annots[0]['colocated_variants'][0]:
                if 'minor_allele' in annots[0]['colocated_variants'][0]:
                    if annots[0]['colocated_variants'][0]['minor_allele']==alt:
                        MAF=annots[0]['colocated_variants'][0]['minor_allele_freq']
        
        return gene_name, variant_consequence, variant_class, VEP_id, MAF

client = EnsemblVepAPI()

#Function for multi-threading
def process_subset(df_subset):
    df_subset[['SYMBOL', 'Consequence', 'VARIANT_CLASS','VEP_ID',"MAF"]] = df_subset.apply(
        lambda row: client.get_VEP_annotation(row['#CHROM'], row['POS'], row['POS'], row['ALT'], row['REF']),
        axis=1, result_type='expand')
    return df_subset

def parallel_apply(df, func, n_threads=4):
    # Split the DataFrame into subsets based on the number of threads
    subsets = np.array_split(df, n_threads)

    # Create a wrapper function that updates the progress bar
    def func_wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        progress_bar.update(len(args[0]))
        return result

    # Process each subset in parallel using ThreadPoolExecutor
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as executor:
        # Create a progress bar using tqdm
        with tqdm(total=len(df), desc="Processing", unit="row") as progress_bar:
            results = list(executor.map(func_wrapper, subsets))

    # Concatenate the results to obtain the final DataFrame
    return pd.concat(results)

#Function to write the output file    
def OutputFile(final_df, outputDir, filename):
    """Given a df, output the new file"""
    final_df.to_csv(os.path.join(outputDir,filename),"\t",index=False)


#Main Function to run the tool
def main():
    """RUNNING THE CODE"""
    
    ## Step 1 - Load the input vcf 
   
    #Find the index of the column header row 
    with open(input_VCF) as file:
        for index, line in enumerate(file):
            if line.startswith("#") and not line.startswith("##"):
                header_row_index = index
                break
    
    #Load the VCF without the ## headers, but keep the column names
    input_vcf = pd.read_csv(input_VCF, skiprows=lambda x: x < header_row_index,sep='\t')


    ## Step 2 - Pre-processing vcf before running VEP API

    #Flag multi-allelic variants 
    input_vcf['ALT_Type']="Biallelic" 
    input_vcf.loc[input_vcf.ALT.str.contains(","),'ALT_Type']='Multi-Allelic'

    #Decompose multi-allelic variants into individual variants
    decomposed_vcf=input_vcf.copy()
    Multi_indices=input_vcf[input_vcf.ALT.str.contains(",")].index

    #Run the decompose function only for the multi-allelic variants
    for i in range(0,len(Multi_indices)):
        multi_per_df=list(decomposed_vcf[decomposed_vcf.ALT.str.contains(",")].index)[0]
        row_to_duplicate = decomposed_vcf.iloc[[multi_per_df]]
        row_decomposed=decompose_multiallelic(decomposed_vcf, row_to_duplicate)
        decomposed_vcf = pd.concat([decomposed_vcf.iloc[:multi_per_df],row_decomposed ,decomposed_vcf.iloc[multi_per_df+1:]], ignore_index=True)
        decomposed_vcf=decomposed_vcf.reset_index(drop=True)

    #Add number of total and supporting reads
    decomposed_vcf['Read_Depth']=decomposed_vcf.apply(lambda row: extract_INFO(row,"TC="),axis=1)
    decomposed_vcf['Supporting_Reads']=decomposed_vcf.apply(lambda row: extract_INFO(row,"TR="),axis=1)

    #Add Variant Allele Frequency (VAF) and Reference Allele Frequency (RAF) as percentage
    decomposed_vcf['VAF']=decomposed_vcf['Supporting_Reads']/decomposed_vcf['Read_Depth']*100
    decomposed_vcf['RAF']=100-decomposed_vcf['VAF']

    ## Step 3 - Running VEP annotation

    print("\n")
    print(f'Running VEP annotation with {number_thread} cores\n')

    #Annotate the decomposed vcf with VEP API
    annotated_VCF = parallel_apply(decomposed_vcf, process_subset, number_thread)
    short_cols=['VEP_ID','SYMBOL','VARIANT_CLASS','Consequence','MAF','VAF','RAF','Read_Depth','Supporting_Reads','ALT_Type','sample']
    final_annotated_VCF=annotated_VCF[annotated_VCF.FILTER=="PASS"][short_cols]   #For final output, only including the variants with PASS under FILTER column

    #Write the output to a tsv file
    OutputFile(final_annotated_VCF, output_dir, "{}_processed.tsv".format(output_filename))

    print("\n")
    print(f"Output file = {output_dir}{output_filename}_processed.tsv")
    print("\n")

    #Also output the table with all raw columns if print_raw_cols is yes
    if print_raw_cols == "yes":
        OutputFile(annotated_VCF, output_dir, "{}_processed_raw.tsv".format(output_filename))

        print("\n")
        print(f"Output file                      = {output_dir}{output_filename}_processed.tsv")
        print(f"Output file with raw vcf columns = {output_dir}{output_filename}_processed_raw.tsv")
        print("\n")
    
    print("Annotation is complete\n")

if __name__ == '__main__':
    main()
