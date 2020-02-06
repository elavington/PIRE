# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:38:17 2020

@author: Erik Lavington


This program caculates Nucleotide diversity (pi) from SNP frequency data. 
Expected input is an Excel or csv file with data from varscan plus metadata.
A metadata column is used to calculate pi by subset from the same data file.
ALL SNVs WILL BE USED IN THE CALCULATION. However, multialleleic sites (more 
than biallelic) are fine. The same reference allele is expected per site and 
one alternate allele is expected per row. 

Returns pi and variance of pi for the queried region, per group
"""

#%%
import pandas as pd
import sys, getopt
#%%
def usage():
	#command line usage message/help
	print("""
	usage: python pi.py -i INPUTFILENAME [options]
	   
	-h, --help	                this help dialog
	
	REQUIRED:
	-i, --input <filename>    	filename to be called for input
	
	OPTIONAL:
	-o, --output <filename>    	filename for output file (Default: prints to stdout when no filename is provided)
	-p, --position <column name> 	column name containing position information (default: -p pos)
	-g, --grouping <column name> 	column name containing group identifiers (default: '-g index)
	-r, --reference <column name> 	column name of reference allele (default: -f ref)
	-a, --alternate <column name> 	column name of alternate allele (default: -a alt)
	-c, --coverage <column name> 	column name of total read depth (default: -c DP)
	-f, --frequency <column name> 	column name of SNP frequency (default: -f freqProp)
	-l, --length <integer> 	 	length of region to calculate pi, default: maximum value in positions
	-s, --per-site 	     	   	if used, calculate per-site pi instead of region pi (default)
	-z, --sample-size-correction 	if used, calculate sample size correction for pi (default: no correction)
	
	""")
	return 
#%%
def get_args(argv):
	
	options = {
			'group':'index',
			'reference':'ref',
			'position':'pos',
			'alternate':'alt',
			'coverage':'DP',
			'frequency':'freqProp',
			'perSite':False,
			'sizeCorrection':False,
			'infileType':'excel',
			'length':0,
			'outputfile':''
			}
	
	infilename = ''
	
	try:
	   opts, args = getopt.getopt(argv[1:],"hi:o:g:p:r:a:c:f:szl:",
			   ["help","input=","output=","grouping=","position=","reference=",
		       "alternate=","coverage=","frequency=","length=","per-site",
			   "sample-size-correction"])
	   
	except getopt.GetoptError as err:
		print("\n<<<",err,">>>")
		usage()
		sys.exit(2)
		
	for opt, arg in opts:
		if opt in ('-h','--help'):
			usage()
			sys.exit()
			
		elif opt in ("-i", "--input"):
			infilename = arg
			
		elif opt in ("-o", "--output"):
			options['outputfile'] = arg
			
		elif opt in ("-g", "--grouping"):
			options['group'] = arg
			
		elif opt in ("-r", "reference"):
			options['reference'] = arg
			
		elif opt in ("-p", "position"):
			options['position'] = arg
			
		elif opt in ("-a", "--alternate"):
			options['alternate'] = arg
			
		elif opt in ("-c", "coverage"):
			options['coverage'] = arg
			
		elif opt in ("-f", "frequency"):
			options['frequency'] = arg
			
		elif opt in ("-s", "per-site"):
			options['perSite'] = True
			
		elif opt in ("-z", "sample-size-correction"):
			options['sizeCorrection'] = True
		
		elif opt in ("-l", "length"):
			options['length'] = arg			
	
	data = open_data(infilename)
	
	return data, options
#%%
def open_data(filename):

	try:
		data = pd.read_excel(filename,header=0)
	
	except:
		try:
			data = pd.read_csv(filename,header=0)
		except:
			print("\nUnsupported input file type or no input file name provided\n")
			sys.exit(2)
	return data	
#%%	
def var_pi(pi,n):	
	#variance of pi from equation 27 in Fu and Li 1992
	varPi = (pi*((n+1)/(3*n-3)))+((pi**2)*((2*(n**2+n+3))/(9*(n**2-n))))

	return varPi

#%%
def site_pi(versions):
	#caculate pi from a single site; DO NOT multiply by 2!
	#takes in a dictionary with A,G,C,and/or T as keys and frequencies (0.0-1.0) as values 
	
	pi = 0
	for i in versions:
		for j in versions:
			if i != j:
				pi += versions[i]*versions[j]	
	return pi
#%%
def pi(pi_data,sites,pos_id,Reference,Alternate,frequency,correction,N_samples,per_site):
	#calculate pi per site from snp frequencies
    #read in a pandas DataFrame object of major or minor allele frequencies,
    #one frequency per row
    #assumes SNPs only, no indels;does NOT assume biallelic sites
    #calculate pi for the region, NOT per-site pi
    #optional sample size correction
    
	subset_data = pi_data.loc[pi_data[pos_id]>=sites[0]].loc[pi_data[pos_id]<sites[1]]
	
	pi = 0
	#extract rows from the ps.DataFrame "data" for this site
	for site in subset_data[pos_id].unique():
		
		site_data = subset_data.loc[subset_data[pos_id]==site]
		#set dictionary of versions and frequency of each: alternates declared, reference calculated
		versions={}
		for i in site_data.index:
			   versions[site_data.loc[i,Alternate]]=site_data.loc[i,frequency]
			   
		total_alt = 0
	
		for j in versions:
		   total_alt += versions[j]
		
		versions[site_data.loc[i,Reference]] = 1-total_alt
		
		pi += site_pi(versions)
		
	if correction:
		pi = (pi*N_samples)/(N_samples-1)
		
	variance_pi = var_pi(pi,N_samples)
	
	
	if per_site:
		
		variance_pi = variance_pi/(sites[1]-sites[0]+1)
		pi = pi/(sites[1]-sites[0]+1)
	
	return (pi, variance_pi)
#%%
def output(output,outfilename,group):
	
	s=group+"\tpi\tvar\n"
	for key,value in output.items(): 
		s+=(str(key)+"\t"+str(value[0])+"\t"+str(value[1])+"\n")
		
	if len(outfilename) > 0:
		with open(outfilename,'a') as f_out:
			f_out.write(s)		
	else:
		print('No output file name provided. Writing to stdout:')
		print(s)	
		
	return
#%%
def main(argv):
	
	data, options = get_args(argv)
	
	pis={}
	#need some sanity checks for existence and values of args if len()
	for ind in data[options['group']].unique():
		pis[ind] = pi(pi_data=data.loc[data[options['group']]==ind],
		sites=(1, data[options['position']].max()),
		pos_id=options['position'],
		Reference=options['reference'],
		Alternate=options['alternate'],
		frequency=options['frequency'],
		correction=options['sizeCorrection'],
		N_samples=int(data[options['coverage']].mean()),
		per_site=options['perSite']
		)
		
	output(output=pis,outfilename=options['outputfile'],group=options['group'])
	
	print("\nDone!")
#%%

if __name__ == "__main__":
   main(sys.argv)
 #%%  