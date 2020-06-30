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
"""explaination of variables in the dictionary "options"

    -p, --position <column name>   column name containing position information (default: -p pos)
    -g, --grouping <column name>   column name containing group identifiers (default: -g index)
    -r, --reference <column name>  column name of reference allele (default: -r ref)
    -a, --alternate <column name>  column name of alternate allele (default: -a alt)
    -v, --sampling_variance        flag to use sampling variance instead of default variance; REQUIRES COVERAGE
    -c, --coverage <column name>   column name of total read depth or False for no coverage data (default: -c False)
    -f, --frequency <column name>  column name of SNP frequency (default: -f freqProp)
    -l, --length <integer>         length of region to calculate pi, default: maximum value in positions
    -s, --per-site                 if used, calculate per-site pi instead of region pi (default)
    -z, --sample-size-correction   if used, calculate sample size correction for pi (default: no correction)
    -y, --covarying-sites          if used, calculate the upper-bound variance of pi
    -d, --ld-correction            if used, adjust covariance by LD decay
    """
import pandas as pd
import numpy as np
#%%
def get_temp_args():
    #standard options for PIRE temperature experiment
    options = {
            'group':'index',
            'reference':'ref',
            'position':'pos',
            'alternate':'alt',
            'coverage': False,
            'frequency':'freqPropMeanNoNA',
            'perSite':True,
            'sizeCorrection':False,
            'length':2800,
            'sampling_variance':False,
            }

    return options
#%%
def get_veg_prop_args():
    #standard options for PIRE vegetative propogation experiment
    options = {
            'group':'index',
            'reference':'ref',
            'position':'pos',
            'alternate':'alt',
            'coverage':False,
            'frequency':'freqPropMeanNoNA',
            'perSite':False,
            'sizeCorrection':False,
            'length':2800,
            'sampling_variance':False,
            }

    return options

#%% variances

def sampling_var(pi,n):
    #variance of pi from equation 27 in Fu and Li 1992
    varPi = (pi*((n+1)/(3*n-3)))+((pi**2)*((2*(n**2+n+3))/(9*(n**2-n))))

    return varPi

def var_pi(p,q):
    #calculate site variance as product of two binomial random variables
    var = ((p*(1-p))*(q*(1-q))) + ((p*(1-p))*(q**2)) + ((q*(1-q))*(p**2))

    return var

def cov_pi(a):
    #calculate variance of a region assuming maximum covariance between sites
    z = a.shape[0]

    for x in range(len(z)):
        for y in range(len(z)):
            if x != y:
                a[x,y] = np.sqrt(a[x,x]*a[y,y])

    return a
"""
#not implemented yet
def ld_adjusted_cov(a,length,function):
    #a is the variance-covariance matrix
    x,y = np.where(a!=0)

    for ind in x.shape[0]:
        les = min(x[ind],y[ind])
        mor = max(x[ind],y[ind])

        ld_dist = min((mor-les),(length+les-mor))

        a[x[ind],y[ind]] = a[x[ind],y[ind]] * function(ld_dist)

    return a
"""
#%% Site_pi
def site_pi(versions):
    #caculate pi from a single site; DO NOT multiply by 2!
    #takes in a dictionary with A,G,C,and/or T as keys and frequencies (0.0-1.0)
        #as values

    pi = 0
    prod_var = 0
    for i in versions:
        for j in versions:
            if i != j:
                pi += versions[i]*versions[j]
                prod_var += var_pi(versions[i],versions[j])
    return pi,prod_var
#%% main pi function
def pi(pi_data,
       sites,
       N_samples,
       pos_id,
       Reference,
       Alternate,
       frequency,
       correction,
       per_site,
       sampling_variance,
       lngth):

    #calculate pi per site from snp frequencies
    #read in a pandas DataFrame object of major or minor allele frequencies,
    #one frequency per row
    #assumes SNPs only, no indels;does NOT assume biallelic sites
    #calculate pi for the region, NOT per-site pi
    #optional sample size correction

    #subset_data = pi_data.loc[pi_data[pos_id]>=sites[0]].loc[pi_data[pos_id]<sites[1]]


    region_var = np.zeros((int(lngth+1),int(lngth+1)))
    pi_a = 0

    #extract rows from the ps.DataFrame "data" for this site
    for site in pi_data[pos_id].unique():

        site_data = pi_data.loc[pi_data[pos_id]==site]
        #set dictionary of versions and frequency of each: alternates declared,
        #    reference calculated
        #check to see if there are more than two alleles
        if site_data.shape[0] > 1:
            #recalculate frequencies
            dp = site_data.DP.sum()
            site_data[frequency] = [x/dp for x in site_data['AD']]


        versions={}
        for i in site_data[Alternate].unique():
            try:
                versions[i] = site_data.loc[site_data[Alternate] == i,frequency][0]
            except:
                versions[i] = site_data.loc[site_data[Alternate] == i,frequency].array[0]

        total_alt = 0

        for j in versions:
           total_alt += versions[j]

        try:
            versions[site_data[Reference][0]] = 1-total_alt
        except:
            versions[site_data[Reference].array[0]] = 1-total_alt

        site_p, site_var = site_pi(versions)

        region_var[int(site),int(site)] = site_var
        pi_a += site_p

    if correction:
        if N_samples:
            pi_a = (pi_a*N_samples)/(N_samples-1)

    if sampling_variance or N_samples:

        if sampling_variance and N_samples:
            variance_pi = sampling_var(pi_a,N_samples)
        else:
            variance_pi = region_var.trace()

    else:
        variance_pi = region_var.trace()

    if per_site:
        if lngth > 0:
            pi_a = pi_a/lngth
            variance_pi = variance_pi/lngth

        else:
            variance_pi = variance_pi/(sites[1]-sites[0]+1)
            pi_a = pi_a/(sites[1]-sites[0]+1)

    return (pi_a, variance_pi)
#%%
def aggregate_groups(data,columns,DP='DP',AD='AD',freq='freqProp'):

    aggregated_data = data.groupby(by=columns).agg(sum)
    aggregated_data[freq] = aggregated_data[AD] / aggregated_data[DP]



    return aggregated_data

def average_groups(data,columns):

    #this function groups the data by `columns` and aggregates other columns by
    #    returning the mean of the rows within each group
    averaged_data = data.groupby(by=columns).agg(np.mean)

    return averaged_data

#%% main
#set up local functions to group our data sets and call the loaded pi_for_temperature_notebook.py functions
def build_pi_df(pis,groupby):

    pis_dict = {}

    for group in groupby:

        pis_dict[group] = []

    pis_dict['pi'] = []

    for key in pis.keys():
        pis_dict['pi'].append(pis[key][0])
        for i in range(len(groupby)):
            pis_dict[groupby[i]].append(key[i])

    df = pd.DataFrame.from_dict(pis_dict, orient='columns')

    return df


def get_group_pis(data,options,group_by):

    pis_func = {}

    grouped_data = data.groupby(by=group_by)

    for group in grouped_data:

        try:
            N_samples = int(group[1][options['coverage']].mean())
        except:
            N_samples = False

        pis_func[group[0]] = pi(pi_data=group[1],
            sites=(1, group[1][options['position']].max()),
            N_samples=N_samples,
            pos_id=options['position'],
            Reference=options['reference'],
            Alternate=options['alternate'],
            frequency=options['frequency'],
            correction=options['sizeCorrection'],
            per_site=options['perSite'],
            sampling_variance=options['sampling_variance'],
            lngth=options['length'])

    return build_pi_df(pis_func,group_by)
