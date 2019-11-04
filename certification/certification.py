#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 14:20:59 2019

@author: alex
"""

"""
Code for analysing recording quality in certification recordings
Current code covers:
Unit yield by channel
Average amplitude  by channel
Requirements:

"""
## Import basic functions

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import csv
import os
import matplotlib.patches as mpatches
from os import path


## Functions

def probe_merger(folder_path):
    '''
    Merges probe data for labs that did not run merger software
    INPUTs: folder with probesephys data
    OUTPUT: merged npy data
    '''
    files  = sorted(os.listdir(folder_path))
    if '.DS_Store' in files:
            files.remove('.DS_Store')  # For mac computers
    a  = np.load(folder_path + '/'+ files[0])
    b  = np.load(folder_path + '/'+ files[1])
    merge  = np.concatenate((a,b))
    np.save(folder_path +'/'+files[0], merge )

def cluster_amps_import(tsv_file):
    '''
    Importing tsv file with average amplitude for each cluster
    INPUTs: tsv files
    OUTPUT: dataframe with cluster amplitudes
    '''
    with open(tsv_file) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        next(rd, None) 
        clusters_amps = pd.DataFrame(columns = ['cluster','amp'] )
        for i, row in enumerate (rd):
            row = {'cluster' : row[0], 'amp' : row[1]}
            clusters_amps = clusters_amps.append(row, ignore_index =  True)
    return clusters_amps

def cluster_amps_from_spike_amps(spike_amps,spike_clusters, cluster_depths):
    '''
    cluster_amps_from_spike_amps
    INPUTs: spike.amps,spike.clusters (alf files as npy)
    OUTPUT: cluster.amps npy
    '''
    
    
    cluster_amp_spike  =  pd.DataFrame({'Cluster': spike_clusters,\
                            'amplitude': spike_amps})
    cluster_map = \
    pd.DataFrame({'amplitude':np.zeros(len(cluster_depths))})
    cluster_amps_raw = cluster_amp_spike.groupby(['Cluster']).mean()
    
    cluster_amps\
    = cluster_map.join(cluster_amps_raw, lsuffix='_caller', rsuffix='')
    cluster_amps = \
            np.array(cluster_amps.sum(axis=1, skipna=True))
            
    return cluster_amps




## Import unit data:

# Import files (ONE in the future):
certification_folder = '/Users/alex/Desktop/Certification/'

labs = os.listdir(certification_folder)
labs.remove('.DS_Store')  # For mac computers



yield_df =  pd.DataFrame(columns =  ['Lab', 'Mouse','Day','Cluster_number', \
                  'Amplitude', \
                  'Depth' , \
                  'Binned_depth', 'total_units','units_per_probe', 'firing_rate'])

for lab in labs:
    mice = os.listdir(certification_folder + '/' + lab + '/')
    if '.DS_Store' in mice:
            mice.remove('.DS_Store')  # For mac computers
    for mouse in mice:
        days = \
        sorted(os.listdir(certification_folder + '/' + lab + '/'\
                          + mouse + '/'))
        if '.DS_Store' in days:
            days.remove('.DS_Store')  # For mac computers
        for d, day in enumerate(days):
            # Load files
            alf = os.listdir(certification_folder + '/' + lab + '/' +  \
                             mouse + '/' +  day + '/')
            amp_file = \
            [files for files in alf if files.startswith('clusters.amps')][0]
            depth_file = \
            [files for files in alf if files.startswith('clusters.depths')][0]
            if amp_file.endswith('tsv'):
                cluster_amps = \
                cluster_amps_import(certification_folder+lab +'/'+amp_file)
            else:
                cluster_amps = \
                np.load(certification_folder + lab + '/' +  \
                             mouse + '/' +  day + '/'+ amp_file)
            cluster_depths = \
            np.load(certification_folder + lab + '/' +  \
                             mouse + '/' +  day + '/' + depth_file)
            
            # Bin data by depth
            bins = np.arange(0,4000,100)
            #Reverse since 0 is bottom of the probe (Not very intuitive)
            binned_depth =  4000 - (np.digitize(cluster_depths, bins) * 100)
            
            #Calculate firing rates
            spikes_clusters_file = \
            [files for files in alf if files.startswith('spikes.clusters')][0]
            spikes_times_file = \
            [files for files in alf if files.startswith('spikes.times')][0]
            
            spikes_clusters = \
                np.load(certification_folder + lab + '/' +  \
                             mouse + '/' +  day + '/'+ spikes_clusters_file)
            spikes_times = \
                np.load(certification_folder + lab + '/' +  \
                             mouse + '/' +  day + '/'+ spikes_times_file)
            
            #the results of this is the firing rate per cluster. Where cluster 
            #did not fire a 0 is assinged
            cluster_rate = \
            pd.DataFrame({'firing_rate':np.zeros(len(cluster_depths))})
            cluster_times  =  pd.DataFrame({'Cluster': spikes_clusters,\
                            'firing_rate': spikes_times})
            cluster_firing = cluster_times.groupby(['Cluster']).count() \
            / max(spikes_times)
            cluster_firing = \
            cluster_rate.join(cluster_firing, lsuffix='_caller', rsuffix='')
            cluster_firing = \
            np.array(cluster_firing.sum(axis=1, skipna=True))
    
    
            
            # Add to dataframe and normalize for 2 probes
            if path.exists(certification_folder + '/' + lab + '/' +  \
                             mouse + '/' +  day + '/' + '2_probes.flag.rtf'):
                
                clusters_lab  = pd.DataFrame({'Lab': lab, \
                                              'Mouse': mouse, \
                                              'Day': d+1, \
                                              'Cluster_number': \
                                              range(len(cluster_amps)), \
                                              'Amplitude' : cluster_amps, \
                                              'Depth' : cluster_depths, \
                                              'Binned_depth' : binned_depth,\
                                              'total_units': len(cluster_amps),
                                              'units_per_probe':\
                                              len(cluster_amps)/2,\
                                              'firing_rate': cluster_firing })
    
            else: 
                clusters_lab  = pd.DataFrame({'Lab': lab, \
                                              'Mouse': mouse, \
                                              'Day': d+1, \
                                              'Cluster_number': \
                                              range(len(cluster_amps)), \
                                              'Amplitude' : cluster_amps, \
                                              'Depth' : cluster_depths, \
                                              'Binned_depth' : binned_depth,\
                                              'total_units': len(cluster_amps),\
                                              'units_per_probe':\
                                              len(cluster_amps),
                                              'firing_rate': cluster_firing })
            
            yield_df = yield_df.append(clusters_lab, ignore_index = True)

yield_df['total_units'] = yield_df['total_units'].astype(int)
yield_df['units_per_probe'] = yield_df['units_per_probe'].astype(float)

## Yield across days and channels

# Average unit per channel pooled the two probes  per lab across days
figure,ax = plt.subplots(1,1, figsize=(10,10))
plt.sca(ax)
sns.lineplot(x = 'Day',  y = 'units_per_probe' ,\
                     hue = 'Lab', data = yield_df)
ax.set_title('Number of units over days')
ax.set_ylabel('Number of units')
ax.set_xlabel('Recording day')
ax.set_ylim([0, 1000])
figure.savefig('units_per_day.pdf')
        

# Average number of units pooled the two probes  per lab across days
figure,ax = plt.subplots(7,3, figsize=(24,100))
for l, lab in enumerate(labs):
    mice = os.listdir(certification_folder + '/' + lab + '/')
    mice.remove('.DS_Store')  # For mac computers
    for m, mouse in enumerate(mice):
        plt.sca(ax[l,m])
        units = yield_df.loc[yield_df['Mouse'] == \
                mouse].groupby(['Binned_depth','Day']).count()\
                             ['units_per_probe'].reset_index()
        units = units.pivot('Binned_depth','Day','units_per_probe')
        sns.heatmap(units, annot=True, cbar=False, vmin=0, vmax=100)
        ax[l,m].set_title(lab +' '+ mouse + ' ' + 'Number of units')
        ax[l,m].set_ylabel('Binned Depth (100 µm bin)')
        ax[l,m].set_xlabel('Recording day')

figure.savefig('units_per_depth.pdf')


# Average amplitude pooled the two probes  per lab across days
figure,ax = plt.subplots(7,3, figsize=(24,100))
for l, lab in enumerate(labs):
    mice = os.listdir(certification_folder + '/' + lab + '/')
    mice.remove('.DS_Store')  # For mac computers
    for m, mouse in enumerate(mice):
        plt.sca(ax[l,m])
        units = yield_df.loc[yield_df['Mouse'] == \
                mouse].groupby(['Binned_depth','Day']).mean()\
                             ['Amplitude'].reset_index()
        units = units.pivot('Binned_depth','Day','Amplitude')
        sns.heatmap(units, annot=True, cbar=False, fmt='.2f', vmin=0, vmax=60)
        ax[l,m].set_title(lab +' '+ mouse + ' ' + 'Amplitude (µV)')
        ax[l,m].set_ylabel('Binned Depth (100 µm bin)')
        ax[l,m].set_xlabel('Recording day')

figure.savefig('amplitude_per_depth.pdf')

# Fraction of units by amplitude  per lab
figure,ax = plt.subplots(1,1, figsize=(10,10))
label_patches = []
colors =  ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
for i, lab in enumerate(labs):
    sns.kdeplot(yield_df.loc[yield_df['Lab']==lab, 'Amplitude'], shade=True, \
                color = colors[i], legend = False)
    label = lab
    label_patch = mpatches.Patch(
        color=colors[i],
        label=label)
    label_patches.append(label_patch)
plt.legend(handles=label_patches, loc='upper left')
ax.set_title('Amplitude distributions')
ax.set_ylabel('Kernel Density Estimate')
ax.set_xlabel('Amplitude(µV)')
ax.set_xlim([0,60])
            
figure.savefig('amplitude_per_lab.pdf')

# Fraction of units by firing rate  per lab
figure,ax = plt.subplots(1,1, figsize=(10,10))
label_patches = []
colors =  ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
for i, lab in enumerate(labs):
    sns.kdeplot(yield_df.loc[yield_df['Lab']==lab, 'firing_rate'], shade=True, \
                color = colors[i], legend = False)
    label = lab
    label_patch = mpatches.Patch(
        color=colors[i],
        label=label)
    label_patches.append(label_patch)
plt.legend(handles=label_patches, loc='upper left')
ax.set_title('Firing rate distributions')
ax.set_ylabel('Kernel Density Estimate')
ax.set_xlabel('Firing Rate(Hz)')   
ax.set_xlim([0,60])

figure.savefig('firing_rate_per_lab.pdf')


#Plot IBL averages

#Units vs Depth
figure,ax = plt.subplots(1,1, figsize=(10,10))
plt.sca(ax)
units = yield_df.groupby(['Binned_depth','Day']).count()\
                             ['total_units'].reset_index()
units = units.pivot('Binned_depth','Day','total_units')
sns.heatmap(units, annot=True, cbar=False)
ax.set_title( 'Number of units')
ax.set_ylabel('Binned Depth (100 µm bin)')
ax.set_xlabel('Recording day')

#Amplitude
figure,ax = plt.subplots(1,1, figsize=(10,10))
sns.kdeplot(yield_df['Amplitude'], shade=True, legend = False)
ax.set_title('IBL amplitude')
ax.set_ylabel('Kernel Density Estimate')
ax.set_xlabel('Amplitude(µV)')       

#Firing rate
figure,ax = plt.subplots(1,1, figsize=(10,10))
sns.kdeplot(yield_df['firing_rate'], shade=True, legend = False)
ax.set_title('IBL firing rate')
ax.set_ylabel('Kernel Density Estimate')
ax.set_xlabel('Firing Rate(Hz)')
ax.set_xlim([0,60])