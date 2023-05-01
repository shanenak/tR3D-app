import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import numpy as np
import os

# from visualizations import make_bar, make_summary_bar
    
# st.title('tRNA-Deacylase Directed Discovery (tR3D) of Noncanonical Amino Acids')
# st.write('AlaX, Interpro Family 18163  |  Douglas Millar, Michelle Chang Lab, Chemical and Biomolecular Engineering, UC Berkeley')

st.write('')
st.header('Evaluate Single Cluster or Subcluster')
st.write('')
st.write('This page provides a deep-dive into a single SSN subcluster. Review which protein families occur with the highest frequency. Use the filter to explore subclusters with the identified protein families.')

# neighbors = pd.read_csv('./data/neighbors.csv')
dir = r"data\20230416\neighbors"
neighbors = pd.DataFrame()
for filename in os.listdir(dir):
    f = os.path.join(dir, filename)
    # checking if it is a file
    if os.path.isfile(f):
        temp_df = pd.read_csv(f)
        temp_df.rename(columns={'cluster_num':'Subcluster Number'}, inplace=True)
        cluster_number_long = filename.split("_")[1].split('.')[0]
        cluster_number = cluster_number_long.replace("cluster", "")
        temp_df['SSN Cluster Number'] = cluster_number
        neighbors = pd.concat([neighbors, temp_df], axis=0)
        
neighbors['transporter'] = neighbors['desc'].apply(lambda x: 'transport' in x if type(x)==str else False)
neighbors['regulator'] = neighbors['desc'].apply(lambda x: 'regulat' in x if type(x)==str else False)
neighbors['EamA'] = neighbors['desc'].apply(lambda x: 'EamA' in x if type(x)==str else False)

attributes = pd.read_csv('./data/attributes.csv')
attributes.rename(columns={'sort_key':'gene_key'},inplace=True)
subcluster_df = pd.read_csv('./data/AlaX_10neighbors.csv')

SSN = st.selectbox(label = 'Select cluster', 
                options = ['Any']+sorted(attributes['SSN Cluster Number'].unique()), 
                index=5
                )
if SSN!='Any':
    # neighbor_df = neighbor_df.loc[neighbor_df['SSN Cluster Number']==SSN]
    attributes = attributes.loc[attributes['SSN Cluster Number']==SSN]
    flag_proteins = attributes.loc[attributes['SSN Cluster Number']==SSN].index
    subcluster_df = subcluster_df.loc[subcluster_df['SSN Cluster Number']==SSN]

subcluster = st.selectbox(label = 'Select subcluster', 
                options = ['Any']+sorted(attributes['Subcluster Number'].unique()),
                index=7
                )
if subcluster!='Any':
    attributes = attributes.loc[attributes['Subcluster Number']==subcluster]
    flag_proteins = attributes.loc[attributes['Subcluster Number']==subcluster].index

st.write('')
st.subheader('Identify High Frequency Protein Families')

THRESHOLD = 10
selected_number_pfams = 5

subset_df = subcluster_df.sort_values('Co-occurrence', ascending=False)
top_shared_names = [x for x in subset_df['shared name'].unique() if x!='none'][:selected_number_pfams]

# flag_proteins = attributes.loc[attributes['SSN Cluster Number']==SSN].index
# subcluster_df = subcluster_df.loc[subcluster_df['SSN Cluster Number']==SSN]

count_neighbors = {}
dict_flag_proteins = {}
filtered_neighbors = pd.DataFrame()
for flag in flag_proteins:
    taxon = attributes.loc[flag, 'taxon_id']
    flag_num = attributes.loc[flag, 'num']
    gene_key = attributes.loc[flag, 'gene_key']
    flag_neighbors = neighbors.loc[(neighbors['taxon_id']==taxon)&(neighbors['gene_key']==gene_key)&(neighbors['num']>=(flag_num-THRESHOLD))&(neighbors['num']<=(flag_num+THRESHOLD))&(neighbors['family']!='none')]
    # flag_neighbors['taxon_id'] = taxon
    flag_neighbors['num'] = flag_num
    # flag_neighbors['gene_key'] = gene_key
    dict_flag_proteins = flag_neighbors['family_desc'].unique()
    filtered_neighbors = pd.concat([filtered_neighbors, flag_neighbors], axis=0)
    for neighbor in flag_neighbors['family_desc'].unique():
        count_neighbors[neighbor] = count_neighbors[neighbor]+1 if neighbor in count_neighbors.keys() else 1
neighbor_list = filtered_neighbors.groupby(['taxon_id', 'num', 'gene_key'])['family_desc'].apply(list).reset_index(name="family_desc")

NUM_PRIMARY = 3
top_neighbors = dict(sorted(count_neighbors.items(), key=lambda x:x[1], reverse=True)[:NUM_PRIMARY])
NUM_SECONDARY = 5
count_secondaryneighbors = pd.DataFrame(columns=['index'])
for neighbor in top_neighbors.keys():
    neighbor_list[neighbor+'_presence'] = neighbor_list.apply(lambda x: neighbor in x['family_desc'], axis=1)
    temp = neighbor_list.loc[neighbor_list[neighbor+'_presence']].explode('family_desc')
    temp_secondary = temp['family_desc'].value_counts(ascending=False).reset_index(name=neighbor).iloc[:NUM_SECONDARY,:]
    temp_secondary.rename(columns={'family_desc':'index'},inplace=True)
    # count_secondaryneighbors = pd.concat([count_secondaryneighbors, temp_secondary[neighbor]], axis=1)
    count_secondaryneighbors = pd.merge(count_secondaryneighbors, temp_secondary, how="outer", on='index')
    # count_secondaryneighbors = count_secondaryneighbors.join(temp_secondary)
count_secondaryneighbors = count_secondaryneighbors.rename(columns={'index':'pair'}).set_index('pair')
neighbor_list.reset_index(drop=True, inplace=True)

st.dataframe(count_secondaryneighbors, use_container_width=True)
st.write('')

st.subheader('Select pair from high frequency protein family table')
primary = st.selectbox('Select primary protein to evaluate', top_neighbors.keys())
secondary = st.selectbox('Select secondary protein to evaluate', [sec_neighbor for sec_neighbor in count_secondaryneighbors.sort_values(by=primary, ascending=False).index[:NUM_SECONDARY] if sec_neighbor!=primary])

neighbor_list[secondary+'_presence'] = neighbor_list.apply(lambda x: secondary in x['family_desc'], axis=1)

pair_df = neighbor_list.loc[(neighbor_list[primary+'_presence']==True)&(neighbor_list[secondary+'_presence']==True)]

st.write('')
st.write('')
for sequence in pair_df.index:
    taxon = pair_df.loc[sequence, 'taxon_id']
    flag_num = pair_df.loc[sequence, 'num']
    gene_key = pair_df.loc[sequence, 'gene_key']
    flag_protein = attributes.loc[(attributes['taxon_id']==taxon)&(attributes['num']==flag_num)]
    st.subheader(flag_protein['accession'].iloc[0])
    st.write(flag_protein['desc'].iloc[0]+' in '+flag_protein['organism'].iloc[0])
    
    flag_protein['selected'] = '' # primary
    flag_protein['family_desc'] = 'AlaX'
    flag_protein['type'] = 'AlaX'
     
    temp_neighbors = neighbors.loc[(neighbors['taxon_id']==taxon)&(neighbors['gene_key']==gene_key)&(neighbors['num']>=(flag_num-THRESHOLD))&(neighbors['num']<=(flag_num+THRESHOLD))&(neighbors['family']!='none')]
    temp_neighbors['selected'] = temp_neighbors.apply(lambda x: '' if (secondary in x['family_desc'])|(primary in x['family_desc']) else '-open', axis=1) # empty if secondary, else -open for marker
    temp_neighbors['type'] = temp_neighbors['desc'].apply(lambda x: 'Transporter' if 'transport' in str(x) else('Regulator' if 'regulat' in str(x) else 'Other'))

    proteins = pd.concat([flag_protein,temp_neighbors], axis=0)
    proteins['direction'] = proteins['direction'].apply(lambda x: 'up' if x=='normal' else 'down')
    proteins.sort_values('num', inplace=True)
    proteins.reset_index(drop=True, inplace=True)
        
    # TODO: check why some EamA aren't categorized as transporter
    
    # TODO: sort by organism name
    
    # TODO: extract chunk and find frequency of pathways
    
    # TODO: color primary and secondary throughout the app
 
    st.write('')
    col1, col2 = st.columns([1,3])
    with col1:
        primary_neighbor = proteins.loc[proteins['family_desc']==primary]
        primary_ave_distance = primary_neighbor['num'].apply(lambda x: abs(x-flag_num)).min()
        st.metric(primary + ' distance', primary_ave_distance)
        secondary_neighbor = proteins.loc[proteins['family_desc']==secondary]
        secondary_ave_distance = secondary_neighbor['num'].apply(lambda x: abs(x-flag_num)).min()
        st.metric(secondary + ' distance', secondary_ave_distance)
    with col2:
        fig = go.Figure()
        type_colors = {'AlaX':'DarkBlue', 'Transporter': 'DarkGreen', 'Regulator':'DarkRed',  'Other':'DarkOrange'}
        for type in type_colors.keys():
            temp_proteins = proteins.loc[proteins['type']==type]
            if len(temp_proteins.index)>0:
                fig.add_trace(go.Scatter(
                    x=[0]*len(temp_proteins.index),
                    y=temp_proteins.index, #proteins.apply(lambda x: x[['start', 'stop']].mean(),axis=1),
                    mode="markers+text",
                    name=type,
                    marker_size = temp_proteins['selected'].apply(lambda x: 10 if x!='' else 30),
                    marker_color = type_colors[type],
                    marker_line_width=3,
                    marker_symbol = temp_proteins.apply(lambda x: 'triangle-'+x['direction']+x['selected'], axis=1),
                    text=temp_proteins['family_desc'], #.apply(lambda x: (str(x['family_desc'])+', type:'+str(x['type'])) if x['type']!=None else x['family_desc'], axis=1),
                    textposition="middle right"
                ))
        fig.update_xaxes(showgrid=False, zeroline=True, showticklabels=False)
        fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False)
        fig.update_layout(xaxis_range=[-1,5])
        fig.update_layout(
            margin=dict(l=0, r=0, t=0, b=0),
            height = 500,
            width = 600,
        )
        
        st.plotly_chart(fig)

