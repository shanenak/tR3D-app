import pandas as pd
import streamlit as st
from visualizations import make_bar, make_summary_bar
    
# st.title('tRNA-Deacylase Directed Discovery (tR3D) of Noncanonical Amino Acids')
# st.write('AlaX, Interpro Family 18163  |  Douglas Millar, Michelle Chang Lab, Chemical and Biomolecular Engineering, UC Berkeley')

st.write('')
st.header('Evaluate Single Subcluster')
st.write('')
st.write('This page provides a deep-dive into a single SSN subcluster. Review which protein families occur with the highest frequency. Use the filter to explore subclusters with the identified protein families.')

subcluster_df = pd.read_csv('./data/subclusters.csv')
subcluster_df['SSN Cluster Number'] = subcluster_df['SSN Cluster Number'].astype(int)
subcluster_df['Subcluster Number'] = subcluster_df['Subcluster Number'].astype(int)

SSN = st.selectbox('Select cluster to evaluate', sorted([int(x) for x in subcluster_df['SSN Cluster Number'].unique() if x>0]), index=3)
df = subcluster_df.loc[subcluster_df['SSN Cluster Number']==SSN]
st.write('')

st.write('')
st.subheader('Identify High Frequency Protein Families')
summary_fig = make_summary_bar(df)
st.plotly_chart(summary_fig)

st.write('')

st.subheader('Explore Co-occurrence in Subclusters')
with st.expander(label='Filter subclusters', expanded=False):
    col1, col2, col3 = st.columns([4, 1, 4])
    with col1:
        # st.write('Minimum Co-occurrence Value')
        cooccurrence_min = st.slider('Select minimum co-occurrence value to filter', min_value=0, max_value=100, value=75, format="%.0f%%")/100 # divide to make decimal
        subclusters_cooccurrence_max = df.groupby(['Subcluster Number'])['Co-occurrence'].max()
        df = df.loc[df.apply(lambda x: subclusters_cooccurrence_max[x['Subcluster Number']]>=cooccurrence_min, axis=1)]
        
        st.write('')
        # st.write('Maximum Neighboring Protein Families')
        selected_number_pfams = st.slider('Number of neighboring protein families to show', 0, 25, 15)
        
        # st.write('Subcluster Number Range')
        ssn_range = st.slider(
                'Select Subcluster Numbers to filter', 
                0, int(df['Subcluster Number'].max()),
                (0, int(sorted(df['Subcluster Number'].unique())[19] if len(df['Subcluster Number'].unique())>20 else df['Subcluster Number'].unique().max())), 
                # disabled = False if len(selected_shared_names)==0 else True
                )
        df = df.loc[(df['Subcluster Number']>=ssn_range[0])&(df['Subcluster Number']<=ssn_range[1])]
    
    with col3:
        
        st.write('')
        # st.write('Shared Names')
        selected_shared_names = st.multiselect('Select shared names to filter', df['shared name'].unique())
        if len(selected_shared_names)>0:
            selected_subclusters = df.loc[df['shared name'].isin(selected_shared_names), 'Subcluster Number'].unique()
            df = df.loc[df['Subcluster Number'].isin(selected_subclusters)]
            
        st.write('')
        # st.write('Pfam Descriptions')
        selected_pfam = st.multiselect('Select protein families to filter', df['Pfam Description'].unique())
        if len(selected_pfam)>0:
            selected_pfam_clusters = df.loc[df['Pfam Description'].isin(selected_pfam), 'SSN Cluster Number'].unique()
            df = df.loc[df['SSN Cluster Number'].isin(selected_pfam_clusters)]
        

st.title('')
for SSN in sorted(df['Subcluster Number'].unique()):
    col1, col2, col3 = st.columns([4, 1, 3])
    with col1:
        st.subheader('Subcluster '+str(SSN))
    
    subset_df = df.loc[(df['Subcluster Number']==SSN)&(df['shared name']!='none')].sort_values('Co-occurrence', ascending=False)
    col1, col2= st.columns([1,4])
    with col1:
        st.title('')
        st.metric(label = 'Sequences in SSN Cluster', value = subset_df.iloc[0]['# of Sequences in SSN Cluster'])
        st.write('')
        st.metric(label = 'Max co-occurrence', value = str(int(subset_df['Co-occurrence'].max()*100))+'%')
        st.write('')
        st.metric(label = 'Unique Protein Families', value = len(subset_df.index.unique()))
    with col2:
        if len(subset_df['shared name'].unique())>=15:
            top_shared_names = subset_df['shared name'].unique()[:selected_number_pfams]
            subset_df = subset_df.loc[subset_df['shared name'].isin(top_shared_names)]
        temp_fig = make_bar(subset_df.sort_values('Co-occurrence'))
        st.plotly_chart(temp_fig)
        
    st.write('')
    st.write('')