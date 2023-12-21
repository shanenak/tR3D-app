import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from streamlit_extras.switch_page_button import switch_page

from visualizations import make_bar, make_summary_bar

#### INSTRUCTIONS
#### 1) Select Terminal from the toolbar > New Terminal
#### 2) In your Terminal console, check if it says 'Power Shell' at the top right. 
#           If so, select the dropdown arrow and click Command Prompt.
#### 3) In your Terminal console, check if the last line starts with (venv). If so, your virtual environment is already activated. 
#           If not, copy the below code into the terminal and press Enter.
#           venv/Scripts/activate
#           venv\Scripts\activate actually
#### 4) Run the code by copying the below code into the terminal and pressing Enter. 
#           Streamlit run app.py
#### 5) Wait for the browser to open with the local version of your app. 
#### 6) Re-run the code by saving this file and selecting "Always re-run" in the browser. 

#this is v important - 
def AlaX_app():
    st.session_state.df = pd.read_csv(r"data/220902_IPR018163_Uniref50_43_480max_Cleaned_GNN_Final.csv") #was previously AlaX_10neighbors
    #reads file, needs to be input twice. But doesn't need to be cleaned?
    df = pd.read_csv(r"data/220902_IPR018163_Uniref50_43_480max_Cleaned_GNN_Final.csv") #was previously AlaX_10neighbors
    unique_pfam = df.groupby(by='Pfam Description')['Pfam Description'].count()
    # df['Pfam Description'] = df.apply(lambda x: x['Pfam Description']+'-'+str(x.index) if unique_pfam[x['Pfam Description']]>1 else x['Pfam Description'], axis=1)
    st.set_page_config(
        layout="wide",
        page_title = 'tR3D',
        page_icon = 'test tube',
        initial_sidebar_state='auto',
        menu_items={
            'Report a Bug':'mailto:dcmillar@berkeley.edu',
            'About': '# This app is being developed in the lab of Michelle Chang at UC Berkeley.'
        })
    with st.sidebar:
        st.subheader('tRNA-Deacylase Directed Discovery (tR3D) of Noncanonical Amino Acids and Enzymes')
        st.caption('AlaX, Interpro Family 18163')
        st.write('This app enables a enzyme family approach for Biosynhetic Gene Cluster (BGC) discovery through efficeint gene neighborhood analysis.')
        st.markdown("Contact [Douglas Millar PhD](mailto:dcmillar@berkeley.edu) to learn more.")
        
    # st.title('tRNA-Deacylase Directed Discovery (tR3D) of Noncanonical Amino Acids')
    # st.write('AlaX, Interpro Family 18163  |  Douglas Millar, Michelle Chang Lab, Chemical and Biomolecular Engineering, UC Berkeley')
    
    st.write('')
    st.header('Overview of all SSN Clusters')
    st.write('')
    st.write('This page provides a summary of all SSN Clusters. Review which protein families occur with the highest frequency. Use the filter to explore clusters with the identified protein families.')
    # st.markdown("[Identify High Frequency Protein Families](#high-frequency-protein-families)")
    # st.markdown("[Explore Co-occurrence in Clusters](#explore-co-occurrence-in-cluster)")

    
    st.write('')
    st.subheader('Identify High Frequency Protein Families')
    summary_fig = make_summary_bar(df)
    st.plotly_chart(summary_fig)
    
    #what is this line for? 
    st.write('')
    
    st.subheader('Explore Co-occurrence in Clusters')
    with st.expander(label='Filter clusters', expanded=False):
        col1, col2, col3 = st.columns([4, 1, 4])
        with col1:
            # st.write('Minimum Co-occurrence Value')
            cooccurrence_min = st.slider('Select minimum co-occurrence value to filter', min_value=0, max_value=100, value=75, format="%.0f%%")/100 # divide to make decimal
            cluster_cooccurrence_max = df.groupby(['SSN Cluster Number'])['Co-occurrence'].max()
            df = df.loc[df['SSN Cluster Number'].notna()]
            df['SSN Cluster Number'] = df['SSN Cluster Number'].astype(int)
            df = df.loc[df.apply(lambda x: cluster_cooccurrence_max[x['SSN Cluster Number']]>=cooccurrence_min, axis=1)]
            
            st.write('')
            # st.write('Maximum Neighboring Protein Families')
            selected_number_pfams = st.slider('Number of neighboring protein families to show', 0, 25, 15)
            
            # st.write('Subcluster Number Range')
            ssn_range = st.slider(
                    'Select SSN Cluster Numbers to filter', 
                    0, int(df['SSN Cluster Number'].max()),
                    (0, int(sorted(df['SSN Cluster Number'].unique())[19] if len(df['SSN Cluster Number'].unique())>20 else df['SSN Cluster Number'].unique().max())), 
                    # disabled = False if len(selected_shared_names)==0 else True
                    )
            df = df.loc[(df['SSN Cluster Number']>=ssn_range[0])&(df['SSN Cluster Number']<=ssn_range[1])]
            
        
        with col3:
            
            st.write('')
            # st.write('Shared Names')
            selected_shared_names = st.multiselect('Select shared names to filter', df['shared name'].unique())
            if len(selected_shared_names)>0:
                selected_sharedname_clusters = df.loc[df['shared name'].isin(selected_shared_names), 'SSN Cluster Number'].unique()
                df = df.loc[df['SSN Cluster Number'].isin(selected_sharedname_clusters)]
            
            st.write('')
            # st.write('Pfam Descriptions')
            selected_pfam = st.multiselect('Select protein families to filter', df['Pfam Description'].unique())
            if len(selected_pfam)>0:
                selected_pfam_clusters = df.loc[df['Pfam Description'].isin(selected_pfam), 'SSN Cluster Number'].unique()
                df = df.loc[df['SSN Cluster Number'].isin(selected_pfam_clusters)]
            

    st.title('')
    for SSN in sorted(df['SSN Cluster Number'].unique()):
        col1, col2, col3 = st.columns([4, 1, 3])
        with col1:
            st.subheader('Cluster '+str(SSN))
        
        subset_df = df.loc[(df['SSN Cluster Number']==SSN)&(df['shared name']!='none')].sort_values('Co-occurrence', ascending=False)
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
            
        # if st.button('See Cluster ' + str(SSN) + ' in more detail', type='primary'):
        #         switch_page("Clusters")
        st.write('')
        st.write('')
        


if __name__ == "__main__":
    AlaX_app()