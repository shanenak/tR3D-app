
import plotly.graph_objects as go
import pandas as pd
import streamlit as st

def get_data():
    st.session_state.df = pd.read_csv(r"data/AlaX_cleaned_data.csv")
    st.session_state.attributes = pd.read_csv('./data/attributes.csv')
    

def make_bar(subset_df: pd.DataFrame()):
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=subset_df['Co-occurrence'],
            y=subset_df['Pfam Description'], 
            # hovertext=subset_df['# of Queries with Pfam Neighbors'],
            # hoverinfo="text",
            customdata=subset_df,
            hovertemplate='<b> %{y}</b><br># of queries with Pfam neighbors: %{customdata[1]: .2f}<extra></extra>',
            orientation='h',
            text = subset_df['Median Distance'],
            texttemplate='%{text:.0s} median distance', 
            textposition='outside',
            )
    )
    fig.update_layout(
        title="Co-occurrences for Cluster "+str(subset_df['SSN Cluster Number'].iloc[0]),
        xaxis = dict(
            title = 'Co-occurrence',
            tickformat = '.0%',
        ),
        yaxis = dict(
            title = 'Neighboring Protein Families',
            tickvals = subset_df['Pfam Description'],
            ticktext = subset_df['Pfam Description'].str[:30]
        ),
        # margin=dict(l=20, r=30, t=20, b=20),
        height = 400,
        # width = 1200,
    )

    fig.update_xaxes(range=[0, 1])
    return fig

def make_summary_bar(df: pd.DataFrame()):
    fig = go.Figure()
    
    pfam_df = df.loc[df['Pfam Description']!='none'].groupby(['Pfam Description'])['# of Queries with Pfam Neighbors'].count().sort_values(ascending=False).iloc[0:15].sort_values(ascending=True).reset_index()
    pfam_df['average_distance'] = pfam_df.apply(lambda x: df.loc[df['Pfam Description']==x['Pfam Description'], 'Median Distance'].median(),axis=1)

    fig.add_trace(
        go.Bar(
            y=pfam_df['Pfam Description'],
            x=pfam_df['# of Queries with Pfam Neighbors'], 
            # hovertext=subset_df['# of Queries with Pfam Neighbors'],
            # hoverinfo="text",
            # customdata=subset_df,
            # hovertemplate='<b> %{y}</b><br># of queries with Pfam neighbors: %{customdata[1]: .2f}<extra></extra>',
            orientation='h',
            text = pfam_df['average_distance'],
            texttemplate='%{text:.0s} median distance', 
            textposition='inside',
            )
    )
    fig.update_layout(
        title="Neighboring protein families by number of occurrences in dataset",
        yaxis = dict(
            title = 'Number of times referenced in data',
            tickvals = pfam_df['Pfam Description'],
            ticktext = pfam_df['Pfam Description'].str[:50]
        ),
        xaxis_title="Neighboring Protein Families",
        height = 500,
        width = 1200,
    )
    return fig