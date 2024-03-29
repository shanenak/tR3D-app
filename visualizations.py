
import plotly.graph_objects as go
import pandas as pd
import streamlit as st
    
# Main Repeated Graph on Clusters Page
def make_bar(subset_df: pd.DataFrame()):
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=subset_df['Co-occurrence'],
            y=subset_df['Pfam Description'], 
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
            titlefont=dict(
                size=18,
                color="black"
            ),
            tickformat = '.0%',
            tickfont = dict(
                color="black"
            )
        ),
        yaxis = dict(
            title = 'Neighboring Protein Families',
            titlefont=dict(
                size=18,
                color="black"
            ),
            tickvals = subset_df['Pfam Description'],
            ticktext = subset_df['Pfam Description'].str[:30],
            tickfont = dict(
                size=15,
                color="black"
            )
        ),
        height = 400,
        font=dict(
            size=18,
            color="black"
        )
    )

    fig.update_xaxes(range=[0, 1])
    return fig

#Summary figure at the top of the page
def make_summary_bar(df: pd.DataFrame()):
    fig = go.Figure()
    
    pfam_df = df.loc[df['Pfam Description']!='none'].groupby(['Pfam Description'])['# of Queries with Pfam Neighbors'].sum().sort_values(ascending=False).iloc[0:15].sort_values(ascending=True).reset_index()
    pfam_df['average_distance'] = pfam_df.apply(lambda x: df.loc[df['Pfam Description']==x['Pfam Description'], 'Median Distance'].median(),axis=1)

    fig.add_trace(
        go.Bar(
            y=pfam_df['Pfam Description'],
            x=pfam_df['# of Queries with Pfam Neighbors'], 
            hovertext=pfam_df['# of Queries with Pfam Neighbors'],
            hoverinfo="text",
            customdata=pfam_df,
            hovertemplate='<b> %{y}</b><br># of queries with Pfam neighbors: %{customdata[1]: .2f}<extra></extra>',
            orientation='h',
            text = pfam_df['average_distance'],
            texttemplate='%{text:.0s} median distance', 
            textposition='inside',
            )
    )

    fig.update_layout(
        title="Neighboring protein families by number of occurrences in dataset",
        xaxis = dict(
            title = 'Total Number of Occurrences',
            titlefont=dict(
                size=16,
                color="black"
            ),
            tickfont = dict(
                size=14,
                color="black"
            )
        ),
        yaxis = dict(
            title = 'Neighboring Protein Families',
            titlefont=dict(
                size=16,
                color="black"
            ),
            tickvals = pfam_df['Pfam Description'],
            ticktext = pfam_df['Pfam Description'].str[:50],
            tickfont = dict(
                size=14,
                color="black"
            )
        ),
        height = 500,
        width = 1000,
    )
    return fig