import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import streamlit as st

# retrieve PPI BioGRID
def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "e5b612d404ccf5884acc4f5230999f64",  
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,
        "searchbiogridids": True,
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    network = response.json()
    network_df = pd.DataFrame.from_dict(network, orient='index')
    return network_df

#  from STRING
def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606
    }
    response = requests.get(string_url, params=params)
    network = response.json()
    network_df = pd.json_normalize(network)
    return network_df

# generate network from DataFrame
def generate_network(dataframe):
    if "OFFICIAL_SYMBOL_A" in dataframe.columns and "OFFICIAL_SYMBOL_B" in dataframe.columns:
        network_graph = nx.from_pandas_edgelist(dataframe, "OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B")
    else:
        network_graph = nx.from_pandas_edgelist(dataframe, "preferredName_A", "preferredName_B")
    return network_graph

# calculate centrality 
def get_centralities(network_graph):
    degree_centrality = nx.degree_centrality(network_graph)
    betweenness_centrality = nx.betweenness_centrality(network_graph)
    closeness_centrality = nx.closeness_centrality(network_graph)
    eigenvector_centrality = nx.eigenvector_centrality(network_graph)
    pagerank_centrality = nx.pagerank(network_graph)
    
    return {
        "Degree Centrality": degree_centrality,
        "Betweenness Centrality": betweenness_centrality,
        "Closeness Centrality": closeness_centrality,
        "Eigenvector Centrality": eigenvector_centrality,
        "PageRank Centrality": pagerank_centrality
    }


st.title("Protein-Protein Interaction (PPI) Network AZHEEM")

# Input 
protein_id = st.text_input("Enter Protein ID:")
database_choice = st.selectbox("Select Database", ("BioGRID", "STRING"))

if st.button("Retrieve PPI Data"):
    # Retrieve PPI data based on user choice
    if database_choice == "BioGRID":
        ppi_data = retrieve_ppi_biogrid(protein_id)
    else:
        ppi_data = retrieve_ppi_string(protein_id)

    # Check if data is retrieved successfully
    if not ppi_data.empty:
        # Generate network graph
        network_graph = generate_network(ppi_data)

        # Split layout into two columns
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("PPI Data Information")
            st.dataframe(ppi_data)  # Use st.dataframe for a more interactive table
            st.write("Number of edges:", network_graph.number_of_edges())
            st.write("Number of nodes:", network_graph.number_of_nodes())
            
            # Visualize network
            plt.figure(figsize=(10, 10))
            slayout = nx.spring_layout(network_graph, seed=123)
            nx.draw(network_graph, slayout, with_labels=True, node_size=500, node_color='lightblue')
            st.pyplot(plt)

        with col2:
            st.subheader("Centrality Measures")
            centralities = get_centralities(network_graph)
            for centrality_name, centrality_values in centralities.items():
                top_5 = sorted(centrality_values.items(), key=lambda x: -x[1])[:5]
                st.write(f"**{centrality_name}**")
                st.write(top_5)
    else:
        st.error("No data found for the given protein ID.")
