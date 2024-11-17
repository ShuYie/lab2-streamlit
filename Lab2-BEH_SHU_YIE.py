import requests
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import streamlit as st

# Layout for title and buttons in a single row
column1, column2 = st.columns([3, 1])

# Display the title in the first column
with column1:
    st.title('Lab 2 - BEH SHU YIE')

# Display the buttons in the second and third columns
with column2:
    choice = st.radio("Select PPI Databases", options=["Biogrid", "String"], index=0)

# Input for Uniprot ID
protein_id = st.text_input('Enter Protein ID')
retrieve = st.button('Retrieve')

# Functions to retrieve data and generate network
def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "50047376b5c8bab7d27cc48417f95b38",
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,  # Example: human (Homo sapiens)
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    network = response.json()
    network_df = pd.DataFrame.from_dict(network, orient='index')
    network_df.OFFICIAL_SYMBOL_A = [gene.upper() for gene in network_df.OFFICIAL_SYMBOL_A]
    network_df.OFFICIAL_SYMBOL_B = [gene.upper() for gene in network_df.OFFICIAL_SYMBOL_B]
    return network_df

def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606  # Example: human (Homo sapiens)
    }
    response = requests.get(string_url, params=params)
    network = response.json()
    network_df = pd.json_normalize(network)
    network_df['preferredName_A'] = network_df['preferredName_A'].str.upper()
    network_df['preferredName_B'] = network_df['preferredName_B'].str.upper()
    return network_df

def generate_network(dataframe):
    network_graph = nx.Graph()
    if 'OFFICIAL_SYMBOL_A' in dataframe.columns and 'OFFICIAL_SYMBOL_B' in dataframe.columns:
        network_graph = nx.from_pandas_edgelist(dataframe, "OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B")
    elif 'preferredName_A' in dataframe.columns and 'preferredName_B' in dataframe.columns:
        network_graph = nx.from_pandas_edgelist(dataframe, "preferredName_A", "preferredName_B")
    else:
        st.warning("Invalid PPI data!")
    return network_graph

def get_centralities(network_graph):
    degree_centrality = nx.degree_centrality(network_graph)
    betweenness_centrality = nx.betweenness_centrality(network_graph)
    closeness_centrality = nx.closeness_centrality(network_graph)
    eigenvector_centrality = nx.eigenvector_centrality(network_graph, 300)
    pagerank_centrality = nx.pagerank(network_graph)
    
    return {
        'Degree Centrality': degree_centrality,
        'Betweenness Centrality': betweenness_centrality,
        'Closeness Centrality': closeness_centrality,
        'Eigenvector Centrality': eigenvector_centrality,
        'PageRank Centrality': pagerank_centrality,
    }


# Retrieve and display data based on button selection
if retrieve:
    
    if protein_id != "":
        
        if choice == "Biogrid":
            data = retrieve_ppi_biogrid(protein_id)
        else:
            data = retrieve_ppi_string(protein_id)

        if not data.empty:
            col1, col2 = st.columns(2)
            network_graph = generate_network(data)
            slayout = nx.spring_layout(network_graph, seed=123)
            node_size = 20 if choice == "Biogrid" else 1000

            with col1:
                st.subheader('\nPPI Data Information')
                st.dataframe(data)
                st.write(f"**Number of edges:** {network_graph.number_of_edges()}")
                st.write(f"**Number of nodes:** {network_graph.number_of_nodes()}")
                plt.figure(figsize=(8, 6)) 
                nx.draw(network_graph, slayout, with_labels=False, node_size=node_size, node_color='lightblue')
                st.pyplot(plt)

            with col2:
                st.subheader('\nCentrality Measures')
                centrality = get_centralities(network_graph)

                for measure, values in centrality.items():
                    st.write(f"**{measure}:**")
                    top_5_nodes = sorted(values.items(), key=lambda x: -x[1])[:5]
                    for index, (node, value) in enumerate(top_5_nodes, start=1):
                        st.write(f"{index}. {node}: {value:.4f}")

                    top_5_node_names = [node for node, _ in top_5_nodes]

                    if measure == 'Degree Centrality':
                        st.write(f" Proteins '{top_5_node_names}' have higher Degree Centrality, indicating they will interact with many other proteins, suggesting they may be a hub proteins. Hub proteins are often crucial in cellular processes and are more likely to be essential for survival or associated with disease pathways.")

                    elif measure == 'Betweenness Centrality':
                        st.write(f" Proteins '{top_5_node_names}' have higher Betweenness Centrality, indicating they often serve as connectors or regulators of information flow between different parts of the network, which can make them important for maintaining network integrity and facilitating communication between modules or clusters.")

                    elif measure == 'Closeness Centrality':
                        st.write(f" Proteins '{top_5_node_names}' have higher Closeness Centrality, indicating they influence many parts of the network efficiently and are often involved in signal transduction or regulatory functions")

                    elif measure == 'Eigenvector Centrality':
                        st.write(f" Proteins '{top_5_node_names}' have higher Eigenvector Centrality, indicating they are the influential proteins that are connected to other influential proteins, potentially marking out proteins that play central roles within highly active or essential functional regions of the network.")

                    elif measure == 'PageRank Centrality':
                        st.write(f" Proteins '{top_5_node_names}' have PageRank Centrality, indicating they are well-connected and linked to other prominent proteins, which helping to highlight functionally essential or central proteins in signaling pathways or disease-related clusters.")

                plt.figure(figsize=(8, 6))
                nx.draw(network_graph, slayout, with_labels=False, node_size=node_size, node_color='lightblue')
                nx.draw_networkx_nodes(network_graph, slayout, nodelist=top_5_node_names, node_size=node_size*5, node_color='orange')  
                st.pyplot(plt)

            st.write('In a biological context, using multiple centrality measures together can help prioritize proteins that might be essential to network stability, information flow, or biological function, potentially identifying targets for further research or therapeutic intervention.')
    else:
        st.warning('Please enter Protein ID')