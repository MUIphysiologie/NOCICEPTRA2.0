import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import altair as alt
import duckdb


@st.cache_resource
def load_data():
    """
    Defines the database connection to the NOCICEPTRA duckdb database
    """
    try:
        print("Try to connect to the database")
        return duckdb.connect(database = "./Data/nociceptra.duckdb", read_only = True)
    except Exception as e:
        print(f"Error: {e}")


def exploratory_data_analysis():
    """ start the analysis
    args:
        data_dict: dict -> dictionary with all dataframe
    """
    con = load_data()
    experimental_description(con)


def experimental_description(con):
    """  initialize the analysis and draw networks for mirnas, mRNAs and trajectories of modules

    Args:
        sup_mrna: pd.DataFrame -> mrna supercluster expression data
        sup_mirna: pd.DataFrame -> mirna supercluster expression data
        mirna_edges: pd.DataFrame -> target spaces for miRNAs
        gene_enr: pd.DataFrame -> enrichments for mRNA clusters
        mirna_enr: pd.DataFrame -> enrichments for miRNA cluster
    """

    tab1, tab2 = st.tabs(["WGCNA mRNA expression","WGCNA miRNA expression"])
    # retrieve the cluster groups


    #set the columns
    columns = ["D0_S1","D0_S2","D0_S3",
               "D5_S1","D5_S2","D5_S3",
               "D9_S1","D9_S2","D9_S3",
               "D16_S1","D16_S2","D16_S3",
               "D26_S1","D26_S2","D26_S3",
               "D36_S1","D36_S2","D36_S3"
               ]

    index= ["D00","D00","D00",
            "D05","D05","D05",
            "D09","D09","D09",
            "D16","D16","D16",
            "D26","D26","D26",
            "D36","D36","D36"
            ]

    #mRNA trajectories
    col1, col2 = tab1.columns(2)
    stages = con.execute("Select hierachical_cluster from sup_mrna").fetchdf()["hierachical_cluster"].unique()

    sel_stage = col1.selectbox("Please choose the differentiation stage:", stages)

    # retrieve the right stage
    modules = con.execute(f"Select cluster from sup_mrna WHERE hierachical_cluster='{sel_stage}'").fetchdf()["cluster"].unique()
    # list the module stages for the frontend
    module_sel = col2.selectbox("Choose WGCNA module: ", modules)

    #retrieve the genes of the module via selection
    sup_mrna_mod = con.execute(f"Select * from sup_mrna WHERE cluster='{module_sel}'").fetchdf()
    sup_mrna_mod = sup_mrna_mod.set_index("gene_name").iloc[:,1:-1]
    sup_mrna_mod.columns = index
    sup_mrna_mod = pd.melt(sup_mrna_mod.reset_index(), id_vars = "gene_name").set_index("gene_name")

    #gene enrichments for the selected cluster
    mod_gene_enr = con.execute(f"Select * from gprofiler_enr WHERE cluster='{module_sel}'").fetchdf()
    #check if user wants to search for hub-genes

    # draw the plot and the network
    fig1 = plot_writing(sup_mrna_mod,"miRNA Module Trajectory") # put the axis within the context
    interactive_figure = draw_network(con, module_sel)

    # correctly order it into the app
    col1.plotly_chart(interactive_figure, use_container_width = True)
    col2.altair_chart(fig1, use_container_width = True)
    tab1.write("---")
    tab1.write("Enrichment Analysis:")
    tab1.dataframe(mod_gene_enr.set_index("description")[["p-value","source"]].iloc[:10,:], use_container_width=True)

    # this should still go into a separate function
    #miRNA trajectories for each module
    modules_mirna = con.execute("Select cluster from sup_mirna").fetchdf()["cluster"].unique()
    tab2.write("Select your Module of Interest:")
    module_mirna_sel = tab2.selectbox("", modules_mirna)
    tab2.write("-----")

    # retrieve the modules for miRNA analysis
    sup_mirna = con.execute(f"Select * from sup_mirna WHERE cluster='{module_mirna_sel}'").fetchdf().set_index("gene_name").iloc[:,1:-1]
    sup_mirna.columns = index
    sup_mirna_mod = pd.melt(sup_mirna.reset_index(), id_vars = "gene_name").set_index("gene_name")
    mod_mirna_enr = con.execute(f"Select * from mirna_enr WHERE mirna_moduel='{module_mirna_sel}'").fetchdf()
    sup_mirna.columns = columns

    # set up the layout
    mirna_1, mirna_2 = tab2.columns(2)
    # fig for trajectories without hub-genes
    mirna_fig = plot_writing(sup_mirna_mod, "miRNA Module Trajectory")

    # write the figurse
    mirna_1.altair_chart(mirna_fig, use_container_width = True)
    mirna_2.dataframe(mod_mirna_enr.set_index("Pathway")[["p_value","Term"]].iloc[:10,:],
                      use_container_width=True)


    # put the dataframe into the expander and style it like a heatmap
    tab2.write("-------")
    heatmap_exp =  tab2.expander("Show miRNA DataFrame:")
    heatmap_exp.write("------")
    heatmap_exp.header("miRNA DataFrame")
    heatmap_exp.write(sup_mirna.style.background_gradient(cmap ='viridis', axis = 1)\
        .set_properties(**{'font-size': '15px'}).set_caption("Hello World"))

def plot_writing(mean_mrna, title):
    """ write a lineplot for the expresssion counts
    mean_mrna: pd.DataFrame -> dataframe that contains the mean expression"""
    return draw_altair_graph(mean_mrna, title)


def draw_network(con, selection, ax = None):
    """Show Top 30 hub-gene networks for each module
    Drawing is based on Correlation above 0.7 and correlation is determined using pearson correlation

    Args:
        hub_genes: list -> list of top hub genes
        expression: pd.DataFrame -> expression matrix for correlation analysis
        selection: st.multiselect list -> gene selection in gui
        ax: where to plot
    returns:
        interactive_figure --> plotly figure
    """

    hub_genes_module = con.execute(f"Select external_gene_name from hub_genes WHERE Module='{selection}'").fetchdf()["external_gene_name"].tolist()
    # retrieve the correlation matrix for futher network analysis using the correlatio as weight
    hub_genes_module = tuple(hub_genes_module)
    genes_corr = con.execute(f"Select * from sup_mrna WHERE gene_name IN {hub_genes_module}").fetchdf().set_index("gene_name").iloc[:,2:-1].T.corr(method = "pearson")
    #genes_corr = expression[expression["external_gene_name"].isin(hub_genes_module)].set_index("external_gene_name").iloc[:,2:-1].T.corr(method = "pearson") #correlation of top 30 hubs
    corr_net = genes_corr.stack() # stack the correlation to get a linkage table
    corr_net.index = corr_net.index.set_names(['gene',"gene2"])
    corr_net = corr_net.reset_index()
    corr_net.columns = ["target","source","corr"]
    corr_net = corr_net.loc[(corr_net["corr"] != 1) & (corr_net["corr"] > 0.7)]

    #draw the network
    network_hub = nx.from_pandas_edgelist(corr_net, "target", "source","corr")
    return draw_interactive_network(network_hub)

@st.cache_data()
def hub_genes():
    """load the hub-genes, this might take a while an will be cached to avoid long loading times """
    return pd.read_parquet("./Data/module_hubs_kme.parquet")

def draw_interactive_network(gene_network):
    """ Function to fill the interactive network graph using plotly from the gene network
    constructed with networkx

    Args:
        gene_network: networkx -> adjacency matrix of network x graph
    Returns:
        interactive_figure -> plotly figure
    """
    module = "skyblue"
    degree = gene_network.degree()
    pos = nx.spring_layout(gene_network)


    edge_x = []
    edge_y = []

    for edge in gene_network.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend((x0, x1, None))
        edge_y.extend((y0, y1, None))
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    node_text = []
    node_hover = []

    for node in gene_network.nodes():
        node_text.append(f"<b>{node}</b>")
        node_hover.append(node)
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    node_link = [f"<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene={i}\">{i}</a>" for i in node_hover]
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        hoverinfo ="text",
        text = "text",
        textposition='top center',
        textfont=dict(color='black',),
        mode='markers+text',
        marker=dict(
            color=module,
            size=30
           ))


    node_adjacencies = [
        len(adjacencies[1]) for adjacencies in gene_network.adjacency()
    ]
    #node_trace.marker.color = node_adjacencies
    node_trace.text = node_link

    fig = go.Figure(data=[edge_trace, node_trace],
        layout=go.Layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False,
        hovermode='closest',
        margin=dict(l=5,r=5,t=10),
        annotations=[ dict(
            showarrow=True,
            xref="paper", yref="paper",
            x=0.05, y=-0.002 ) ],
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, ticks = ""),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, ticks = ""))
    )

    fig.update_layout(
            hoverlabel=dict(
                bgcolor="darkgrey",
                font_size=16,
                align = "auto",
            )
        )

    return fig


def draw_altair_graph(data_draw,title, gene_annotation = None):
    """ Draw the altai graph boxplot with lines
    args:
        data_draw: pd.DataFrame -> dataframe to draw from
        title: str -> tilte of the plot
        gene_annotation -> title of the legend
    returns:
        altair graph
    """
    data_draw.columns = ["Timepoint","z-score variance stabilized counts"]
    chart = (
                alt.Chart(data_draw)
                .mark_boxplot(size = 30)
                .encode(x=alt.X("Timepoint"), y=alt.Y("z-score variance stabilized counts"))
                .interactive()
                .properties(title = title)
                )

    chart_line = (
            alt.Chart(data_draw)
            .mark_line(interpolate = "natural")
            .encode(x="Timepoint", y="mean(z-score variance stabilized counts)")
            .interactive()
            .properties())

    final = chart_line + chart
    final  = final.configure_legend(padding=10,
                                    cornerRadius=5,
                                    orient='bottom').configure_axis(grid = False, labelFontSize = 13).configure_view(strokeOpacity = 0)
    return final



if __name__ == "__main__":
    exploratory_data_analysis()

