import contextlib
import streamlit as st
import pandas as pd
import plotly.express as px
import http
import urllib
import re
import statsmodels.api as sm
import holoviews as hv
from holoviews import opts, dim
import requests
import duckdb
import logging
# Create a logger
logger = logging.getLogger('my_logger')

# Set the level of the logger. This can be DEBUG, INFO, WARNING, ERROR, or CRITICAL
logger.setLevel(logging.DEBUG)

# Create a stream handler for the logger
stream_handler = logging.StreamHandler()

# Create a formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Set the formatter for the stream handler
stream_handler.setFormatter(formatter)

# Add the stream handler to the logger
logger.addHandler(stream_handler)
hv.extension('bokeh', logo=False) #draw interactive visualization using holoviews



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


def kegg_disease_analysis():
    con = load_data()
    tab1, tab2, tab3 = st.tabs(["KEGG Network","Disease Network", "miRNA Networks (BETA)"])
    selected_kegg, disease_selection, selected_mirnas = execute_fill_list_pathways(con)
    kegg_decision = tab1.selectbox("Please select your KEGG pathway:",selected_kegg)
    kegg_id = con.execute(f"Select gene_name, ID from kegg_pathway WHERE gene_name='{kegg_decision}'").fetchdf()["ID"].tolist()[0]
    resulting_gene_list = get_kegg(kegg_id)
    super_cluster_statistics = con.execute("Select * from super_cluster_statistics").fetchdf()
    super_cluster_counts = con.execute("Select * from super_cluster_counts").fetchdf()
    statistics, p_value = gene_set_kegg_enrichment(resulting_gene_list,
                                                   super_cluster_statistics,
                                                   super_cluster_counts)


    st.sidebar.info("Please select you KEGG/Disease pathway or the miRNA of interest to derive enriched genes throughout iPSC-derived sensory neuron development")

    make_chart(resulting_gene_list,
               con,
               statistics = statistics,
               p_value = p_value,
               tab = tab1
               )


    # also evaluate the disease gene interaction
    diesase_selected = tab2.selectbox("Please select the Disease of Interest:", disease_selection)
    show_labels = tab2.checkbox("Show labels in rich plots", value = True)
    run_disease_analysis(con, diesase_selected, show_labels,super_cluster_statistics,super_cluster_counts, tab2)


    #
    mirna = tab3.selectbox("Choose your miRNA:", selected_mirnas)
    target_score = tab3.slider("Target-Score Treshold (considers only values above the threshold)", min_value = 0, max_value = 200, step = 1, value = 75)

    if tab3.button("Start miNRA analysis:"):
        get_mirna_information(mirna, con,target_score, tab3)


def run_disease_analysis(con,diesase_selected, show_labels, super_cluster_statistics, super_cluster_counts, tab):

    if tab.button("Start Analysis:"):

        try:
            disease_genes = con.execute(
                                "SELECT geneSymbol FROM diseases WHERE diseaseName=?", 
                                (disease_selected,)
                            ).fetchnumpy()["geneSymbol"].tolist()
        except Exception as e:
            logger.error(e)
            tab.warning("Disease currently not found")
            return None

        #disease_genes = dataframe_dictionary["disease"][dataframe_dictionary["disease"]["diseaseName"] == diesase_selected]["geneSymbol"].unique()
        disstats, dis_pval = gene_set_kegg_enrichment(disease_genes, super_cluster_statistics,
                                                super_cluster_counts)
        make_chart(disease_genes,con,
                statistics = disstats,
                p_value = dis_pval,
                tab = tab,
                checkbox = show_labels,
                threshold = 0.95
                )



def get_kegg(kegg_pathway):
    """ REST API:
    query a kegg pathway to obtain all the genes in the pathway
    kegg_pathway -> ID of the pathway
    """
    #get the patway
    pathway = kegg_pathway
    url = f"http://rest.kegg.jp/get/{pathway}"
    genes_liste =[] # make an empty gene list that will be filled with the genes
    #connect to KEGG URL

    http.client.HTTPConnection._http_vsn = 10
    http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0' # avoid version problems arising between http and python


    with urllib.request.urlopen(url) as f:
        lines = f.read().decode('utf-8').splitlines()
        want = 0

        #get the kegg pathways and all the genes belonging to the kegg pathways
        for line in lines:
            fields = line.split()
            ## The list of genes starts here
            if fields[0] == 'GENE':
                want = 1
                ## The line with GENE is different
                genes_liste.append(fields[2].rstrip(';'))
            ## We reached the next section of the file
            elif want == 1 and re.match('^\S', line):
                return genes_liste;
            ## We're still in the list of genes
            if want == 1 and len(fields)>1:
                genes_liste.append((fields[1].rstrip(';')))
        return genes_liste

def gene_set_kegg_enrichment(genes_list, genes_query, genes_baseline):

    """ Here we write a function to define contigency table of observed genes and expected genes"""
    genes_diagram = genes_baseline.set_index("stage")
    genes_expected = genes_query[genes_query["external_gene_name"].isin(genes_list)] # queried genes
    genes_expected = pd.DataFrame(genes_expected["supercluster_gene"].value_counts())

    #merging of genes from the kegg enrichments
    gesamt = pd.merge(genes_diagram, genes_expected, how = "left", left_index = True, right_index = True)
    gesamt = gesamt.fillna(0)

    #table statistics and summary
    table = sm.stats.Table(gesamt)
    rslt = table.test_nominal_association()
    enrichments_residuals = pd.DataFrame(table.resid_pearson)
    enrichments_residuals.columns = ["observed_genes","expected_genes"]

    #get the p-value
    p_value = rslt.pvalue # get the p-value
    bars = enrichments_residuals.iloc[:,1] # retrieve the pearson residuals
    return bars, p_value

def make_chart(genes, con, disease_mirna = None,  mirna = None, statistics = None, p_value = None, mirna_name = None, tab = None, checkbox = True, threshold = 0.4):

    """
    miRNA and mRNA network analyis, build a bokeh plot
    get PPI and miRNA::mRNA interaction table --> index
    get all members of PPI, miRNA::mRNA --> node.data

    """

    links = get_interaction(genes, con, threshold).drop("score", axis = 1)
    node = get_node_table_hv(con, links, genes)

    links.columns = ["source","target","value"]
    links = links[links["source"].isin(node["gene_name"])]
    links = links[links["target"].isin(node["gene_name"])]
    color_liste = {"pluripotency":"green","differentiation_induction/maturation":"orange",
                    "early differentiation":"skyblue", "early neural progentior":"pink"
                    ,"neural progenitor 1": "purple", "neuroal progenitor late" : "rose","nociceptor":"red"}
    node["color"] = [color_liste.get(i) for i in node["supercluster_gene"]]
    node = node[node["gene_name"].isin(links["source"].tolist() + links["target"].tolist())]
    nodes = hv.Dataset(pd.DataFrame(node.sort_values("supercluster_gene")),"index")

    links_index = get_indeces_hv(links, nodes)
    col_enr1, col_enr2 = tab.columns(2)
    # here single miRNA interactions are determined

    chord_diagram = chart_plot(nodes, links_index, checkbox = checkbox)
    enrichment_diagram = enrichments_genes(nodes.data["gene_name"].tolist(),"hsapiens")

    col_enr1.bokeh_chart(hv.render(chord_diagram, backend = "bokeh"), use_container_width = True)
    draw_table_info(p_value,col_enr2)

    col_enr2.markdown(r""" **Table 1**:  Standardized Pearson residuals""")
    col_enr2.dataframe(statistics, use_container_width = True)

    if enrichment_diagram.shape[0] > 0:
        col_enr2.markdown(r""" **Table 2**: G:Profiler enrichments""")
        col_enr2.table(enrichment_diagram.set_index("description"))

    else:
        tab.warning("No enrichments found check your internet connection")

def get_node_table_hv(con, links, genes):
    """
    """
    genes = tuple(genes)
    mirna_genes = con.execute(f"Select gene_name, supercluster_gene, score from chord_diagram WHERE gene_name IN {genes} AND score > 50").fetchdf()
    mirna_genes = mirna_genes.groupby(["gene_name","supercluster_gene"])["score"].agg("sum").reset_index()
    node = mirna_genes[mirna_genes["gene_name"].isin(links["preferred_name_x"])]
    node = node[node["gene_name"].isin(links["preferred_name_y"])]

    return node

def get_indeces_hv(links, nodes):
    """
    """

    index_data = nodes.iloc[:,:2]
    links_index = pd.merge(links,index_data.data, how = "left", left_on = "source",right_on = "gene_name")
    links_index = pd.merge(links_index, index_data.data, how = "left", left_on = "target", right_on = "gene_name")
    links_index = links_index.drop(["source","target","gene_name_x","gene_name_y"], axis = 1)
    links_index = links_index[["index_x","index_y","value"]]
    links_index.columns = ["source","target","value"]
    links_index = pd.merge(links_index,nodes.data, how = "left", left_on = "source", right_on = "index").drop(["index","gene_name","score"],axis=1)
    links_index.columns = ["source","target","value","trial","color"]
    links_index = links_index.sort_values("trial")
    links_index = links_index[~links_index[['source', 'target']].apply(frozenset, 1).duplicated()]

    return links_index

def draw_table_info(p_value: float, colum_sel: str):
    """
    Should draw the appropriate Table Information choosen by the p-value detected

    input:
        p_value: float <- the tested p-value for the network enrichemt
        column_sel <- st.column <- column where the information should be drawn with
    """
    if p_value < 0.05:
        colum_sel.markdown(r''' **Figure 1:** Chord-plot of custom Gene interaction, colors define the  differentiation stages.
        Size of the Gene dots is defined by the cummulative miRNA target-score and edges weight are defined by the StringDB database confidence score.
        The $chi-squared$ p-value of $\textbf{%f}$, indicates an significant enrichment of genes''' %(p_value), unsafe_allow_html=True)
    else:
        colum_sel.markdown(r''' **Figure 1** Chord-plot of custom Gene interaction, colors define the  differentiation stages.
        Size of the Gene dots is defined by the cummulative miRNA target-score and edges weight are defined by the StringDB database confidence score. <br>
        The $chi-squared$ p-value of $\textbf{%f}$, indicates no significance''' %(p_value), unsafe_allow_html=True)

def get_interaction(genes: pd.DataFrame,con, threshold:float = 0.4):

    """ use the string DB database

    - frozen set
    - remove duplicates that are just exchanges
    - prefered_name_x --> source
    - preferred_name_y --> target
    - genes -> input
    """
    genes = tuple(genes)
    return con.execute(
        f'Select * from string_interaction_high WHERE preferred_name_x IN {genes} AND preferred_name_y IN {genes} AND score>{threshold};'
    ).fetchdf()

def chart_plot(node: pd.DataFrame, index: list, mirna: bool = None, checkbox: bool = True):

    """ two different networks
    mirna -> if miRNA is queried or not
    node -> a dataset of all nodes involved in the network
    index -> a network adjacency matrix of all interactions
    """

    return (
        hv.Chord((index, node)).opts(
            opts.Chord(
                hooks=[hook],
                cmap="Paired",
                edge_color=dim("trial").str(),
                node_color=dim('supercluster_gene').str(),
                edge_line_width=1,
                node_size=dim("score") / 400,
                fontsize={"labels": 0.5},
                node_line_width=1,
                height=500,
                labels=dim("gene_name").str(),
            )
        )
        if checkbox is True
        else hv.Chord((index, node)).opts(
            opts.Chord(
                hooks=[hook],
                cmap="Paired",
                edge_color=dim("trial").str(),
                node_color=dim('supercluster_gene').str(),
                edge_line_width=1,
                fontsize={"labels": 0.5},
                node_line_width=1,
                height=500,
            )
        )
    )


def set_toolbar_autohide(plot, element):
    bokeh_plot = plot.state
    bokeh_plot.toolbar.autohide = True

def enrichments_genes(genes: list, species: str) -> pd.DataFrame:

    """ use the G:Profiler enrichment analysis REST-api
    genes -> queried genes
    species --> hsapiens but can be changed have a look gprofiler
    """
    import json
    go_profiler = {}
    df_go_end = pd.DataFrame()  # select a table to append
    genes = genes #reference to the list of genes

    # request the gprofiler website
    r = requests.post(
        url= 'https://biit.cs.ut.ee/gprofiler/api/gost/profile',
        json={
        'organism':species,
        'query': genes,
        'sources' :["GO:BP","GO:MF","GO:CC","KEGG"]},
        headers={
        'User-Agent':'FullPython'
        },
    )
    #st.write(r.raise_for_status())
    if r.status_code == 200: # check if connections worked
        with contextlib.suppress(KeyError):
            data = r.json()["result"]
            parents_list = []
            go_list = []
            # Check for the best per parent pathways the most specialized pathways with the most detailed description
            for n in data:
                go_list.append(n["native"])
                parents_list.extend(iter(n["parents"]))
            end_list = [i for i in go_list if i not in parents_list]

            p_value = [m["p_value"] for m in data if m["native"] in end_list]
            desc_value = [l["name"] for l in data if l["native"] in end_list]
            source_list = [l["source"] for l in data if l["native"] in end_list]
            # update the dictionary
            go_profiler |= {
                "p-value": p_value,
                "go-terms": end_list,
                "description": desc_value,
                "source": source_list,
            }

            df_go = pd.DataFrame(columns = ["go-terms","description","source","p-value"])
            df_go["go-terms"] = list(end_list)
            df_go["description"] = desc_value
            df_go["source"] = source_list
            df_go["p-value"] = p_value
            df_mf = df_go[df_go["source"] == "GO:MF"].sort_values("p-value").iloc[:10,:]
            df_bp = df_go[df_go["source"] == "GO:BP"].sort_values("p-value").iloc[:10,:]
            df_cc = df_go[df_go["source"] == "GO:CC"].sort_values("p-value").iloc[:10,:]
            df_kegg = df_go[df_go["source"] == "KEGG"].sort_values("p-value").iloc[:10,:]
            df_go = pd.concat([df_mf,df_bp,df_cc,df_kegg])
            df_go = df_go.drop(["go-terms"], axis = 1)
            df_go_end = df_go_end.append(df_go)

    else:
        print("enrichment is not working currently")

    return df_go_end

def get_mirna_information(mirna, con, target_score, tab):

    """
    positive --> above positive correlation (default = false, means below the queried correlation)
    mirna --> queried miRNA
    mirna_genes --> all miRNAs miRNA_edges Database
    target_score --> all miRNA::mRNAs targets above the target score
    """


    #statistics_mirna = mirna_genes[(mirna_genes["score"] >= target_score) & (mirna_genes["correlation"] < -0.7)]
    mirnas_stats = con.execute(f"Select * from chord_diagram WHERE score>{target_score} AND correlation<-0.7").fetchdf()
    #mirna_genes = mirna_genes[(mirna_genes["score"] >= target_score) & (mirna_genes["correlation"] < -0.7)]

    targeted_mirna = mirnas_stats[mirnas_stats["mirna"]==mirna]

    mirna_target_genes = set(targeted_mirna["gene_name"].tolist())
    statistics,p_value = mirna_enrichments_statistics(mirnas_stats, mirna)
    statistics = statistics.iloc[:,1]
    #draw a plot for the statistics

    #check how lenghty the genes list
    if len(mirna_target_genes) > 10:

        make_chart(mirna_target_genes,con,
                   disease_mirna = None,
                   mirna = True,
                   statistics = statistics,
                   p_value = float(p_value),
                   mirna_name = mirna,
                   tab = tab)
    else:
        #else the number of targets is to loo
        tab.warning("Try different settings or a different miRNA, the number of miRNA targets is to low")

def mirna_enrichments_statistics(mirna_targets, mirna):
    """ Function to detect gene_enrichments of miRNA targetign

    mirna_targets -> mirna_prediction,correlation,validation table
    mirna -> queried mirna

    returns: enrichment for the 7 different stages as table, float(p_value) """


    #detect the number of genes belonging to each cluster, remove duplicates
    genes_diagram = mirna_targets[["gene_name","supercluster_gene"]]
    genes_diagram = genes_diagram.drop_duplicates(keep = "first")["supercluster_gene"].value_counts()

    #search for the queried miRNA and count cluster occurences
    mirna_diagram = mirna_targets.loc[mirna_targets["mirna"].str.contains(mirna, case = False)]
    mirna_diagram = mirna_diagram[["gene_name","supercluster_gene"]]
    mirna_diagram = mirna_diagram.drop_duplicates(keep = "first")["supercluster_gene"].value_counts()

    #merge both tables
    end = pd.merge(pd.DataFrame(genes_diagram),mirna_diagram, how = "left", left_index = True, right_index = True)
    end = end.fillna(0)
    ending = end
    #ending = ending[["supercluster_gene_y", "supercluster_gene_x"]]
    #statistical analysis via contigency table
    table = sm.stats.Table(ending)
    rslt = table.test_nominal_association()
    enrichments_residuals = pd.DataFrame(table.resid_pearson)
    enrichments_residuals.columns = ["observed_genes","expected_genes"]
    p_value = rslt.pvalue


    return enrichments_residuals, p_value

def hook(plot, example):
    """_summary_: This should handle the plot background

    Args:
        plot (_type_): _description_
        example (_type_): _description_
    """
    print("function is working")
    print('plot.state:   ', plot.state)
    print('plot.handles: ', sorted(plot.handles.keys()))
    plot.handles["plot"].background_fill_alpha= 0
    plot.handles["plot"].border_fill_alpha= 0
    plot.handles["text_1_glyph"].text_color = "grey"

@st.cache_resource
def execute_fill_list_pathways(_con):
    """_summary_: Loads the selectooxes with the appropiate data

    Args:
        _con (DuckDBConnection): _description_

    Returns:
        _type_: _description_
    """
    kegg = _con.execute("Select gene_name from kegg_pathway").fetchnumpy()["gene_name"].tolist()
    disease = _con.execute("Select diseaseName from diseases").fetchdf()["diseaseName"].unique()
    mirna = _con.execute("Select mirna from chord_diagram").fetchdf()["mirna"].unique()
    return (kegg,disease, mirna)

if __name__ == "__main__":
    kegg_disease_analysis()