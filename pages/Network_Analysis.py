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

hv.extension('bokeh', logo=False) #draw interactive visualization using holoviews

@st.cache
def load_data(data_path):

    """
    """

    kegg_pathway = pd.read_parquet(data_path + "kegg_pathway.parquet")
    super_cluster_statistics = pd.read_parquet(data_path + "gene_general.parquet")
    super_cluster_counts = pd.read_parquet(data_path + "gene_general_counting.parquet")
    chord_diagram = pd.read_parquet(data_path + "mirna_edges.parquet")
    string_interaction = pd.read_parquet(data_path +"string_hsa_preprocess.parquet")
    string_interaction["score"] = string_interaction["combined_score"].apply(lambda x: x/1000)
    string_interaction_high = string_interaction.loc[string_interaction["score"] >= 0.4]
    diseases = pd.read_parquet(data_path + "all_gene_disease_associations.parquet")

    dataframe_dictionary = {"kegg_pathway": kegg_pathway,
                            "super_cluster_statistics":super_cluster_statistics,
                            "super_cluster_counts": super_cluster_counts,
                            "string_interaction_high": string_interaction_high,
                            "chord_diagram":chord_diagram,
                            "disease": diseases}
  
    return dataframe_dictionary

def kegg_disease_analysis(dataframe_dictionary):
    tab1, tab2, tab3 = st.tabs(["KEGG Network","Disease Network", "miRNA Networks (BETA)"])
    kegg_decision = tab1.selectbox("Please select your KEGG pathway:", dataframe_dictionary["kegg_pathway"]["gene_name"].tolist())
    kegg_id = dataframe_dictionary["kegg_pathway"][dataframe_dictionary["kegg_pathway"]["gene_name"] == kegg_decision]["ID"].tolist()[0]
    resulting_gene_list = get_kegg(kegg_id)
    statistics, p_value = gene_set_kegg_enrichment(resulting_gene_list,
                                                   dataframe_dictionary["super_cluster_statistics"],
                                                   dataframe_dictionary["super_cluster_counts"])


    st.sidebar.info("Please select you KEGG/Disease pathway or the miRNA of interest to derive enriched genes throughout iPSC-derived sensory neuron development")

    make_chart(resulting_gene_list, 
               dataframe_dictionary["chord_diagram"],
               dataframe_dictionary["string_interaction_high"],
               statistics = statistics,
               p_value = p_value,
               tab = tab1
               )


    # also evaluate the disease gene interaction
    disease_selection = dataframe_dictionary["disease"]["diseaseName"].unique()
    diesase_selected = tab2.selectbox("Please select the Disease of Interest:", disease_selection)
    show_labels = tab2.checkbox("Show labels in rich plots", value = True)

    # here we check the label size
    if diesase_selected:
        if tab2.button("Create Disease Network"):
            tab2.write("------")
            disease_genes = dataframe_dictionary["disease"][dataframe_dictionary["disease"]["diseaseName"] == diesase_selected]["geneSymbol"].unique()
            disstats, dis_pval = gene_set_kegg_enrichment(disease_genes, 
                                                        dataframe_dictionary["super_cluster_statistics"],
                                                        dataframe_dictionary["super_cluster_counts"])
            make_chart(disease_genes, 
                    dataframe_dictionary["chord_diagram"],
                    dataframe_dictionary["string_interaction_high"][dataframe_dictionary["string_interaction_high"]["score"] > 0.95],
                    statistics = disstats,
                    p_value = dis_pval,
                    tab = tab2,
                    checkbox = show_labels
                    ) 

    # 
    mirna = tab3.selectbox("Choose your miRNA:", dataframe_dictionary["chord_diagram"]["mirna"].unique())
    target_score = tab3.slider("Target-Score Treshold (considers only values above the threshold)", min_value = 0, max_value = 200, step = 1, value = 75)

    get_mirna_information(mirna,dataframe_dictionary["chord_diagram"], dataframe_dictionary["string_interaction_high"],target_score, tab3)

def get_kegg(kegg_pathway):
    """ REST API: 
    query a kegg pathway to obtain all the genes in the pathway
    kegg_pathway -> ID of the pathway 
    """
    #get the patway
    pathway = kegg_pathway
    url = "http://rest.kegg.jp/get/%s"%(pathway) #retrieve the url from the kegg database
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
    genes_diagram = genes_baseline.set_index("Unnamed: 0.1")
    genes_expected = genes_query[genes_query["external_gene_name"].isin(genes_list)] # queried genes
    genes_expected = pd.DataFrame(genes_expected["supercluster_gene"].value_counts())

    #merging of genes from the kegg enrichments
    gesamt = pd.merge(genes_diagram, genes_expected, how = "left", left_index = True, right_index = True)
    gesamt = gesamt.fillna(int(0))

    #table statistics and summary
    table = sm.stats.Table(gesamt)
    rslt = table.test_nominal_association()
    enrichments_residuals = pd.DataFrame(table.resid_pearson)
    enrichments_residuals.columns = ["observed_genes","expected_genes"]

    #get the p-value
    p_value = rslt.pvalue # get the p-value
    bars = enrichments_residuals.iloc[:,1] # retrieve the pearson residuals
    return bars, p_value

def make_chart(genes, mirna_genes, string_interaction_high, disease_mirna = None,  mirna = None, statistics = None, p_value = None, mirna_name = None, tab = None, checkbox = True):

    """
    miRNA and mRNA network analyis, build a bokeh plot
    get PPI and miRNA::mRNA interaction table --> index
    get all members of PPI, miRNA::mRNA --> node.data

    """

    links = get_interaction(genes, string_interaction_high).drop("score", axis = 1)
    cluster = mirna_genes[["gene_name","supercluster_gene"]].drop_duplicates()
    node = get_node_table_hv(mirna_genes, links, genes, cluster)


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
    
    col_enr1.bokeh_chart(hv.render(chord_diagram, backend = "bokeh"))
    draw_table_info(p_value,col_enr2)
    
    col_enr2.markdown(r""" **Table 1**:  Standardized Pearson residuals""")
    col_enr2.table(statistics)

    if enrichment_diagram.shape[0] > 0:
        col2.markdown(r""" **Table 2**: G:Profiler enrichments""")
        col2.table(enrichment_diagram.set_index("description"))

    else:
        tab.warning("No enrichments found check your internet connection")
        
def get_node_table_hv(mirna_genes, links, genes, cluster):
    """
    """
    mirna_genes = mirna_genes[mirna_genes["score"]>50]
    mirna_genes = mirna_genes[mirna_genes["gene_name"].isin(genes)]
    mirna_degree = mirna_genes.copy(deep = True)
    mirnas = mirna_genes.copy(deep = True)
    gene_mirna = mirna_genes.groupby(["mirna"])["score"].mean()
    count_mirna = mirna_genes["mirna"].value_counts()
    mirna_genes = mirna_genes.groupby(["gene_name"])["score"].agg("sum").reset_index()
    mirna_genes = pd.merge(mirna_genes,cluster, how = "left", left_on = "gene_name",right_on = "gene_name")
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

def draw_table_info(p_value, colum_sel):
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

def get_interaction(genes,string):

    """ use the string DB database

    - frozen set 
    - remove duplicates that are just exchanges
    - prefered_name_x --> source
    - preferred_name_y --> target 
    - genes -> input 
    """

    string_network = string[string["preferred_name_x"].isin(genes)]
    string_end = string_network[string_network["preferred_name_y"].isin(genes)]
    return string_end

def chart_plot(node, index, mirna = None, checkbox = True):
 
    """ two different networks
    mirna -> if miRNA is queried or not
    node -> a dataset of all nodes involved in the network
    index -> a network adjacency matrix of all interactions
    """

    if checkbox is True:
    #draw th echord plot
        chord = hv.Chord((index, node)).opts(opts.Chord(hooks=[set_toolbar_autohide],
                                            cmap="Paired",
                                            edge_cmap="tab20",
                                            edge_color=dim("trial").str(),
                                            node_color=dim('supercluster_gene').str(),
                                            edge_line_width = 1,
                                            width = 500, 
                                            height = 500, 
                                            node_size = dim("score")/400, 
                                            fontsize = {"labels":0.5},
                                            node_line_width = 1, 
                                            labels = dim("gene_name").str()))
        return chord

    else:
        chord = hv.Chord((index, node)).opts(opts.Chord(hooks=[set_toolbar_autohide],
                                            cmap="Paired",
                                            edge_cmap="tab20",
                                            edge_color=dim("trial").str(),
                                            node_color=dim('supercluster_gene').str(),
                                            edge_line_width = 1,
                                            width = 500, 
                                            height = 500, 
                                            fontsize = {"labels":0.5},
                                            node_line_width = 1
                                            ))
        return chord


def set_toolbar_autohide(plot, element):
    bokeh_plot = plot.state
    bokeh_plot.toolbar.autohide = True    

def enrichments_genes(genes, species):

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
        try: # catch errors if 404
            data = r.json()["result"]
            parents_list = []
            go_list = []
            p_value = []
            desc_value = []
            source_list = []

            # Check for the best per parent pathways the most specialized pathways with the most detailed description
            for n in data:
                go_list.append(n["native"]) 
                for t in n["parents"]:
                    parents_list.append(t)
            end_list = [i for i in go_list if i not in parents_list]

            for m in data:
                # cehck for the 
                if m["native"] in end_list:
                    p_value.append(m["p_value"])
            for l in data:
                if l["native"] in end_list:
                    desc_value.append(l["name"])
            for l in data:
                if l["native"] in end_list:
                    source_list.append(l["source"])

            # update the dictionary
            go_profiler.update({"p-value": p_value, "go-terms":end_list, "description":desc_value, "source": source_list})

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

        except KeyError as e:
            pass
    else:
        print("enrichment is not working currently")

    return df_go_end

def get_mirna_information(mirna, mirna_genes,string_interaction_high, target_score, tab):

    """ 
    positive --> above positive correlation (default = false, means below the queried correlation)
    mirna --> queried miRNA
    mirna_genes --> all miRNAs miRNA_edges Database
    target_score --> all miRNA::mRNAs targets above the target score
    """


    statistics_mirna = mirna_genes[(mirna_genes["score"] >= target_score) & (mirna_genes["correlation"] < -0.7)]
    mirna_genes = mirna_genes[(mirna_genes["score"] >= target_score) & (mirna_genes["correlation"] < -0.7)]

    targeted_mirna = mirna_genes[mirna_genes["mirna"]==mirna]

    mirna_target_genes = set(targeted_mirna["gene_name"].tolist())
    statistics,p_value = mirna_enrichments_statistics(statistics_mirna, mirna)
    statistics = statistics.iloc[:,1]
    #draw a plot for the statistics

    #check how lenghty the genes list
    if len(mirna_target_genes) > 10:
        
        make_chart(mirna_target_genes, mirna_genes,
                   string_interaction_high,
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
    end = end.fillna(int(0)) #fill the nan with 0 for further analysis
    ending = end
    #ending = ending[["supercluster_gene_y", "supercluster_gene_x"]]
    #statistical analysis via contigency table
    table = sm.stats.Table(ending)
    rslt = table.test_nominal_association()
    enrichments_residuals = pd.DataFrame(table.resid_pearson)
    enrichments_residuals.columns = ["observed_genes","expected_genes"]
    p_value = rslt.pvalue
    
   
    return enrichments_residuals, p_value

if __name__ == "__main__":
    data_path = "./Data/"
    dataframe_dictionary = load_data(data_path)
    kegg_disease_analysis(dataframe_dictionary)