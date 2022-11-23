import streamlit as st 
import pandas as pd
import plotly.express as px
from scipy.stats import zscore
import altair as alt
import duckdb
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly.express as px


@st.experimental_singleton
def load_data():
    """
    should load the parquet files to retrieve the tables
    """
    try:
        con = duckdb.connect(database = "./Data/nociceptra.duckdb", read_only = True)
        return con
    except Exception as e:
        print(f"Error: {e}")
    
def trajectory_start(con):

    """
    Make the selection of the genes and draw the trajectories that belong to the genes
    
    input:
        supercluster_mrna: pd.DataFrame -> all vst scores
        tpm_end: pd.DataFrame -> all gene tpms
        
    """
    tab1, tab2, tab3, tab4 = st.tabs(["mRNA and miRNA Trajectories",
                                      "miRNA with Multimapper",
                                      "ncRNA Trajectories (Excerpt)",
                                      "lncRNA Trajectories"]) 
    
    selected_draw = sorted(con.execute("SELECT gene_name from vsd_counts").fetchnumpy()["gene_name"])

    #retrieve the genes for selection for searching
    gene_queried = tab1.multiselect("Gene or miRNA:", selected_draw)
    tab1.write("---")

    #retrieve the mirnas for selection for searching
    mirna_select = con.execute("SELECT gene_name from fivep_counts").fetchnumpy()["gene_name"].tolist() + con.execute("SELECT gene_name from threep_counts").fetchnumpy()["gene_name"].tolist()
    mirna_select = sorted(list(set(mirna_select)))
    mirna_queried = tab2.multiselect("Select miRNA: ", mirna_select)
    tab2.write("---")
    #ncRNA 
    nc_select = con.execute("SELECT gene_name from ncrna_counts").fetchnumpy()["gene_name"].tolist()
    nc_queried = tab3.multiselect("Select ncRNA: ", nc_select)
    tab3.write("---")
    
    lnc_select = con.execute("Select external_gene_name from lnc_counts").fetchdf()
    lnc_queried = tab4.multiselect("Select your lncRNA of interest: ", lnc_select["external_gene_name"].tolist())

    # here the layout of the sidebar should be added
    draw_information_layout()

    if gene_queried:
        preprocess_tpm_vsd(gene_queried,con,tab1) 

    if mirna_queried:
        mirna_multimap_drawing(mirna_queried,con, tab2)

    if nc_queried:
        nc_multimap_drawing(nc_queried,con,tab3)
        
    if lnc_queried:
        lnc_query_genes(lnc_queried, con, tab4)
        
        
        
def preprocess_tpm_vsd(genes_liste,con,tab):

    """ Preprocess the Tables of genes for normalized counts and tpm before drawing
    args:
        df: pd.DataFrame -> normalized counts data for mRNAs 
        tpm_end: pd.DataFrame -> transcript per million for mRNAs
        tab: selected tab to draw figures
    """

    # make a tabbed sheet  # set name column as index
    columns = ["DAY00","DAY00","DAY00",
               "DAY05","DAY05","DAY05",
               "DAY09","DAY09","DAY09",
               "DAY16","DAY16","DAY16",
               "DAY26","DAY26","DAY26",
               "DAY36","DAY36","DAY36"]
    genes = tuple(genes_liste)
    df_curves = con.execute(f"SELECT * from vsd_counts WHERE gene_name IN {genes}").fetchdf().set_index("gene_name")
    df_curves.columns = columns
    
    df_curves = pd.melt(df_curves.reset_index(), id_vars = "gene_name") # melt the table
    df_curves.columns = ["Gene Name","Timepoint","z-scored variance stabilized counts"]

    # preprocess the tpm table
    tpm_curves = con.execute(f"SELECT * from tpm_counts WHERE gene_name IN {genes}").fetchdf().set_index("gene_name")
    tpm_curves.columns = columns
    tpm_curves = pd.melt(tpm_curves.reset_index(),id_vars = "gene_name")
    tpm_curves.columns = ["Gene Name","Timepoint","TPM (Transcript per Million)"]
    
    # check the shapes of the dataframe 
    if tpm_curves.shape[0] > 0:
        make_trajectories(df_curves, tpm_curves, tab)

    else:
        make_trajectories(df_curves, None, tab)
    
    
    #make the cellline plot
    
    col1, col2 = tab.columns(2)
    col1.write("---")
    col2.write("---")
    genes_queried = con.execute(f"Select * from all_counts WHERE gene_name IN {tuple(genes_liste)}").fetchdf().set_index("gene_name")
    meta_data_table = con.execute("Select * from mrna_metadata").fetchdf()
    cell_line_specific_printing(genes_queried, meta_data_table, col1)
    genes_liste = tuple(genes_liste)
    if len(genes_liste) >= 3:
        correlation_matrix_analysis(genes_queried, genes_liste, col2)
        
    elif len(genes_liste) == 2:
        scatter_comparison(genes_queried, genes_liste, col2)
    
    
def cell_line_specific_printing(genes_queried,metadata_table, col1):
    columns = [3* [i] for i in ["Day00", "Day05", "Day09", "Day16", "Day26", "Day36"]] *3
    columns = [t for i in columns for t in i]
    
    metadata_table["Timepoint"] = columns
    gene_cell_line = genes_queried.copy().reset_index()
    melted_queried = pd.melt(gene_cell_line, id_vars="gene_name")
    melted_queried = pd.merge(melted_queried, metadata_table, left_on = "variable", right_on = "samples",how = "left")
    cell_select = col1.selectbox("Select Cell-Line:", ["AD2","AD3","840"])
    selected_table = melted_queried[melted_queried["cell_line"] == cell_select]
    fig = draw_altair_graph(selected_table, "value", "gene_name")
    col1.write(fig)
    
  
def scatter_comparison(genes_queried, genes_liste, col2):
    """"""
    genes_queried = genes_queried.T
    max = genes_queried.to_numpy().max()
    min = genes_queried.to_numpy().min()
    fig = alt.Chart(genes_queried).mark_circle(size=60).encode(
    x= alt.X(str(genes_liste[0]), scale=alt.Scale(domain=[min-0.2, max+0.2])),
    y=alt.Y(str(genes_liste[1]), scale=alt.Scale(domain=[min-0.2, max+0.2])),
        ).interactive().properties(width=550, height = 400)
    
    final  = fig.configure_legend(padding=2,
                                    cornerRadius=5,
                                    orient='bottom').configure_axis(grid = False, 
                                                                      labelFontSize = 13).configure_view(strokeOpacity = 0)
    #final_figure = fig + fig.transform_regression(str(genes_liste[0]),str(genes_liste[1])).mark_line()
    col2.write(final)
    
    
def correlation_matrix_analysis(genes_queried, genes_liste, col2):
    """"""
    
    genes_queried = genes_queried.T.corr()
    mask = np.triu(np.ones_like(genes_queried, dtype=bool))
    fig = px.imshow(genes_queried, 
                    color_continuous_scale='Viridis',
                    text_auto=True,
                    aspect="auto",
                    title="Correlation Matrix between queried Genes")
    col2.write(fig)
                       
def draw_information_layout():
    st.sidebar.subheader("InfoBox")
    st.sidebar.info("Please Search for your gene/miRNA or ncRNA of interest by using the multiSelect panel shown in each Tab. The Panels will be automatically updated!")

def make_trajectories(df_curves, tpm_curves = None, tab = None):

    """ this function is to draw the gene and miRNA trajectories """
    
    col1, col2 = tab.columns(2)
    if tpm_curves is not None:
        #draw lineplots for the data with standard errors
            #vsd_fig = px.box(df_curves,x = "variable", y = "value", color = "Unnamed: 0")
        vsd_figure = draw_altair_graph(df_curves, "z-scored variance stabilized counts", "Gene Name")
        tpm_figure = draw_altair_graph(tpm_curves, "TPM (Transcript per Million)","Gene Name")
        # write the figure into the ap
        col1.empty()
        col1.write(vsd_figure)
        col2.write(tpm_figure)

    else:
        st.markdown("---")
        vsd_figure = draw_altair_graph(df_curves,"z-scored variance stabilized counts", "Gene Name")
        col1.write(vsd_figure)

def mirna_multimap_drawing(searched_mirna, con,tab):
    """ Make drawings for the miRNA multimappers
    args:
        five: pd.DataFrame -> data for five prime miRNAs
        three: pd.DataFrame -> data for three prime miRNAs
        tab: selected tab to draw plots
    """

    col1, col2 = tab.columns(2)
    searched_mirna = tuple(searched_mirna)
    columns_zelline = ["gene_name","Day00","Day00","Day00",
                       "Day05","Day05","Day05",
                       "Day09","Day09","Day09",
                       "Day16","Day16","Day16",
                       "Day26","Day26","Day26",
                       "Day36","Day36","Day36"]
    
    
    five = con.execute(f"SELECT * from fivep_counts WHERE gene_name IN {searched_mirna}").fetchdf()
    three = con.execute(f"SELECT * from threep_counts WHERE gene_name IN {searched_mirna}").fetchdf()
    five.columns = columns_zelline
    three.columns = columns_zelline
    
    five_p_mirna = pd.melt(five, id_vars=["gene_name"])
    five_p_mirna.columns = ["miRNA","Timepoint","Raw Counts"]
    three_p_mirna = pd.melt(three, id_vars=["gene_name"])
    three_p_mirna.columns = ["miRNA","Timepoint","Raw Counts"]
    
    if (five_p_mirna.shape[0] > 0 and three_p_mirna.shape[0] > 0):
        five_figure = draw_altair_graph(five_p_mirna,"Raw Counts", "miRNA")
        three_figure = draw_altair_graph(three_p_mirna,"Raw Counts","miRNA")
    
        col1.write(five_figure)
        col2.write(three_figure)
        
    elif (five_p_mirna.shape[0] == 0 and three_p_mirna.shape[0] > 0):
        col1.info("No counts identified for five prime miRNA")
        three_figure = draw_altair_graph(three_p_mirna,"Raw Counts","miRNA")
        col2.write(three_figure)
        
    elif (five_p_mirna.shape[0] > 0 and three_p_mirna.shape[0] == 0):
        col2.info("No counts identified for three prime miRNA")
        five_figure = draw_altair_graph(five_p_mirna,"Raw Counts", "miRNA")
        col1.write(five_figure)
         
    else:
        st.error("Connection Timeout, or selected miRNA not present in list!")
        
def nc_multimap_drawing(nc_queried, con,tab):
    """
    """

    col1, col2 = tab.columns(2)
    columns = [3* [i] for i in ["Day00", "Day05", "Day09", "Day16", "Day26", "Day36"]] *3
    columns = [t for i in columns for t in i]
    columns = ["gene_name"] + columns
    nc_queried = tuple(nc_queried)
    nc_data = con.execute(f"SELECT * from ncrna_counts WHERE gene_name IN {nc_queried}").fetchdf()
    nc_data.columns = columns
    nc_data = nc_data.set_index("gene_name")
    nc_data_query = nc_data.apply(zscore, axis = 1)
    nc_data_query = nc_data_query.reset_index()
    nc_melt = pd.melt(nc_data_query, id_vars=["gene_name"])
    nc_melt.columns = ["ncRNA", "Timepoint", "z-scored variance stabilized counts"]
    nc_figure = draw_altair_graph(nc_melt, "z-scored variance stabilized counts", gene_annotation = "ncRNA")
    col1.write(nc_figure)
    
    
    nc_sig = con.execute(f"SELECT * from significance WHERE ncRNA IN {nc_queried}").fetchdf()
    
    #nc_sig.applymap(lambda x: 'background-color: red' if x > 0.01 else 'background-color: green', subset='padj')
    col2.write("ncRNA Differential Gene Expression Information:")
    col2.dataframe(nc_sig)
    
def draw_altair_graph(data_draw, value, gene_annotation = None):
    selection = alt.selection_multi(fields=[gene_annotation], bind='legend')
    chart = (
                alt.Chart(data_draw)
                .mark_boxplot(size = 30)
                .encode(x=alt.X("Timepoint"),
                        y=alt.Y(value, scale=alt.Scale(domain=[data_draw[value].min()-1, data_draw[value].max()+1])),
                        color=gene_annotation,
                        opacity=alt.condition(selection, alt.value(1), alt.value(0.2)))
                .interactive()
                .properties(width=550)
                .add_selection(selection)
                )
        
    chart_line = (
            alt.Chart(data_draw)
            .mark_line(interpolate = "natural")
            .encode(x=alt.X("Timepoint",
                            scale=alt.Scale(padding=1)),
                            y=alt.Y(f"mean({value})",scale=alt.Scale(domain=[data_draw[value].min()-1, data_draw[value].max()+1])),
                            color=gene_annotation,
                            opacity=alt.condition(selection, alt.value(1), alt.value(0.2)),
                            size=alt.condition(~selection, alt.value(1), alt.value(3)))
            .interactive()
            .properties(width=550, height = 400)
            .add_selection(selection))
                                                
    final = chart_line + chart
    final  = final.configure_legend(padding=2,
                                    cornerRadius=5,
                                    orient='bottom').configure_axis(grid = False, 
                                                                      labelFontSize = 13).configure_view(strokeOpacity = 0)
    return final

def lnc_query_genes(genes_liste, con, tab):
    columns = [9*[i] for i in ["Day00", "Day05", "Day09", "Day16", "Day26", "Day36"]]
    columns = [t for i in columns for t in i]
    columns = ["gene_name"] + columns
   
    selected_lncs = con.execute(f"Select * from lnc_counts WHERE external_gene_name IN {tuple(genes_liste)}").fetchdf()
    selected_lncs.columns = columns
    melted_selected_lncs = pd.melt(selected_lncs, id_vars="gene_name")
    melted_selected_lncs.columns = ["gene_name","Timepoint", "vsd counts"]
    fig = draw_altair_graph(melted_selected_lncs, "vsd counts", "gene_name" )
    tab.write("---")
    tab.write(fig)



if __name__ == "__main__":
    connection = load_data()
    trajectory_start(connection) 