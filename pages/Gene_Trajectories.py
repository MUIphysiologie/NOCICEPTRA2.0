import streamlit as st 
import pandas as pd
import plotly.express as px
from scipy.stats import zscore
import altair as alt

@st.cache
def load_data(data_path):
    """
    should load the parquet files to retrieve the tables
    """

    columns = [[i]*3 for i in ["Day00", "Day05","Day09","Day16", "Day26", "Day36"]] * 3
    columns = [t for i in columns for t in i]

    columns_to_skip = ["Unnamed: 0"]
    dataframe_dictionary = {}
    supercluster_mrna = pd.read_parquet(data_path + "all_mirna_mRNA_trajectories.parquet").iloc[:21467,:]

    # transcript per million loading
    tpm_end = pd.read_parquet(data_path + "cell-line_average_mean_tpm.parquet")
    tpm_end = tpm_end.set_index("gene_name")

    # load multimapper miRNA
    five_p = pd.read_parquet(data_path + "miRNA_nc_multimappers_5p.parquet").drop("mean", axis = 1).set_index("Name")
    three_p = pd.read_parquet(data_path + "miRNA_nc_multimappers_3p.parquet").drop("mean", axis = 1).set_index("Name")

    nc_rna_counts = pd.read_parquet(data_path + "vsd_counts_ncRNA_excerpt.parquet").set_index("Unnamed: 0")
    nc_rna_counts.columns = columns
    
    
    significance = pd.read_parquet(data_path +"excerpt_sig.parquet").drop("Unnamed: 0", axis = 1)
    significance.columns = ["ncRNA", "baseMean", "FDR-correct p-value"]
    significance["significant?"] = ["significant" if i < 0.05 else "not significant" for i in significance["FDR-correct p-value"]]

    dataframe_dictionary.update({"supercluster_mrna":supercluster_mrna,
                                 "tpm_end": tpm_end, 
                                 "five_p": five_p, 
                                 "three_p": three_p, 
                                 "ncRNA_counts": nc_rna_counts,
                                 "significance": significance})

    return dataframe_dictionary

def trajectory_start(dataframe_dict):

    """
    Make the selection of the genes and draw the trajectories that belong to the genes
    
    input:
        supercluster_mrna: pd.DataFrame -> all vst scores
        tpm_end: pd.DataFrame -> all gene tpms
        
    """
    
    tab1, tab2, tab3 = st.tabs(["mRNA and miRNA Trajectories", "miRNA with Multimapper", "ncRNA Trajectories (Excerpt)"]) 
    selected_draw = dataframe_dict["supercluster_mrna"]["Unnamed: 0"].tolist()

    #retrieve the genes for selection for searching
    gene_queried = tab1.multiselect("Gene or miRNA:", selected_draw)
    tab1.write("---")

    #retrieve the mirnas for selection for searching
    mirna_select = dataframe_dict["five_p"].index.tolist() + dataframe_dict["three_p"].index.tolist()
    mirna_select = sorted(list(set(mirna_select)))
    mirna_queried = tab2.multiselect("Select miRNA: ", mirna_select)
    tab2.write("---")
    #ncRNA 
    nc_select = dataframe_dict["ncRNA_counts"].index.tolist()
    nc_queried = tab3.multiselect("Select ncRNA: ", nc_select)
    tab3.write("---")

    # here the layout of the sidebar should be added
    draw_information_layout()

    if gene_queried:
        preprocess_tpm_vsd(gene_queried, dataframe_dict["supercluster_mrna"],dataframe_dict["tpm_end"],tab1) 

    if mirna_queried:
        mirna_multimap_drawing(mirna_queried, dataframe_dict["five_p"], dataframe_dict["three_p"], tab2)

    if nc_queried:
        nc_multimap_drawing(nc_queried, dataframe_dict["ncRNA_counts"],dataframe_dict["significance"],tab3)

def preprocess_tpm_vsd(genes_liste, df, tpm_end,tab):

    """ Preprocess the Tables of genes for normalized counts and tpm before drawing
    args:
        df: pd.DataFrame -> normalized counts data for mRNAs 
        tpm_end: pd.DataFrame -> transcript per million for mRNAs
        tab: selected tab to draw figures
    """

    # make a tabbed sheet 
    df_curv = df.set_index("Unnamed: 0") # set name column as index
    columns = ["DAY00","DAY00","DAY00",
               "DAY05","DAY05","DAY05",
               "DAY09","DAY09","DAY09",
               "DAY16","DAY16","DAY16",
               "DAY26","DAY26","DAY26",
               "DAY36","DAY36","DAY36"]
    df_curv.columns = columns #unified Data columns for making the line and scatterplots
    df_curves = df_curv[df_curv.index.isin(genes_liste)].sort_index().reset_index()
    df_curves = pd.melt(df_curves, id_vars = "Unnamed: 0") # melt the table
    df_curves.columns = ["Gene Name","Timepoint","z-scored variance stabilized counts"]

    # preprocess the tpm table
    tpm_end.columns = columns
    tpm_curves = tpm_end[tpm_end.index.isin(genes_liste)].sort_index().reset_index()
    tpm_curves = pd.melt(tpm_curves,id_vars = "gene_name")
    tpm_curves.columns = ["Gene Name","Timepoint","TPM (Transcript per Million)"]
    
    # check the shapes of the dataframe 
    if tpm_curves.shape[0] > 0:
        make_trajectories(df_curves, tpm_curves, tab)

    else:
        make_trajectories(df_curves, None, tab)
  
def draw_information_layout():
    st.sidebar.subheader("InfoBox")
    st.sidebar.info("Please Search for your gene/miRNA or ncRNA of interest by using the multiSelect panel shown in each Tab. The Panels will be automatically updated!")

def make_trajectories(df_curves, tpm_curves = None, tab = None):

    """ this function is to draw the gene and miRNA trajectories """
    
    col1, col2 = tab.columns(2)
    if tpm_curves is not None:
        #draw lineplots for the data with standard errors
            #vsd_fig = px.box(df_curves,x = "variable", y = "value", color = "Unnamed: 0")
        print(tpm_curves, df_curves)
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

def mirna_multimap_drawing(searched_mirna, five, three,tab):
    """ Make drawings for the miRNA multimappers
    args:
        five: pd.DataFrame -> data for five prime miRNAs
        three: pd.DataFrame -> data for three prime miRNAs
        tab: selected tab to draw plots
    """

    col1, col2 = tab.columns(2)
    columns_zelline = ["Day00","Day00","Day00","Day05","Day05","Day05","Day09","Day09","Day09","Day16","Day16","Day16","Day26","Day26","Day26","Day36","Day36","Day36"]
    five.columns = columns_zelline
    three.columns = columns_zelline
    five_p_mirna = five[five.index.isin(searched_mirna)].reset_index()
    three_p_mirna = three[three.index.isin(searched_mirna)].reset_index()
    five_p_mirna = pd.melt(five_p_mirna, id_vars=["Name"])
    five_p_mirna.columns = ["miRNA","Timepoint","Raw Counts"]
    three_p_mirna = pd.melt(three_p_mirna, id_vars=["Name"])
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
        
def nc_multimap_drawing(nc_queried, nc_data, significance,tab):
    """
    """

    col1, col2 = tab.columns(2)
    nc_data_query = nc_data.apply(zscore, axis = 1)
    nc_data_query = nc_data_query[nc_data_query.index.isin(nc_queried)].reset_index()
    nc_melt = pd.melt(nc_data_query, id_vars=["Unnamed: 0"])
    nc_melt.columns = ["ncRNA", "Timepoint", "z-scored variance stabilized counts"]
    nc_figure = draw_altair_graph(nc_melt, "z-scored variance stabilized counts", gene_annotation = "ncRNA")
    col1.write(nc_figure)
    
    
    nc_sig = significance[significance["ncRNA"].isin(nc_queried)]
    
    #nc_sig.applymap(lambda x: 'background-color: red' if x > 0.01 else 'background-color: green', subset='padj')
    col2.write("ncRNA Differential Gene Expression Information:")
    col2.dataframe(nc_sig)
    
    


def draw_altair_graph(data_draw, value, gene_annotation = None):
    selection = alt.selection_multi(fields=[gene_annotation], bind='legend')
    chart = (
                alt.Chart(data_draw)
                .mark_boxplot(size = 30)
                .encode(x=alt.X("Timepoint"),
                        y=alt.Y(value),
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
                            y=alt.Y(f"mean({value})"),
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




if __name__ == "__main__":

    data_path = "./Data/"   
    dataframe_dictionary = load_data(data_path)
    trajectory_start(dataframe_dictionary) 