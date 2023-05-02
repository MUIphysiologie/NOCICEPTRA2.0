import streamlit as st
import pandas as pd
import plotly.express as px
from scipy.stats import zscore
import altair as alt
import duckdb
import numpy as np
import plotly.express as px

@st.cache_resource
def load_data():
    """
    Defines the database connection to the NOCICEPTRA duckdb database
    """
    try:
        return duckdb.connect(database = "./Data/nociceptra.duckdb", read_only = True)
    except Exception as e:
        print(f"Error: {e}")

def trajectory_start():

    """
    Make the selection of the genes and draw the trajectories that belong to the genes

    input:
        supercluster_mrna: pd.DataFrame -> all vst scores
        tpm_end: pd.DataFrame -> all gene tpms

    """
    con = load_data()
    st.header("Expression Signatures of iPSC-derived sensory Neuron (iDNs)")
    tab1, tab2, tab3, tab4 = st.tabs(["mRNA and miRNA Trajectories",
                                      "miRNA with Multimapper",
                                      "ncRNA Trajectories (Excerpt)",
                                      "lncRNA Trajectories"])

    selected_draw, mirna_select,nc_select, lnc_select = execute_fill_list(con)

    tab1.write("---")
    #retrieve the genes for selection for searching
    gene_queried = tab1.multiselect("Select Gene:", selected_draw)
    tab1.write("---")

    #retrieve the mirnas for selection for searching
    mirna_select = sorted(list(set(mirna_select)))
    mirna_queried = tab2.multiselect("Select miRNA: ", mirna_select)
    tab2.write("---")
    #ncRNA
    nc_queried = tab3.multiselect("Select ncRNA: ", nc_select)
    tab3.write("---")
    lnc_queried = tab4.multiselect("Select lncRNA: ", lnc_select["external_gene_name"].tolist())

    # here the layout of the sidebar should be added

    if gene_queried:
        preprocess_tpm_vsd(gene_queried,con,tab1)

    if mirna_queried:
        mirna_multimap_drawing(mirna_queried,con, tab2)

    if nc_queried:
        nc_multimap_drawing(nc_queried,con,tab3)

    if lnc_queried:
        lnc_query_genes(lnc_queried, con, tab4)

    draw_information_layout()

def preprocess_tpm_vsd(genes_liste: list,con: duckdb, tab: st.tabs) -> None:

    """ Preprocess the Tables of genes for normalized counts and tpm before drawing
    args:
        df: pd.DataFrame -> normalized counts data for mRNAs
        tpm_end: pd.DataFrame -> transcript per million for mRNAs
        tab: selected tab to draw figures
    """

    # make a tabbed sheet  # set name column as index
    genes = tuple(genes_liste)
    df_curves, df_original = get_counts_tables(con, genes)
    metadata_table = con.execute("Select * from MetaDataTable").fetchdf()
    df_curves = pd.merge(df_curves, metadata_table, how = "left", left_on = "samples", right_on = "samples")
    print(df_curves)
    tpm_curves = get_tpm_tables(con,genes)
    # check logic if heatmap or boxplot

    if tpm_curves.shape[0] > 0:
        if len(df_curves["Gene Name"].unique()) <= 5:
            make_trajectories(df_curves, tpm_curves, tab)
        else:
            fig = make_heatmap(df_curves, "Gene Name", "z-scored variance stabilized counts", "Timepoint", title = "Gene Signatures")
            tab.altair_chart(fig, use_container_width = True)

    elif len(df_curves["Gene Name"].unique()) <= 5:
        make_trajectories(df_curves, None, tab)
    else:
        fig = make_heatmap(df_curves, "Gene Name", "z-scored variance stabilized counts", "Timepoint", title = "Gene Signatures")
        tab.altair_chart(fig, use_container_width = True)


    #make the cellline plot
    col1, col2 = tab.columns(2)
    col1.write("---")
    col2.write("---")


    # here we check per celline
    cell_line_specific_printing(df_curves, metadata_table, col1)
    genes_liste = tuple(genes_liste)

    # correlation matrix will be performed whenever more than 3 genes are select
    # if there only two genes draw a scatter plot
    if len(genes_liste) >= 3:
        correlation_matrix_analysis(df_original, genes_liste, col2)

    elif len(genes_liste) == 2:
        scatter_comparison(df_original, genes_liste, col2)


def get_counts_tables(con, genes_liste: list)-> pd.DataFrame:
    """_summary_: Retrieves the count tables from the Database

    Args:
        con (DuckDBConnection): Connection to the database
        genes_list (list): The selected genes to query

    Returns:
        pd.DataFrame: DataFrame long format with the zscored counts
    """
    df_original = con.execute(f"Select * from all_counts WHERE gene_name IN {tuple(genes_liste)}").fetchdf().set_index("gene_name")
    df_curves = df_original.apply(zscore, axis = 1).reset_index()
    df_curves = pd.melt(df_curves, id_vars = "gene_name") # melt the table
    df_curves.columns = ["Gene Name","samples","z-scored variance stabilized counts"]
    return df_curves, df_original


def get_tpm_tables(con: duckdb, genes: list) -> pd.DataFrame:
    """_summary_: Retrieves the transcript per million tables from the database
    using the selected genes

    Args:
        con (DuckDBConnection): Connection to the database
        genes (list): The selected genes to query

    Returns:
        pd.DataFrame: DataFrame long format with the tpm values
    """
    columns_tpm = ["DAY00_1","DAY00_2","DAY00_3",
               "DAY05_1","DAY05_2","DAY05_3",
               "DAY09_1","DAY09_2","DAY09_3",
               "DAY16_1","DAY16_2","DAY16_3",
               "DAY26_1","DAY26_2","DAY26_3",
               "DAY36_1","DAY36_2","DAY36_3"]

    tpm_curves = con.execute(f"SELECT * from tpm_counts WHERE gene_name IN {genes}").fetchdf().set_index("gene_name")
    tpm_curves.columns = columns_tpm
    tpm_curves = pd.melt(tpm_curves.reset_index(),id_vars = "gene_name")
    tpm_curves.columns = ["Gene Name","Timepoint","TPM (Transcript per Million)"]
    return tpm_curves

def cell_line_specific_printing(df_curves,metadata_table, col1):
    """ Print the trajectoreis for the selected genes for the cell line specifically
    Args:
        genes_queried: list -> selected genes
        metadata_table: pd.DataFrame --> table that hold the sample and timepoint information
        col1: st.tabs.col --> tab and column where this figure should be drawn
    """

    cell_select = st.sidebar.selectbox("Select Cell-Line:", ["AD2","AD3","840"])
    selected_table = df_curves[df_curves["cell_line"] == cell_select]

    if len(selected_table["Gene Name"].unique()) <= 5:
        fig = draw_altair_graph(selected_table, "z-scored variance stabilized counts", "Gene Name", "Time aggregated")
    else:
        fig = make_heatmap(selected_table, "Gene Name", "z-scored variance stabilized counts", "Timepoint","Cell-Type specific Signatures")
    col1.altair_chart(fig, use_container_width = True)
    st.sidebar.markdown("---")

def scatter_comparison(genes_queried: list, genes_liste:list, col2):
    # sourcery skip: avoid-builtin-shadow
    """ two selected genes scatter correlation analysis
    Args:
        genes_queried: list -> selected genes from mulitselect
        genes_list: list
        cols2: st.tabs.cols
    """
    #this removes the duplicated
    genes_queried = genes_queried.groupby(level=0)
    genes_queried = genes_queried.last() # this retrieves the last value per group
    genes_queried = genes_queried.T

    # retrieve the max and minimum values
    max = genes_queried.to_numpy().max()
    min = genes_queried.to_numpy().min()
    fig = alt.Chart(genes_queried).mark_circle(size=60).encode(
    x= alt.X(str(genes_liste[0]), scale=alt.Scale(domain=[min-0.2, max+0.2])),
    y=alt.Y(str(genes_liste[1]), scale=alt.Scale(domain=[min-0.2, max+0.2])),
        ).interactive().properties()

    final  = fig.configure_legend(padding=2,
                                    cornerRadius=5,
                                    orient='bottom').configure_axis(grid = False,
                                                                      labelFontSize = 13).configure_view(strokeOpacity = 0)
    #final_figure = fig + fig.transform_regression(str(genes_liste[0]),str(genes_liste[1])).mark_line()
    col2.altair_chart(final, use_container_width = True)


def correlation_matrix_analysis(genes_queried: list, genes_liste:list, col2):
    """
    Draws the correlation matrix form the selected data
    Args:
        genes_queried: list -> list of selected genes
    """

    corrMatrix = prepare_correlation_matrix(genes_queried)
    chart = alt.Chart(corrMatrix).mark_rect().encode(
    x=alt.X('var1', title=None),
    y=alt.Y('var2', title=None),
    color=alt.Color('correlation', legend=None),
                    ).properties(
                        width=alt.Step(40),
                        height=alt.Step(40),
                        title = "Correlation Matrix"
                    )

    chart += chart.mark_text(size=12).encode(
        text=alt.Text('correlation', format=".2f"),
        color=alt.condition(
            "datum.correlation > 0.5",
            alt.value('white'),
            alt.value('black')
        )
    )
    col2.altair_chart(chart, use_container_width = True)

def prepare_correlation_matrix(genes_queried: list):
    """
    Prepares the correlation matrix for the selected data
    Args:
        genes_queried: list -> list of selected genes
    returns:
        corrMatrix: pd.DataFrame -> correlation matrix

    """
    corrMatrix = genes_queried.T.corr().stack()
    corrMatrix.index.names = ["var1","var2"]
    corrMatrix = corrMatrix.reset_index()
    corrMatrix.columns = ['var1', 'var2', 'correlation']
    corrMatrix = corrMatrix.dropna()
    return corrMatrix

def draw_information_layout():
    """_summary_ Should draw the Info Box
    """
    st.sidebar.subheader("InfoBox")
    st.sidebar.info("Please Search for your gene/miRNA or ncRNA of interest by using the multiSelect panel shown in each Tab. The Panels will be automatically updated!")

def make_trajectories(df_curves: pd.DataFrame, tpm_curves: pd.DataFrame = None, tab: st.tabs = None):
    """_summary_: Draws the trajectories for the selected genes using altair

    Args:
        df_curves (pd.DataFrame): Holding vst counts
        tpm_curves (pd.Dataframe, optional): Holding tpm values. Defaults to None.
        tab (st.tabs, optional): Tab where the figure should be drawn. Defaults to None.
    """
    col1, col2 = tab.columns(2)

    if tpm_curves is not None:
        #draw lineplots for the data with standard errors
            #vsd_fig = px.box(df_curves,x = "variable", y = "value", color = "Unnamed: 0")
        draw_trajectories(tpm_curves, df_curves, col1, col2)
    else:
        st.markdown("---")
        vsd_figure = draw_altair_graph(df_curves,"z-scored variance stabilized counts", "Gene Name", "Time aggregated")
        col1.altair_chart(vsd_figure, use_container_width = True)
        col2.warning("No Data found for TPM here!")


# TODO Rename this here and in `make_trajectories`
def draw_trajectories(tpm_curves: pd.DataFrame,
                      df_curves: pd.DataFrame,
                      col1: st.columns,
                      col2: st.columns):
    """_summary_: Draws the trajectories for the selected genes using altair"""
    tpm_curves["Timepoint"] = [i.split("_")[0] for i in tpm_curves["Timepoint"]]
    vsd_figure = draw_altair_graph(df_curves, "z-scored variance stabilized counts", "Gene Name", "Time aggregated")
    tpm_figure = draw_altair_graph(tpm_curves, "TPM (Transcript per Million)","Gene Name")
    # write the figure into the ap
    col1.empty()
    col1.altair_chart(vsd_figure,use_container_width = True)
    col2.altair_chart(tpm_figure,use_container_width= True)

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

        col1.altair_chart(five_figure,use_container_width = True)
        col2.altair_chart(three_figure, use_container_width = True)

    elif (five_p_mirna.shape[0] == 0 and three_p_mirna.shape[0] > 0):
        col1.info("No counts identified for five prime miRNA")
        three_figure = draw_altair_graph(three_p_mirna,"Raw Counts","miRNA")
        col2.altair_chart(three_figure, use_container_width = True)

    elif (five_p_mirna.shape[0] > 0 and three_p_mirna.shape[0] == 0):
        col2.info("No counts identified for three prime miRNA")
        five_figure = draw_altair_graph(five_p_mirna,"Raw Counts", "miRNA")
        col1.altair_chart(five_figure, use_container_width = True)

    else:
        st.error("Connection Timeout, or selected miRNA not present in list!")

def nc_multimap_drawing(nc_queried: list, con: duckdb,tab):
    """
    Draws the ncRNAs
    nc_queried: list -> selected ncRNAs
    con: duckdb -> duckdb connection
    tab: st.tabs -> tab where data should drawn in
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
    col1.altair_chart(nc_figure, use_container_width = True)


    nc_sig = con.execute(f"SELECT * from significance WHERE ncRNA IN {nc_queried}").fetchdf()

    #nc_sig.applymap(lambda x: 'background-color: red' if x > 0.01 else 'background-color: green', subset='padj')
    col2.write("ncRNA Differential Gene Expression Information:")
    col2.dataframe(nc_sig, use_container_width = True)

def draw_altair_graph(data_draw, value, gene_annotation = None, timepoint = "Timepoint"):
    """ Draws the trajectories as altair graph
    Args:
        data_draw: pd.DataFrame -> preselected data of genes
        value: str -> should be the column that holds the normalized data
        gene_annotation: str -> columns that is holding the hue
    """
    selection = alt.selection_multi(fields=[gene_annotation], bind='legend')
    chart = (
                alt.Chart(data_draw)
                .mark_boxplot(size = 30)
                .encode(x=alt.X(timepoint),
                        y=alt.Y(value, scale=alt.Scale(domain=[data_draw[value].min()-1, data_draw[value].max()+1])),
                        color=gene_annotation,
                        opacity=alt.condition(selection, alt.value(1), alt.value(0.2)))
                .interactive()
                .properties()
                .add_selection(selection)
                )

    chart_line = (
            alt.Chart(data_draw)
            .mark_line(interpolate = "natural")
            .encode(x=alt.X(timepoint,
                            scale=alt.Scale(padding=1)),
                            y=alt.Y(f"mean({value})",scale=alt.Scale(domain=[data_draw[value].min()-1, data_draw[value].max()+1])),
                            color=gene_annotation,
                            opacity=alt.condition(selection, alt.value(1), alt.value(0.2)),
                            size=alt.condition(~selection, alt.value(1), alt.value(3)))
            .interactive()
            .properties()
            .add_selection(selection))

    final = chart_line + chart
    final  = final.configure_legend(padding=2,
                                    cornerRadius=5,
                                    orient='bottom').configure_axis(grid = False,
                                                                      labelFontSize = 13).configure_view(strokeOpacity = 0)
    return final

def lnc_query_genes(genes_liste, con, tab):
    """ Queries lncRNA genes from the duckdb database
    Args:
        genes_liste: list -> list of selected lncRNAs
        con: duckdb -> duckdb database
        tab: st.tab -> streamlit tab
    """
    columns = [9*[i] for i in ["Day00", "Day05", "Day09", "Day16", "Day26", "Day36"]]
    columns = [t for i in columns for t in i]
    columns = ["gene_name"] + columns

    selected_lncs = con.execute(f"Select * from lnc_counts WHERE external_gene_name IN {tuple(genes_liste)}").fetchdf()
    selected_lncs.columns = columns
    melted_selected_lncs = pd.melt(selected_lncs, id_vars="gene_name")
    melted_selected_lncs.columns = ["gene_name","Timepoint", "vsd counts"]
    fig = draw_altair_graph(melted_selected_lncs, "vsd counts", "gene_name" )
    tab.write("---")
    tab.altair_chart(fig, use_container_width = True)


def make_heatmap(df_curves,  gene_name, value, timepoint, title = None):
    """_summary_

    Args:
        df_curves (pd.DataFrame): Variance stabilized counts of mRNAs
        tab (st.tabs): Streamlit Tab to draw in

    Returns:
        _type_: _description_
    """
    return (
        alt.Chart(df_curves)
        .mark_rect()
        .encode(
            alt.X(timepoint),
            alt.Y(gene_name),
            alt.Color(value, scale=alt.Scale(scheme='greenblue')),
        )
        .properties(title=title)
    )

@st.cache_resource
def execute_fill_list(_con):
    """_summary_: This function should retrieve the list for the multiselect box
    This should be hold in cache since these list will not change

    Args:
        _con (DuckDB): Database connection to duckdb

    Returns:
        tuple(list): A tuple holding the lists (names) for the different ncRNA species
    """
    draw = sorted(list(
        {
            i
            for i in sorted(
                _con.execute("SELECT gene_name from all_counts").fetchnumpy()[
                    "gene_name"
                ]
            )
            if "5" not in i
        }
    ))
    mirna = _con.execute("SELECT gene_name from fivep_counts").fetchnumpy()["gene_name"].tolist() + _con.execute("SELECT gene_name from threep_counts").fetchnumpy()["gene_name"].tolist()
    nc = _con.execute("SELECT gene_name from ncrna_counts").fetchnumpy()["gene_name"].tolist()
    lnc = _con.execute("Select external_gene_name from lnc_counts").fetchdf()
    return (draw, mirna, nc, lnc)

if __name__ == "__main__":
    trajectory_start()