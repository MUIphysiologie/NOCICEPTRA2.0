# NOCICEPTRA_Tool

# Beta Version:

Tool in Python that should visualize the Analysis of miRNAs, mRNAs, ncRNAs.
<br>
<img src="Images/Readme.png" alt="StartingPage">
<br>

## Details:
Data is provided using the DuckDB OLAP Database.
Tool was updated to use interactive libraries instead of Matplotlib such as Altair and Plotly
In addition also further updates were made on the underlying data structure. Now all data is stored in the duckdb database. Which can be downloaded here on github.


<p> Tool can be started using the following command:
<code>
streamlit run Main_page.py
</code>


## Structure of the app
Currently there exist three different sections, including a Trajectory Analysis, WGNCA analysis and Network analysis module:

<ul>
  <li> Gene, ncRNA, miRNA Trajectory Analysis  <br>
  <img src="Images/Trajectories.png" alt="StartingPage">
  <br>
  </li>
  <li> WGCNA Analysis </li>
  <li> Network Analysis </li>

</ul>

The app is also available at nociceptra.streamlitapps.io, but might lack Memory if parallel access is saturating the memory.

## Contributor to the Applicatin

Maximilian Zeidler PhD

