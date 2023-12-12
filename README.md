# NOCICEPTRA_Tool

# Beta Version:

Exploratory Tool that should visualize the Analysis of miRNAs, mRNAs and ncRNAs.
<br>
<img src="Images/Readme.png" alt="StartingPage">
<br>

## Details:
Data ETL was performed using the DuckDB OLAP Database.
The original tool was updated to use interactive libraries (e.g Plotly and Altair) instead of Matplotlib.
In addition also further updates were made on the underlying data structure. Now all data is stored in a comman duckdb database which only allows for read operations. The database can be downloaded here on Github.


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

The app is also available at nociceptra.streamlitapps.io, but might lack Memory if parallel access is saturating the memory, since it is only using the Free Tier here.

## Contributors

Maximilian Zeidler PhD

