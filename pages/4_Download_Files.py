import streamlit as st
import pandas as pd
import duckdb


def download_files(database):
    """_summary_: Streamlit page to provide download links for the
    NOCICEPTRA database and the NOCICEPTRA database schema
    Args:
        database (DuckDBConnection): _description_
    Returns:
        None
    """
    st.header("Download the NOCICEPTRA database")
    tables = database.execute("SHOW TABLES").fetchdf()


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


if __name__ == "__main__":
    database = load_data()
    download_files(database)
