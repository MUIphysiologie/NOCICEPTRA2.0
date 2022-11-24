# Import all the modules necessary for the Tool
import streamlit as st 
import duckdb


st.set_page_config(page_title = "NOCICEPTRA",layout="wide") 
# set the tab name #draw interactive visualization using holoviews
st.session_state.update(st.session_state)
@st.experimental_singleton
def load_data():
    """
    Defines the database connection to the NOCICEPTRA duckdb database
    """
    try:
        con = duckdb.connect(database = "./Data/nociceptra.duckdb", read_only = True)
        st.session_state.connection = con
    except Exception as e:
        print(f"Error: {e}")

def main():
    """
    Landing Page of the Nociceptra Trajectory Tool
    """

    #set the logo
    st.sidebar.image("./Images/Logo_iPSC.png", output_format = "png", width = 200)
    load_data()
    start_page()


        
    
def start_page():
    """ layout of the Webpage for the susceptibility windows"""
    st.image("./Images/Logo_iPSC.png", output_format = "png", use_column_width = True)
    st.image("./Images/Web_landing_page.png", output_format = "png", use_column_width = True)
    st.write(r"""***Articles:***""")
    st.markdown(r""" ***NOCICEPTRA: Gene and microRNA signatures and their trajectories characterizing human iPSC-derived nociceptor maturation.***  <br>
    Zeidler M, Kummer K, Schoepf C, Theodora Kalpachidou, Kern G, Cader Z, Kress M, Advanced Science (2021)""", unsafe_allow_html=True)
    st.markdown(r""" ***Selected Ionotropic Receptors and Voltage-Gated Ion Channels: More Functional Competence for Human Induced Pluripotent Stem Cell (iPSC)-Derived Nociceptors.***
    by Clemens L. Schoepf, Maximilian Zeidler, Lisa Spiecker, Georg Kern, Judith Lechner, Kai K. Kummer and Michaela Kress, (2020) Brain Sciences""", unsafe_allow_html=True)
    
  
    # add 
    
if __name__ == "__main__":
    main() 