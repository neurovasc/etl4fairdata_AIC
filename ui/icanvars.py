import streamlit as st
import time
import requests
import utils as ut

#### Page configuration
st.set_page_config(
    page_title="Aggregated variants Knowledge Graph demo",
    layout="centered",
    initial_sidebar_state="collapsed",
)

st.info(
    "Check [related work on Semantic Beacons]() and the [documentation](https://github.com/SemanticBeacon/SemanticBeacon) of the project"
)

st.title('Aggregated variants from the iCAN cohort')
st.markdown(
    '<div style="text-align: justify;">This page demonstrates the integration of aggregated genomic variants \
    from the iCAN cohort into a knowledge graph. \
    The iCAN cohort is the result of the French nationwide collaborative project to study and understand \
    the physiopathology of intracranial aneurysms. \
    See <a href="https://doi.org/10.1093/neuros/nyw135">Bourcier et al. 2017</a> for more information \
    about inclusion and exclusion criteria. <br><br>\
    The knowledge graph is build following the data model described below. It contains genomic informations \
    obtained from the exome sequencing of patients carrying intracranial aneurysms that also have a detailed\
    clinical profile documented in the medical database GAIA. \
    The knowledge graph contains aggregated variants, meaning there is no individual level information. \
    What is represented is the presence of a variant in the cohort and the allele count and frequency in \
    respect to phenotypes of interest in the light of intracranial aneurysms analysis. \
    </div>', unsafe_allow_html=True
)

st.subheader("Summary")
st.markdown(
    """
    1. Knowledge graph schema and ontology choice
    2. Execute prewritten queries from our catalog
"""
)

# KNOWLEDGE GRAPH STRUCTURE DESCRIPTION

st.subheader("1. Knowledge graph schema and ontology choice")
st.markdown(
    """
    We use several domain expert ontologies to describe different aspects of the genomic variant \
    (location, alleles, identifiers, zygosity). The genomic variant is liked to an Observation that \
    is the combination of a phenotype, a zygosity, a frequency and/or a count value. \
    The variant is also associated with the Cohort. The Cohort is yet to be described in ontological terms. <br>\
    """, 
    unsafe_allow_html=True
)
st.image("images/data-model-aggregated-variantion-in-cohort_v1_2025-02-21_nocolors.png", caption="Aggregated variant data model")


st.subheader("2. Execute prewritten queries from our catalog")

# Question by default
question = ut.querycatalog[0]["descriptor"]
thequery = ut.querycatalog[0]["query"]
# Select question
questions = [item["descriptor"] for item in ut.querycatalog]
select_question = st.selectbox("", questions)
# Find the query according to the question
thequery = next(item for item in ut.querycatalog if item["descriptor"] == select_question)

# Tabs to display the query, and execute it
tabQuery, tabResp = st.tabs(["Sparql query", "Execute query"])

with tabQuery:
    st.info(
        thequery["descriptor"]
    )
    st.code(thequery["query"], language="sparql", line_numbers=True)

with tabResp:
    st.info(
        thequery["descriptor"]
    )
    #st.code(thequery["query"], language="sparql", line_numbers=True)
    with st.spinner("This may take a while."):
        start_time = time.time()  # Record start time
        comres_stdout, comres_stderr = ut.executequery(thequery["query"])
        end_time = time.time()  # Record end time
        execution_time = end_time - start_time  # Calculate duration
        #
        st.success(f"query complete -> (Time taken: {execution_time:.2f} seconds)")
        #
        prettydf = ut.resultlayout(comres_stdout, comres_stderr)
        st.dataframe(prettydf)
        #
        st.code(comres_stdout, language="sparql", line_numbers=True)