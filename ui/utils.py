import streamlit as st
import subprocess
import json 
import pandas as pd

query1 = """ s-query --service http://localhost:3030/ican '
PREFIX so: <http://purl.obolibrary.org/obo/SO_>

SELECT COUNT(?variant)
WHERE {
  ?variant a so:0001059 .
}
'
"""

query2 = """ s-query --service http://localhost:3030/ican '
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>

SELECT COUNT(?variant)
WHERE {
  ?variant a so:0001059 .
  ?variant faldo:location ?region .
  ?region faldo:reference ?chr .
  ?chr rdf:label "6" .
  #FILTER (STR(?label) = "6") 
}
'
"""

query3 = """ s-query --service http://localhost:3030/ican '
PREFIX xs: <http://www.w3.org/2001/XMLSchema#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>
PREFIX sio: <http://semanticscience.org/resource/SIO_>  
prefix geno: <http://purl.obolibrary.org/obo/GENO_> 
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

SELECT ?populationlabel ?dellen ?chromosome ?start ?end ?refseq ?altseq ?variantID
WHERE {
  ?variant a so:0001059 .
  ?variant sio:000671 ?videntifier .
  ?videntifier sio:000300 ?variantID .
  # location
  ?variant faldo:location ?region .
  ?region faldo:reference ?chr .
  ?chr rdf:label ?chromosome .
  #?chr rdf:label "6" .
  ?region faldo:begin ?beginEP .
  ?region faldo:end ?endEP .
  ?beginEP faldo:position ?start .
  ?endEP faldo:position ?end .
  # alleles
  ?variant geno:0000385 ?ref .
  ?variant geno:0000382 ?alt .
  ?ref sio:000300 ?refseq .
  ?alt sio:000300 ?altseq .
  # observation
  ?variant sio:001403 ?observation .
  ?observation sio:001403 ?subpopulation .
  ?observation geno:000608 geno:0000135 . 
  ?observation sio:000216 ?count .
  ?count sio:000300 ?countvalue .
  ?subpopulation rdf:label ?populationlabel .

  BIND((STRLEN(?refseq) - STRLEN(?altseq)) AS ?dellen)
  FILTER (?dellen > 20 && ?countvalue>5)
  FILTER (?populationlabel = "earlyonset")
}
'
"""

query4 = """ s-query --service http://localhost:3030/ican '
PREFIX xs: <http://www.w3.org/2001/XMLSchema#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>
PREFIX sio: <http://semanticscience.org/resource/SIO_>  
prefix geno: <http://purl.obolibrary.org/obo/GENO_> 
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

SELECT ?variantID ?countvalue ?subpopulation
WHERE {
  ?variant a so:0001059 .
  ?variant sio:000671 ?videntifier .
  ?videntifier sio:000300 ?variantID .
  # location
  ?variant faldo:location ?region .
  ?region faldo:reference ?chr .
  ?chr rdf:label ?chromosome .
  ?region faldo:begin ?beginEP .
  ?region faldo:end ?endEP .
  ?beginEP faldo:position ?start .
  ?endEP faldo:position ?end .
  # alleles
  ?variant geno:0000385 ?ref .
  ?variant geno:0000382 ?alt .
  ?ref sio:000300 ?refseq .
  ?alt sio:000300 ?altseq .
  # observation
  ?variant sio:001403 ?observation . 
  ?observation sio:000216 ?count .
  ?count sio:000300 ?countvalue .
  ?observation sio:001403 ?subpopulation .
  ?subpopulation rdf:label ?populationlabel .

  BIND((STRLEN(?refseq) - STRLEN(?altseq)) AS ?dellen)
  FILTER (?dellen > 20 && ?countvalue>5)
  FILTER (?populationlabel = "whole")
}'
"""

query5 = """ s-query --service http://localhost:3030/ican '
PREFIX si: <http://sisteminformasi.com/>
PREFIX xs: <http://www.w3.org/2001/XMLSchema#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>
PREFIX sio: <http://semanticscience.org/resource/SIO_>  
prefix geno: <http://purl.obolibrary.org/obo/GENO_> 
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

SELECT ?variantID ?altseq ?populationlabel ?zygosity ?frequencyValue
WHERE {
  # variantID
  ?variant a so:0001059 .
  ?variant sio:000671 ?videntifier .
  ?videntifier sio:000300 ?variantID .
  # alleles
  ?variant geno:0000385 ?ref .
  ?variant geno:0000382 ?alt .
  ?ref sio:000300 ?refseq .
  ?alt sio:000300 ?altseq .
  # location
  ?variant faldo:location ?region .
  ?region faldo:reference ?chr .
  ?chr rdf:label ?chromosome .
  ?region faldo:begin ?beginEP .
  ?region faldo:end ?endEP .
  ?beginEP faldo:position ?v_start .
  ?endEP faldo:position ?v_end .
  # observation
  ?variant sio:001403 ?observation .
  ?observation sio:000900 ?frequency .
  ?frequency sio:000300 ?frequencyValue .
  ?observation sio:000216 ?count .
  ?count sio:000300 ?countValue .
  ?observation sio:001403 ?subpop .
  ?subpop rdf:label ?populationlabel .
  ?observation geno:000608 ?zygosity .
  
  # subpopulation filter
  FILTER(xsd:integer(?countValue) != 0)
  FILTER (?populationlabel = "obese")
  FILTER (?altseq = "GC")
  FILTER (?zygosity = geno:0000137)
 }
LIMIT 10
'
"""

query6 = """ s-query --service http://localhost:3030/ican '
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX si: <http://sisteminformasi.com/>
PREFIX xs: <http://www.w3.org/2001/XMLSchema#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>
PREFIX sio: <http://semanticscience.org/resource/SIO_>  
prefix geno: <http://purl.obolibrary.org/obo/GENO_> 
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
#
PREFIX p: <http://www.wikidata.org/prop/>
PREFIX pq: <http://www.wikidata.org/prop/qualifier/>
PREFIX ps: <http://www.wikidata.org/prop/statement/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
PREFIX GO: <http://purl.obolibrary.org/obo/GO_>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX wdp: <http://www.wikidata.org/prop/>
PREFIX wdpq: <http://www.wikidata.org/prop/qualifier/>
PREFIX wdps: <http://www.wikidata.org/prop/statement/>

SELECT ?variantID ?wikiGene ?geneName ?fullProteinName ?goTermLabel ?countValue ?frequencyValue
WHERE {
    # get proteins involved in angiogenesis from uniprot
    SERVICE <https://sparql.uniprot.org/sparql> {
      ?protein a up:Protein ;
          up:organism taxon:9606 ;      
          up:classifiedWith ?goTerm ;
    	  up:encodedBy ?gene .
      ?gene  skos:prefLabel ?geneName .
      #?protein up:mnemonic "NPAS3_HUMAN" .
      ?protein up:recommendedName ?recommendedName .
  	  ?recommendedName up:fullName ?fullProteinName .
      ?goTerm rdfs:subClassOf* GO:0001525 .
      ?goTerm rdfs:label ?goTermLabel .
    }
    BIND(SUBSTR(STR(?protein), STRLEN(STR(up:)) + 4) AS ?proteinID2) .
	SERVICE <https://query.wikidata.org/sparql> {
    ?wp wdt:P352 ?proteinID2 ;
        wdt:P702 ?wikiGene . 
    ?wikiGene wdp:P644 ?wgss ;
        wdp:P645 ?wgse .
    ?wgss wdps:P644 ?startcoordinate ;
        wdpq:P1057/wdt:P1813 ?wchromosome ;
        wdpq:P659/rdfs:label ?assembly .
    ?wgse wdps:P645 ?endcoordinate ;
        wdpq:P1057/wdt:P1813 ?wchromosome ;
        wdpq:P659/rdfs:label ?assembly .
    FILTER(lang(?assembly) = "en")
    FILTER(STR(?assembly) = "genome assembly GRCh38")
  }
  # variantID
  ?variant a so:0001059 .
  ?variant sio:000671 ?videntifier .
  ?videntifier sio:000300 ?variantID .
  # location
  ?variant faldo:location ?region .
  ?region faldo:reference ?chr .
  ?chr rdf:label ?chromosome .
  ?region faldo:begin ?beginEP .
  ?region faldo:end ?endEP .
  ?beginEP faldo:position ?v_start .
  ?endEP faldo:position ?v_end .
  # observation
  ?variant sio:001403 ?observation .
  ?observation sio:000900 ?frequency .
  ?frequency sio:000300 ?frequencyValue .
  ?observation sio:000216 ?count .
  ?count sio:000300 ?countValue .
  ?observation sio:001403 ?subpop .
  ?subpop rdf:label "earlyonset" .
  ?observation geno:000608 geno:0000137 . # heterozygous
  
  # variant filter by frequency in iCAN
  FILTER (xsd:float(?frequencyValue) > 0.01)
  FILTER (xsd:float(?frequencyValue) < 0.10)
  FILTER (xsd:integer(?countValue) > 3)
  FILTER( xsd:string(?chromosome) = xsd:string(?wchromosome) && ((((?v_start >= xsd:integer(?startcoordinate)) && 
            (?v_start <= xsd:integer(?endcoordinate)) )) 
        || ((?v_end >= xsd:integer(?startcoordinate)) && 
        (?v_end <= xsd:integer(?endcoordinate))) ))
 }
LIMIT 10
'
"""

querycatalog = [{"descriptor" : "[local] How many variants are there in the KG?",
                 "query" : query1}, 

                {"descriptor" : "[local] How many variants are there in chromosome 6?",
                "query" : query2}, 

                {"descriptor" : "[local] List the deletions larger than 20 nucleotides in any chromosome that have a count larger than 5 in early-onset subpopulation?",
                "query" : query3},

                {"descriptor" : "[local] List deletions larger than 20 nucleotide occur in more than 5 individuals (any chromosome, any (sub)population)?",
                 "query" : query4}, 

                {"descriptor" : "[local] Retrieve 10 variants with a alternative heterozygous allele 'GC' among individuals with obese phenotype?",
                 "query" : query5},

                {"descriptor" : "[federated] Retrieve 10 variants located in genes coding for proteins involved in angiogenesis* among individuals with early onset phenotype?",
                 "query" : query6} ]

def executequery(command):
    try:
        # Run the command and capture the output
        result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
        # Save output to a JSON file

        #st.write(result.stdout)
        #st.write(result.stderr)
        return result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        #st.write(e.stdout)
        #st.write(e.stderr)
        return e.stdout, e.stderr

def resultlayout(stdout, stderr):
    """
    Parses SPARQL JSON output into a Pandas DataFrame.
    """
    output = stdout
    try:
        data = json.loads(output)  # Convert output string to JSON
        headers = data["head"]["vars"]  # Column names
        rows = data["results"]["bindings"]  # Query results

        # Convert to a list of dictionaries
        formatted_rows = [
            {col: row[col]["value"] if col in row else "" for col in headers} for row in rows
        ]
        return pd.DataFrame(formatted_rows)  # Convert to Pandas DataFrame
    except Exception as e:
        return f"Error parsing SPARQL output: {e}"