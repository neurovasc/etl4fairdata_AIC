from rdflib import Namespace
AGGRVAR = Namespace('https://ican.univ-nantes.io/variants-kg-schema/')
ICAN = Namespace('http://ican.ressource.org/')
FALDO = Namespace('http://biohackathon.org/resource/faldo#')
GENO = Namespace('http://purl.obolibrary.org/obo/GENO_')
HPO = Namespace('http://purl.obolibrary.org/obo/HP_')
LINKML = Namespace('https://w3id.org/linkml/')
MONDO = Namespace('http://purl.obolibrary.org/obo/MONDO_')
NCIT = Namespace('http://purl.obolibrary.org/obo/NCIT_')
OWL = Namespace('http://www.w3.org/2002/07/owl#')
RDF = Namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#')
RDFS = Namespace('http://www.w3.org/2000/01/rdf-schema#')
SIO = Namespace('http://semanticscience.org/resource/SIO_')
SKOS = Namespace('http://www.w3.org/2004/02/skos/core#')
SO = Namespace('http://purl.obolibrary.org/obo/SO_')
XSD = Namespace('http://www.w3.org/2001/XMLSchema#')

namespace = {'faldo' : FALDO, 'aggrvar' : AGGRVAR, 'aic' : ICAN,
             'geno' : GENO, 'hpo' : HPO, 'linkml' : LINKML, 'mondo' : MONDO,
             'ncit' : NCIT, 'owl' : OWL, 'rdf' : RDF, 'rdfs' : RDFS,
             'sio' : SIO, 'skos' : SKOS, 'so' : SO, 'xsd' : XSD}