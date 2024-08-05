# varibales for BeaconFriendlyFormat metadata
# classes for BeaconFriendlyFormat categories

# Analyses and Runs: use biosampleid list
# Cohorts and Datasets : 
# Individuals and Biosamples
# GenomicVariations

import pandas as pd

#
# Metadata
#
# what is the metadata that is true for all analyses?
analyses = {"aligner" : 'dragen,Version="SW: 07.021.624.3.10.8, HW: 07.021.624"', 
            "analysisDate" : "2024-02-13", 
            "pipelineName" : "dragen-variant-caller-for-AIC", 
            "pipelineRef" : "(Blanchet R et al. 2024)", 
            "variantCaller" : 'dragen,Version="SW: 07.021.624.3.10.8, HW: 07.021.624"',
            "idbasename" : "AIC_analyisid",
            "runbasename" : "SRR00"
            }
# what is the metadata that is true for all runs?
librarySource = {"id": "GENEPIO:0001966", "label": "genomic Source"}
runs = {'libraryLayout' : 'PAIRED', 
        'librarySelection' : 'RANDOM', 
        'librarySource' : librarySource, 
        'platform' : 'unknown', 
        'platformModel' : 'unknown'}
# what is the metadata that is true for the cohort, no matters its subsampling?
# https://docs.genomebeacons.org/schemas-md/cohorts_defaultSchema/
# DOI:10.1093/neuros/nyw135
cohortDesign = {'id': 'OMIABIS:0001020', 'label': 'prospective cohort study'}
inclusionCriteria = {
                    "genders" : [{"id": "NCIT:C16576","label": "female"},
                                {"id": "NCIT:C20197","label": "male"}],

                    "ageRange" : [{"end": {"iso8601duration": "P120Y"}, 
                                   "start": {"iso8601duration": "P0Y"}}],

                    "locations": [{"id": "SNOMED:223666001", "label": "France"}]
                    }
# if needed, we can expand inclusion criteria, add exclusion criteria, and refine geographic
# origin. See publication (Bourcier R et al. 2017)
cohorts = {'cohortDesign' : cohortDesign,
           'cohortId' : 'ICAN', 
           'cohortType' : 'user-defined', 
           'cohortName' : 'The French ICAN project',
           'id' : 'ICAN', 
           'inclusionCriteria' : inclusionCriteria
           }
# to build: cohortSize int, collectionEvents [], ids {biosamples : [], individuals: []}
#
# Functions 
#
def load_samples(samples):
    ''' Load samples from a file containing phenotypes
    '''
    df = pd.read_csv(samples, sep=',')
    samples = df['N°ADN IRT 1'].tolist()

    return samples
#
# Classes
#
class Analysis:
    #
    def __init__(self, aligner, analysisDate, pipelineName, pipelineRef, variantCaller,
                biosampleId, id, individualId, runId):
        self.aligner = aligner
        self.analysisDate = analysisDate
        self.biosampleId = biosampleId
        self.id = id
        self.individualId = individualId
        self.pipelineName = pipelineName
        self.pipelineRef = pipelineRef
        self.runId = runId
        self.variantCaller = variantCaller
    #
    def __str_(self):
        return (f"Analysis(aligner={self.aligner}, \
                analysisDate={self.analysisDate}, \
                biosampleId={self.biosampleId}, \
                id={self.id}, \
                individualId={self.individualId}, \
                pipelineName={self.pipelineName}, \
                pipelineRef={self.pipelineRef}, \
                runId={self.runId}, \
                variantCaller={self.variantCaller})")
    #
    def to_dict(self):
        return {
            'aligner': self.aligner,
            'analysisDate': self.analysisDate,
            'biosampleId': self.biosampleId,
            'id': self.id,
            'individualId': self.individualId,
            'pipelineName': self.pipelineName,
            'pipelineRef': self.pipelineRef,
            'runId': self.runId,
            'variantCaller': self.variantCaller
        }
#
class Run:
    #
    def __init__(self, libraryLayout, librarySelection, librarySource, LibraryStrategy,
                 platform, platformModel,
                 biosampleId, id, individualId, runDate):
        self.libraryLayout = libraryLayout
        self.librarySelection = librarySelection
        self.librarySource = librarySource
        self.LibraryStrategy = LibraryStrategy
        self.platform = platform
        self.platformModel = platformModel
        self.biosampleId = biosampleId
        self.id = id
        self.individualId = individualId
        self.runDate = runDate
        
    #
    def __str_(self):
        return (f"Run(libraryLayout={self.libraryLayout},  \
            librarySelection={self.librarySelection}, \
            librarySource={self.librarySource}, \
            libraryStrategy={self.LibraryStrategy}, \
            platform={self.platform}, \
            platformModel={self.platformModel}, \
            biosampleId={self.biosampleId}, \
            id={self.id}, \
            individualId={self.individualId}, \
            runDate={self.runDate})")
    #
    def to_dict(self):
        return {
            'libraryLayout': self.libraryLayout,
            'librarySelection': self.librarySelection,
            'librarySource': self.librarySource,
            'libraryStrategy': self.LibraryStrategy,
            'platform': self.platform,
            'platformModel': self.platformModel,
            'biosampleId': self.biosampleId,
            'id': self.id,
            'individualId': self.individualId,
            'runDate': self.runDate
        }
#
class Cohort:
    #
    def __init__(self, cohortDesign, cohortId, cohortType, cohortName, id, inclusionCriteria,
                 cohortSize, collectionEvents, ids):
        self.cohortDesign = cohortDesign
        self.cohortId = cohortId
        self.cohortName = cohortName
        self.cohortType = cohortType
        self.id = id
        self.inclusionCriteria = inclusionCriteria
        self.cohortSize = cohortSize
        self.collectionEvents = collectionEvents
        self.ids = ids
    #
    def __str__(self):
        return (f"Cohort(cohortDesign={self.cohortDesign}, \
                cohortId={self.cohortId}, \
                cohortName={self.cohortName}, \
                cohortType={self.cohortType}, \
                id={self.id}, \
                inclusionCriteria={self.inclusionCriteria}, \
                cohortSize={self.cohortSize}, \
                collectionEvents={self.collectionEvents}, \
                ids={self.ids})")
    #
    def to_dict(self):
        return {
            'cohortDesign' : self.cohortDesign,
            'cohortId': self.cohortId,
            'cohortName': self.cohortName,
            'cohortType': self.cohortType,
            'id': self.id,
            'inclusionCriteria': self.inclusionCriteria,
            'cohortSize': self.cohortSize,
            'collectionEvents': self.collectionEvents,
            'ids': self.ids
        }
    
        
