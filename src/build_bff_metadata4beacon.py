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
analyses = {"aligner" : 'dragen,Version="SW: 07.021.624.3.10.8, HW: 07.021.624"', 
            "analysisDate" : "2024-02-13", 
            "pipelineName" : "dragen-variant-caller-for-AIC", 
            "pipelineRef" : "(Blanchet R et al. 2024)", 
            "variantCaller" : 'dragen,Version="SW: 07.021.624.3.10.8, HW: 07.021.624"',
            "idbasename" : "AIC_analyisid",
            "runbasename" : "SRR00"
            }
librarySource = {"id": "GENEPIO:0001966", "label": "genomic Source"}
runs = {'libraryLayout' : 'PAIRED', 
        'librarySelection' : 'RANDOM', 
        'librarySource' : librarySource, 
        'platform' : 'unknown', 
        'platformModel' : 'unknown'}
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