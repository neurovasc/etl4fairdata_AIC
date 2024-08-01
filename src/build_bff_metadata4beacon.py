# varibales for BeaconFriendlyFormat metadata
# classes for BeaconFriendlyFormat categories

# Metadata
analyses = {"aligner" : 'dragen,Version="SW: 07.021.624.3.10.8, HW: 07.021.624"', 
            "analysisDate" : "2024-02-13", 
            "pipelineName" : "dragen-variant-caller-for-AIC", 
            "pipelineRef" : "(Blanchet R et al. 2024)", 
            "variantCaller" : 'dragen,Version="SW: 07.021.624.3.10.8, HW: 07.021.624"',
            "idbasename" : "AIC_analyisid_",
            "runbasename" : "SRR00"
            }
# Classes
# Class
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
        return (f"Analysis(aligner={self.aligner},  \
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