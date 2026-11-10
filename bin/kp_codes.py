PLASMID_CODE_SYSTEM = "https://kp-pipeline-docs.readthedocs.io/fhir/CodeSystem/kp-plasmid"

PLASMID_ANALYSIS_CODE = "plasmid-analysis"
PLASMID_CODE = "plasmid"
GENOMIC_LOCATION_CODE = "genomic-location"

PLASMID_CODE_DISPLAYS = {
    PLASMID_ANALYSIS_CODE: "Plasmid analysis",
    PLASMID_CODE: "Plasmid identified",
    GENOMIC_LOCATION_CODE: "Genomic location",
}

PLASMID_OBSERVATION_CODES = [PLASMID_ANALYSIS_CODE, PLASMID_CODE]

# Placeholder until a kemkes SP code is issued
QC_CODE_SYSTEM = "https://kp-pipeline-docs.readthedocs.io/fhir/CodeSystem/kp-qc"
SEQUENCING_QC_CODE = "sequencing-qc"
SEQUENCING_QC_DISPLAY = "Sequencing quality control"

SEQUENCE_ONTOLOGY_SYSTEM = "http://www.sequenceontology.org"

LOCATION_CODINGS = {
    "plasmid": {"system": SEQUENCE_ONTOLOGY_SYSTEM, "code": "SO:0000155", "display": "plasmid"},
    "chromosome": {"system": SEQUENCE_ONTOLOGY_SYSTEM, "code": "SO:0000340", "display": "chromosome"},
}
