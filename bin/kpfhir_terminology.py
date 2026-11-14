
LOINC = "http://loinc.org"
SNOMED = "http://snomed.info/sct"
UCUM = "http://unitsofmeasure.org"
HGVS = "http://varnomen.hgvs.org"
HGNC_V1 = "http://www.genenames.org/geneId"
SEQUENCE_ONTOLOGY = "http://www.sequenceontology.org"
KAPTIVE_CAPSULE_V1 = "http://kaptive.holtlab.net/capsule"
KAPTIVE_O_ANTIGEN_V1 = "http://kaptive.holtlab.net/o-antigen"
KAPTIVE_MLST_V1 = "http://kaptive.holtlab.net/mlst"
CGMLST_NCS = "https://www.cgmlst.org/ncs"

CANONICAL_BASE = "http://terminology.spheres.id"
CS_KP_PLASMID = f"{CANONICAL_BASE}/CodeSystem/kp-plasmid"
CS_OBSERVATION_CODE = f"{CANONICAL_BASE}/CodeSystem/kpfhir-observation-code"
CS_KLEBSIELLA_ST = f"{CANONICAL_BASE}/CodeSystem/klebsiella-mlst-st"
CS_KAPTIVE_K = f"{CANONICAL_BASE}/CodeSystem/kaptive-k-locus"
CS_KAPTIVE_O = f"{CANONICAL_BASE}/CodeSystem/kaptive-o-locus"
CS_PIPELINE_PROPERTY = f"{CANONICAL_BASE}/CodeSystem/pipeline-property"

PLASMID_ANALYSIS_CODE = "plasmid-analysis"
PLASMID_CODE = "plasmid"
GENOMIC_LOCATION_CODE = "genomic-location"
PLASMID_CODE_DISPLAYS = {
    PLASMID_ANALYSIS_CODE: "Plasmid analysis",
    PLASMID_CODE: "Plasmid identified",
    GENOMIC_LOCATION_CODE: "Genomic location",
}

SEQUENCING_QC_CODE = "sequencing-qc"
SEQUENCING_QC_DISPLAY = "Sequencing quality control"

OBSERVATION_CODE_DISPLAYS = {
    SEQUENCING_QC_CODE: SEQUENCING_QC_DISPLAY,
    "capsule-type": "Klebsiella pneumoniae capsular type by genomic analysis",
    "virulence-score": "Klebsiella pneumoniae virulence score by genomic analysis",
    "cgmlst-profile": "Core genome MLST profile by sequencing",
    "cgmlst-loci-called": "Core genome MLST loci called",
    "cgmlst-loci-total": "Core genome MLST loci in scheme",
    "cgmlst-call-rate": "Core genome MLST call rate",
}

LOCATION_CODINGS = {
    "plasmid": {"system": SEQUENCE_ONTOLOGY, "code": "SO:0000155", "display": "plasmid"},
    "chromosome": {"system": SEQUENCE_ONTOLOGY, "code": "SO:0000340", "display": "chromosome"},
}

DRUG_CLASS_ALIASES = {
    'agly': 'aminoglycoside', 'flq': 'fluoroquinolone', 'sul': 'sulfonamide',
    'tmt': 'trimethoprim', 'col': 'colistin', 'tet': 'tetracycline',
}

COLISTIN_MUTATION_GENES = ['Colistin', 'MgrB', 'PhoP', 'PhoQ', 'PmrA', 'PmrB', 'CrrA', 'CrrB']

_BETA_LACTAM = {'code': '765422000', 'display': 'Product containing beta-lactam (product)'}
_AMINOGLYCOSIDE = {'code': '324116004', 'display': 'Product containing aminoglycoside (product)'}
_FLUOROQUINOLONE = {'code': '1010205001',
                    'display': 'Medicinal product containing fluoroquinolone and acting as antibacterial agent'}
_MACROLIDE = {'code': '763878009', 'display': 'Product containing macrolide (product)'}
_TETRACYCLINE = {'code': '66261008', 'display': 'Product containing tetracycline (medicinal product)'}
_SULFONAMIDE = {'code': '763875007', 'display': 'Product containing sulfonamide (product)'}
_TRIMETHOPRIM = {'code': '32792001', 'display': 'Product containing trimethoprim (medicinal product)'}
_CHLORAMPHENICOL = {'code': '57191001', 'display': 'Product containing chloramphenicol (medicinal product)'}
_COLISTIN = {'code': '73074003', 'display': 'Product containing colistin (medicinal product)'}
_FOSFOMYCIN = {'code': '387065000', 'display': 'Product containing fosfomycin (medicinal product)'}
_RIFAMPICIN = {'code': '77891000', 'display': 'Product containing rifampicin (medicinal product)'}

DRUG_CLASS_SNOMED = {
    'beta_lactam': _BETA_LACTAM, 'bla': _BETA_LACTAM,
    'aminoglycoside': _AMINOGLYCOSIDE, 'agly': _AMINOGLYCOSIDE,
    'fluoroquinolone': _FLUOROQUINOLONE, 'flq': _FLUOROQUINOLONE,
    'macrolide': _MACROLIDE, 'mls': _MACROLIDE,
    'tetracycline': _TETRACYCLINE, 'tet': _TETRACYCLINE,
    'sulfonamide': _SULFONAMIDE, 'sul': _SULFONAMIDE,
    'trimethoprim': _TRIMETHOPRIM, 'tmt': _TRIMETHOPRIM, 'dfr': _TRIMETHOPRIM,
    'chloramphenicol': _CHLORAMPHENICOL, 'phenicol': _CHLORAMPHENICOL,
    'colistin': _COLISTIN, 'col': _COLISTIN,
    'carbapenem': {'code': '350134005', 'display': 'Product containing carbapenem (product)'},
    'fosfomycin': _FOSFOMYCIN, 'fos': _FOSFOMYCIN,
    'rifampicin': _RIFAMPICIN, 'rif': _RIFAMPICIN,
}

DRUG_LOINC = {
    "Ampicillin": "18864-9",
    "Amoxicillin-Clavulanate": "18861-5",
    "Piperacillin-Tazobactam": "18970-4",
    "Cefotaxime": "18876-3",
    "Ceftazidime": "18879-7",
    "Cefepime": "18879-7",
    "Meropenem": "18943-1",
    "Ciprofloxacin": "18906-8",
    "Gentamicin": "18928-2",
    "Amikacin": "18860-7",
    "Trimethoprim-Sulfamethoxazole": "18998-5",
    "Colistin": "18912-6",
}

AMINO_ACIDS = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln',
    'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr', '*': 'Ter',
}

GENE_DESCRIPTIONS = {
    'SHV-1': 'narrow-spectrum β-lactamase',
    'CTX-M': 'ESBL family',
    'TEM': 'β-lactamase family',
    'KPC': 'class A carbapenemase',
    'NDM': 'class B metallo-β-lactamase (carbapenemase)',
    'VIM': 'class B metallo-β-lactamase (carbapenemase)',
    'IMP': 'class B metallo-β-lactamase (carbapenemase)',
    'OXA-48': 'class D carbapenemase',
}


def drug_class_snomed(drug_class):
    normalized = drug_class.lower().strip().replace('-', '_').replace(' ', '_')
    return DRUG_CLASS_SNOMED.get(normalized)
