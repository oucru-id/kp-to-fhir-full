
from kpfhir_resistance import conclusion_coding
from kpfhir_terminology import (CS_KAPTIVE_K, CS_KAPTIVE_O, CS_KLEBSIELLA_ST, CS_OBSERVATION_CODE, HGNC_V1,
                                KAPTIVE_CAPSULE_V1, KAPTIVE_MLST_V1, KAPTIVE_O_ANTIGEN_V1, LOINC,
                                OBSERVATION_CODE_DISPLAYS, UCUM)

GENOMICS_SD = "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition"
FHIR_OBSERVATION_SD = "http://hl7.org/fhir/StructureDefinition/Observation"

KEMKES_SD = "https://fhir.kemkes.go.id/r4/StructureDefinition"
KEMKES_SP = "http://terminology.kemkes.go.id/sp"
KEMKES_TAG = {"system": KEMKES_SP, "code": "genomics", "display": "Genomics"}
KEMKES_SYS_IDS = "http://sys-ids.kemkes.go.id"
KEMKES_NIK = "https://fhir.kemkes.go.id/id/nik"

REPORT_CODE_V1 = {"system": LOINC, "code": "81247-9", "display": "Master HL7 genetic variant reporting panel"}
REPORT_CODE_IG3 = {"system": LOINC, "code": "51969-4", "display": "Genetic analysis report"}
REPORT_TEXT = "Klebsiella pneumoniae Genomic Analysis Report"

_KEMKES_OBSERVATION_CODES = {
    "capsule-type": ({"system": KEMKES_SP, "code": "SP000678",
                      "display": "Klebsiella pneumoniae capsular type [Identifier] by genomic analysis"},
                     "Capsule Type"),
    "virulence-score": ({"system": KEMKES_SP, "code": "SP000680",
                         "display": "Klebsiella pneumoniae virulence score [Numeric] by genomic analysis"},
                        "Virulence score"),
    "cgmlst-profile": ({"system": KEMKES_SP, "code": "SP000682",
                        "display": "Core genome MLST [Type] in Isolate by Sequencing"},
                       "Core Genome MLST Profile"),
    "cgmlst-loci-called": ({"system": KEMKES_SP, "code": "SP000690", "display": "Loci called"}, None),
    "cgmlst-loci-total": ({"system": KEMKES_SP, "code": "SP000691", "display": "Total loci in scheme"}, None),
    "cgmlst-call-rate": ({"system": KEMKES_SP, "code": "SP000692", "display": "Call rate"}, None),
}
_OBSERVATION_TEXT = {"capsule-type": "Capsule Type", "virulence-score": "Virulence score",
                     "cgmlst-profile": "Core Genome MLST Profile",
                     "sequencing-qc": OBSERVATION_CODE_DISPLAYS["sequencing-qc"]}


class Profile:

    name = "generic"
    include_nik = False
    kemkes_extensions = False
    specimen_defaults = False

    def __init__(self, identifier_system_base=None):
        self.identifier_system_base = (identifier_system_base or "").rstrip("/") or None

    def meta(self, kind):
        profiles = self._profiles(kind)
        meta = {}
        if profiles:
            meta["profile"] = profiles
        tags = self._tags(kind)
        if tags:
            meta["tag"] = tags
        return meta or None

    def _profiles(self, kind):
        return {
            "variant": [f"{GENOMICS_SD}/variant"],
            "panel": [FHIR_OBSERVATION_SD],
            "typing": [FHIR_OBSERVATION_SD],
            "report": [f"{GENOMICS_SD}/genomic-report"],
        }.get(kind, [])

    def _tags(self, kind):
        return []

    def bundle_meta(self, last_updated):
        return {"lastUpdated": last_updated}

    def observation_code(self, code):
        coding = {"system": CS_OBSERVATION_CODE, "code": code, "display": OBSERVATION_CODE_DISPLAYS[code]}
        concept = {"coding": [coding]}
        text = _OBSERVATION_TEXT.get(code)
        if text:
            concept["text"] = text
        return concept

    def component_code(self, code):
        return {"coding": [{"system": CS_OBSERVATION_CODE, "code": code,
                            "display": OBSERVATION_CODE_DISPLAYS[code]}]}

    def gene_concept(self, name):
        return {"text": name}

    def mutation_concept(self, name):
        return {"text": name}

    def sequence_type_coding(self, sequence_type):
        number = sequence_type.replace('ST', '')
        return {"system": CS_KLEBSIELLA_ST, "code": number, "display": f"Klebsiella pneumoniae MLST {sequence_type}"}

    def capsule_coding(self, code, display):
        return {"system": CS_KAPTIVE_K, "code": code, "display": display}

    def o_antigen_coding(self, o_type):
        return {"system": CS_KAPTIVE_O, "code": o_type, "display": f"Klebsiella pneumoniae {o_type}"}

    def score_quantity(self, value):
        return {"value": value, "unit": "score", "system": UCUM, "code": "{score}"}

    def conclusion_coding(self, resistance_class):
        return conclusion_coding(resistance_class, include_kemkes=False)

    def report_code(self):
        return {"coding": [dict(REPORT_CODE_IG3), dict(REPORT_CODE_V1)], "text": REPORT_TEXT}

    def identifier(self, kind, value, org_id=None):
        if not self.identifier_system_base or value in (None, ""):
            return None
        return {"system": f"{self.identifier_system_base}/{kind}", "value": str(value)}


class SatusehatProfile(Profile):
    name = "satusehat"
    include_nik = True
    kemkes_extensions = True
    specimen_defaults = True

    _KEMKES = {
        "patient": "Patient", "specimen": "Specimen", "organization": "Organization",
        "practitioner": "Practitioner", "practitioner_role": "PractitionerRole",
        "service_request": "ServiceRequest", "encounter": "Encounter",
    }

    def _profiles(self, kind):
        if kind in self._KEMKES:
            return [f"{KEMKES_SD}/{self._KEMKES[kind]}"]
        if kind == "panel":
            return [f"{GENOMICS_SD}/therapeutic-implication"]
        if kind == "report":
            return [f"{GENOMICS_SD}/genomics-report"]
        return super()._profiles(kind)

    def _tags(self, kind):
        if kind in ("variant", "panel", "typing", "report", "device", "provenance"):
            return [dict(KEMKES_TAG)]
        return []

    def bundle_meta(self, last_updated):
        return {"lastUpdated": last_updated, "profile": [f"{KEMKES_SD}/Bundle"]}

    def observation_code(self, code):
        if code in _KEMKES_OBSERVATION_CODES:
            coding, text = _KEMKES_OBSERVATION_CODES[code]
            concept = {"coding": [dict(coding)]}
            if text:
                concept["text"] = text
            return concept
        return super().observation_code(code)

    def component_code(self, code):
        if code in _KEMKES_OBSERVATION_CODES:
            return {"coding": [dict(_KEMKES_OBSERVATION_CODES[code][0])]}
        return super().component_code(code)

    def gene_concept(self, name):
        return {"coding": [{"system": HGNC_V1, "code": name, "display": name}], "text": name}

    def mutation_concept(self, name):
        return {"coding": [{"system": LOINC, "code": name, "display": name}], "text": name}

    def sequence_type_coding(self, sequence_type):
        return {"system": KAPTIVE_MLST_V1, "code": sequence_type.replace('ST', ''),
                "display": f"Klebsiella pneumoniae MLST {sequence_type}"}

    def capsule_coding(self, code, display):
        return {"system": KAPTIVE_CAPSULE_V1, "code": code, "display": display}

    def o_antigen_coding(self, o_type):
        return {"system": KAPTIVE_O_ANTIGEN_V1, "code": o_type, "display": f"Klebsiella pneumoniae {o_type}"}

    def score_quantity(self, value):
        return {"value": value, "unit": "score", "system": UCUM, "code": "score"}

    def conclusion_coding(self, resistance_class):
        return conclusion_coding(resistance_class, include_kemkes=True)

    def report_code(self):
        return {"coding": [dict(REPORT_CODE_V1)], "text": REPORT_TEXT}

    def identifier(self, kind, value, org_id=None):
        if value in (None, ""):
            return None
        if kind == "organization":
            return {"use": "official", "system": f"{KEMKES_SYS_IDS}/organization", "value": str(value)}
        if not org_id:
            return None
        system_kind = {"service_request": "servicerequest"}.get(kind, kind)
        return {"system": f"{KEMKES_SYS_IDS}/{system_kind}/{org_id}", "value": str(value)}


PROFILES = {"generic": Profile, "satusehat": SatusehatProfile}


def get_profile(name, identifier_system_base=None):
    try:
        return PROFILES[name](identifier_system_base=identifier_system_base)
    except KeyError:
        raise SystemExit(f"unknown --fhir_profile {name!r}; choose one of {sorted(PROFILES)}")
