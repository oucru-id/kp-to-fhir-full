
import copy
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from kpfhir_terminology import CS_KP_PLASMID, CS_OBSERVATION_CODE  # noqa: E402

V1_PLASMID_CS = "https://kp-pipeline-docs.readthedocs.io/fhir/CodeSystem/kp-plasmid"
V1_QC_CS = "https://kp-pipeline-docs.readthedocs.io/fhir/CodeSystem/kp-qc"
NIK_SYSTEM = "https://fhir.kemkes.go.id/id/nik"
NEW_RESOURCE_TYPES = {"Encounter", "Device", "Provenance"}
TIMESTAMP = "<TIMESTAMP>"


def _replace_strings(node, old, new):
    if isinstance(node, dict):
        return {k: _replace_strings(v, old, new) for k, v in node.items()}
    if isinstance(node, list):
        return [_replace_strings(v, old, new) for v in node]
    if isinstance(node, str):
        return node.replace(old, new)
    return node


def _patient(resource):
    names = []
    for name in resource.get("name", []):
        name = {k: v for k, v in name.items() if v != "Unknown"}
        if "given" in name:
            name["given"] = [g for g in name["given"] if g != "Unknown"] or None
            if name["given"] is None:
                del name["given"]
        if set(name) - {"use"}:
            names.append(name)
    if names:
        resource["name"] = names
    else:
        resource.pop("name", None)
    if resource.get("birthDate") == "Unknown":
        del resource["birthDate"]
    identifiers = []
    for identifier in resource.get("identifier", []):
        if identifier.get("system") == NIK_SYSTEM:
            value = str(identifier.get("value", ""))
            if not value.isdigit():
                continue
            identifier["value"] = value.zfill(16)
        identifiers.append(identifier)
    resource["identifier"] = identifiers


def expected_v2(v1_bundle):
    bundle = copy.deepcopy(v1_bundle)
    bundle = _replace_strings(bundle, "https://varnomen.hgvs.org", "http://varnomen.hgvs.org")
    bundle = _replace_strings(bundle, V1_PLASMID_CS, CS_KP_PLASMID)
    bundle = _replace_strings(bundle, V1_QC_CS, CS_OBSERVATION_CODE)

    for resource in (e["resource"] for e in bundle["entry"]):
        rtype = resource["resourceType"]
        if rtype == "DiagnosticReport":
            sid = resource["id"][: -len("-genomic-report")]
            resource["text"] = {"status": "generated",
                                "div": ('<div xmlns="http://www.w3.org/1999/xhtml">KP genomic analysis report '
                                        f'for sample {sid}. {resource["conclusion"]}</div>')}
            resource.get("subject", {}).pop("display", None)
            for spec in resource.get("specimen", []):
                spec.pop("display", None)
            for form in resource.get("presentedForm", []):
                form["data"] = re.sub(r"<p><strong>Patient:</strong> .*?</p>\n", "", form["data"])
        elif rtype == "Specimen":
            sid = resource["id"][: -len("-specimen")]
            resource["subject"].pop("display", None)
            resource["note"] = [{"text": f"Collected clinical specimen ({sid}) for Klebsiella pneumoniae "
                                         "genomic testing"}]
            if resource.get("receivedTime") == TIMESTAMP:
                del resource["receivedTime"]
            if resource.get("collection", {}).get("collectedDateTime") == TIMESTAMP:
                del resource["collection"]["collectedDateTime"]
        elif rtype == "ServiceRequest":
            resource["subject"].pop("display", None)
        elif rtype == "Patient":
            _patient(resource)
    return bundle


def without_new_resources(bundle):
    bundle = copy.deepcopy(bundle)
    bundle["entry"] = [e for e in bundle["entry"] if e["resource"]["resourceType"] not in NEW_RESOURCE_TYPES]
    return bundle
