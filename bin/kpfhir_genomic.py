#!/usr/bin/env python3

import argparse
import base64
import html
import json
import os
import re
import sys
from collections import defaultdict

from kpfhir_common import (SOURCE_REPOSITORY, assembly_sample_id, bundle_id, fhir_id, full_url, now, now_iso,
                           point_references_at_entries, sanitize_id, today_compact, xhtml_text)
from kpfhir_profiles import get_profile
from kpfhir_resistance import classify_resistance_profile, predict_panel
from kpfhir_terminology import (AMINO_ACIDS, CGMLST_NCS, COLISTIN_MUTATION_GENES, CS_KP_PLASMID,
                                CS_PIPELINE_PROPERTY, GENE_DESCRIPTIONS, GENOMIC_LOCATION_CODE, HGVS,
                                LOCATION_CODINGS, LOINC, PLASMID_ANALYSIS_CODE, PLASMID_CODE,
                                PLASMID_CODE_DISPLAYS, SNOMED, UCUM, DRUG_LOINC, drug_class_snomed)

XHTML = "http://www.w3.org/1999/xhtml"
CATEGORY_LAB = {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory", "display": "Laboratory"}]}
CATEGORY_GE = {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                           "code": "GE", "display": "Genetics"}]}
REPORT_TITLE = "Klebsiella pneumoniae Genomic Analysis Report"


class InputError(Exception):
    pass


def div(inner):
    return f"<div xmlns=\"{XHTML}\">{inner}</div>"


def with_meta(resource, meta):
    if not meta:
        return resource
    out = {}
    for key, value in resource.items():
        out[key] = value
        if key == "id":
            out["meta"] = meta
    return out


def log(message):
    print(f"DEBUG: {message}", file=sys.stderr)


class SampleContext:
    def __init__(self, sample_id, fid, profile, performer, timestamp):
        self.sample_id = sample_id
        self.fid = fid
        self.profile = profile
        self.performer = performer
        self.timestamp = timestamp

    def refs(self, specimen=True):
        refs = {"subject": {"reference": f"Patient/{self.fid}-patient"}}
        if specimen:
            refs["specimen"] = {"reference": f"Specimen/{self.fid}-specimen"}
        refs["effectiveDateTime"] = self.timestamp
        if self.performer:
            refs["performer"] = [{"reference": self.performer["reference"]}]
        return refs

    def observation(self, kind, obs_id, body, specimen=True):
        resource = {"resourceType": "Observation", "id": obs_id}
        resource.update(body)
        resource.update(self.refs(specimen))
        return with_meta(resource, self.profile.meta(kind))


def text_component(label, value):
    if isinstance(value, list):
        value = ', '.join(str(v) for v in value)
    if value in [None, '']:
        return None
    return {"code": {"text": label}, "valueString": str(value)}


def virulence_components(typing_data):
    return [{"code": {"text": f"Virulence factor {name}"}, "valueString": value}
            for name, value in typing_data.get('virulence', {}).items()
            if name != 'virulence_score' and value not in ['Unknown', 'None', '']]


def gene_annotation(gene, drug_class, typing_data, intrinsic):
    clean = gene.replace('^', '').strip()
    annotation = f"Resistance gene {gene} - {GENE_DESCRIPTIONS.get(clean, f'{drug_class} resistance gene')}"
    if intrinsic:
        annotation += " (intrinsic chromosomal)"
    sequence_type = typing_data.get('mlst', {}).get('sequence_type', 'Unknown')
    k_type = typing_data.get('capsule', {}).get('k_type', 'Unknown')
    if sequence_type not in ['Unknown', 'Undetermined']:
        annotation += f" - KP MLST: {sequence_type}"
    if k_type not in ['Unknown', 'Capsule null']:
        annotation += f" - Capsule: {k_type}"
    return annotation


def plasmid_clusters_for(location_info):
    clusters = []
    for copy in location_info.get('copies', []):
        cluster_id = copy.get('cluster_id')
        if copy.get('molecule') == 'plasmid' and cluster_id and cluster_id not in clusters:
            clusters.append(cluster_id)
    return clusters


def describe_plasmid(plasmid):
    details = [', '.join(plasmid.get('replicon_types') or []) or 'replicon not typed']
    if plasmid.get('predicted_mobility'):
        details.append(plasmid['predicted_mobility'])
    return f"{plasmid.get('cluster_id')} ({', '.join(details)})"


def describe_location(location_info, plasmids_by_cluster):
    location = location_info.get('location', 'unresolved')
    if location == 'plasmid':
        plasmids = [describe_plasmid(plasmids_by_cluster.get(c, {'cluster_id': c}))
                    for c in plasmid_clusters_for(location_info)]
        return f" - located on plasmid {'; '.join(plasmids)}"
    if location == 'chromosome':
        return " - located on chromosome"
    return " - genomic location unresolved"


def location_components(location_info, plasmids_by_cluster):
    location = location_info.get('location', 'unresolved')
    coding = LOCATION_CODINGS.get(location)
    display = PLASMID_CODE_DISPLAYS[GENOMIC_LOCATION_CODE]
    components = [{
        "code": {"coding": [{"system": CS_KP_PLASMID, "code": GENOMIC_LOCATION_CODE, "display": display}],
                 "text": display},
        "valueCodeableConcept": {"coding": [coding], "text": location} if coding else {"text": "Unresolved"},
    }]
    for cluster_id in plasmid_clusters_for(location_info):
        replicons = ', '.join(plasmids_by_cluster.get(cluster_id, {}).get('replicon_types') or [])
        components.append({"code": {"text": "Plasmid"}, "valueString": f"{cluster_id} {replicons}".strip()})
    return components


def protein_change_component(mutation):
    match = re.search(r'([A-Z\*])(\d+)([A-Z\*]?)', mutation)
    if not match:
        return None
    ref, pos, alt = match.groups()
    if ref not in AMINO_ACIDS or (alt and alt not in AMINO_ACIDS):
        return None
    ref3, alt3 = AMINO_ACIDS[ref], AMINO_ACIDS.get(alt, '?')
    return {"code": {"coding": [{"system": LOINC, "code": "48005-3", "display": "Amino acid change (pHGVS)"}]},
            "valueCodeableConcept": {"coding": [{"system": HGVS, "code": f"p.({ref3}{pos}{alt3})",
                                                 "display": f"{ref3}{pos}{alt3}"}]}}


def build_variant_observation(ctx, idx, source_type, feature, typing_data, gene_locations, plasmids_by_cluster):
    components = []
    gene = feature.get('gene', 'unknown')
    if source_type == 'acquired':
        drug_class = feature.get('drug_class', 'unknown')
        if gene == 'unknown' or not gene:
            return None
        clean_gene = gene.replace('^', '').strip()
        concept = ctx.profile.gene_concept(clean_gene)
    else:
        mutation = feature.get('mutation', 'unknown')
        drug_class = 'unknown'
        if gene in ['GyrA', 'ParC']:
            drug_class = 'fluoroquinolone'
        elif gene in ['OmpK35', 'OmpK36']:
            drug_class = 'carbapenem'
        elif gene in COLISTIN_MUTATION_GENES:
            drug_class = 'colistin'
        clean_gene = f"{gene} {mutation}"
        concept = ctx.profile.mutation_concept(clean_gene)
        phgvs = protein_change_component(mutation)
        if phgvs:
            components.append(phgvs)

    components.append({"code": {"coding": [{"system": LOINC, "code": "48018-6", "display": "Gene studied [ID]"}]},
                       "valueCodeableConcept": concept})
    if drug_class and drug_class != 'unknown':
        snomed = drug_class_snomed(drug_class)
        if snomed:
            components.append({
                "code": {"coding": [{"system": LOINC, "code": "51961-1",
                                     "display": "Genetic variation's effect on drug efficacy"}]},
                "valueCodeableConcept": {"coding": [{"system": SNOMED, "code": snomed['code'],
                                                     "display": snomed['display']}], "text": drug_class}})

    location_info = gene_locations.get(gene if source_type == 'acquired' else clean_gene)
    if location_info:
        components.extend(location_components(location_info, plasmids_by_cluster))
    components.extend(virulence_components(typing_data))

    if source_type == 'acquired':
        text = gene_annotation(gene, drug_class, typing_data, feature.get('intrinsic', False))
    else:
        text = f"Point mutation {mutation} in {gene} conferring resistance"
    if location_info:
        text += describe_location(location_info, plasmids_by_cluster)

    obs_id = sanitize_id(f"{ctx.fid}-resistance-{idx + 1}")
    observation = ctx.observation("variant", obs_id, {
        "text": {"status": "generated", "div": div(xhtml_text(text))},
        "status": "final",
        "category": [CATEGORY_LAB, CATEGORY_GE],
        "code": {"coding": [{"system": LOINC, "code": "69548-6", "display": "Genetic variant assessment"}],
                 "text": f"Resistance Variant: {clean_gene}"},
        "valueCodeableConcept": {"coding": [{"system": LOINC, "code": "LA9633-4", "display": "Present"}],
                                 "text": "Present"},
    })
    observation["component"] = components
    clusters = plasmid_clusters_for(location_info) if location_info else []
    return observation, clusters


def build_no_resistance_observation(ctx, typing_data):
    sequence_type = typing_data.get('mlst', {}).get('sequence_type', 'Unknown')
    k_type = typing_data.get('capsule', {}).get('k_type', 'Unknown')
    text = "Klebsiella pneumoniae isolate with no antimicrobial resistance genes detected"
    if sequence_type not in ['Unknown', 'Undetermined']:
        text += f" - MLST: {sequence_type}"
    if k_type not in ['Unknown', 'Capsule null']:
        text += f" - Capsule: {k_type}"
    observation = ctx.observation("variant", f"{ctx.fid}-susceptible", {
        "text": {"status": "generated", "div": div(xhtml_text(text))},
        "status": "final",
        "category": [CATEGORY_LAB, CATEGORY_GE],
        "code": {"coding": [{"system": LOINC, "code": "69548-6", "display": "Genetic variant assessment"}],
                 "text": "Antimicrobial Susceptibility Assessment"},
        "valueCodeableConcept": {
            "coding": [{"system": LOINC, "code": "LA9634-2", "display": "Absent"},
                       {"system": SNOMED, "code": "444018003", "display": "No antimicrobial resistance detected"}],
            "text": "No resistance genes detected"},
    })
    observation["component"] = virulence_components(typing_data)
    return observation


def typing_body(text, code, value):
    return {"text": {"status": "generated", "div": div(text)}, "status": "final",
            "category": [CATEGORY_LAB, CATEGORY_GE], "code": code, **value}


def build_mlst_observation(ctx, typing_data):
    mlst = typing_data.get('mlst', {})
    sequence_type = mlst.get('sequence_type', 'Unknown')
    alleles = {locus: allele for locus, allele in mlst.get('alleles', {}).items()
               if allele not in ['Unknown', '-', '', 'NA']}
    if sequence_type in ['Unknown', 'Undetermined', 'NA']:
        if not alleles:
            return None
        if len(alleles) < 7 or any('*' in a or '?' in a for a in alleles.values()):
            sequence_type = 'ST not determined (inexact or missing alleles)'
        else:
            sequence_type = 'Novel ST (new allele combination)'
        st_value = {"text": sequence_type}
    else:
        st_value = {"coding": [ctx.profile.sequence_type_coding(sequence_type)], "text": sequence_type}

    allele_text = ', '.join(f"{locus} {allele}" for locus, allele in alleles.items())
    text = f"Klebsiella pneumoniae MLST {xhtml_text(sequence_type)}"
    text += f". Alleles: {xhtml_text(allele_text)}" if allele_text else ""
    body = typing_body(text, {"coding": [{"system": LOINC, "code": "612-2",
                                          "display": "Bacterial strain [Type] in Isolate by Bacteria subtyping"}]},
                       {"valueCodeableConcept": st_value})
    if alleles:
        body["component"] = [{"code": {"text": f"MLST allele {locus}"}, "valueString": allele}
                             for locus, allele in alleles.items()]
    return ctx.observation("typing", sanitize_id(f"{ctx.fid}-mlst"), body)


def build_virulence_observation(ctx, typing_data):
    try:
        score = int(typing_data.get('virulence', {}).get('virulence_score', 'Unknown'))
    except ValueError:
        return None
    body = typing_body(f"Virulence score: {score}", ctx.profile.observation_code("virulence-score"),
                       {"valueQuantity": ctx.profile.score_quantity(score)})
    return ctx.observation("typing", sanitize_id(f"{ctx.fid}-virulence-score"), body)


def build_capsule_observation(ctx, typing_data):
    capsule = typing_data.get('capsule', {})
    k_type = capsule.get('k_type', 'Unknown')
    k_locus = capsule.get('k_locus', 'Unknown')
    if k_type in ['Unknown', 'NA'] and k_locus in ['Unknown', 'NA']:
        return None
    if k_type in ['Capsule null', 'Unknown', 'NA']:
        phenotype = 'capsule null' if k_type == 'Capsule null' else 'type unknown'
        k_type = f"{k_locus} ({phenotype})"
        coding = ctx.profile.capsule_coding(k_locus, f"Klebsiella capsule locus {k_locus} ({phenotype})")
    else:
        coding = ctx.profile.capsule_coding(k_type, f"Klebsiella capsule type {k_type}")
    body = typing_body(f"Klebsiella pneumoniae capsule type {xhtml_text(k_type)}",
                       ctx.profile.observation_code("capsule-type"),
                       {"valueCodeableConcept": {"coding": [coding], "text": k_type}})
    return ctx.observation("typing", sanitize_id(f"{ctx.fid}-capsule"), body)


def build_o_antigen_observation(ctx, typing_data):
    o_type = typing_data.get('o_locus', {}).get('o_type', 'Unknown')
    if o_type in ['Unknown', 'NA', '-', '']:
        return None
    body = typing_body(f"Klebsiella pneumoniae O-antigen type {xhtml_text(o_type)}",
                       {"coding": [{"system": LOINC, "code": "612-2",
                                    "display": "Bacterial strain [Type] in Isolate by Bacteria subtyping"}]},
                       {"valueCodeableConcept": {"coding": [ctx.profile.o_antigen_coding(o_type)],
                                                 "text": o_type}})
    return ctx.observation("typing", sanitize_id(f"{ctx.fid}-o-antigen"), body)


def susceptibility_interpretation(status):
    if status == 'Resistant':
        return {"coding": [{"system": LOINC, "code": "LA6676-6", "display": "Resistant"}], "text": "Resistant"}
    return {"coding": [{"system": LOINC, "code": "LA24225-7", "display": "Susceptible"}], "text": "Susceptible"}


def build_panel_observation(ctx, typing_data):
    drug_status = predict_panel(typing_data, log)
    components = [{"code": {"coding": [{"system": LOINC, "code": DRUG_LOINC.get(drug),
                                        "display": f"{drug} [Susceptibility]"}], "text": drug},
                   "valueCodeableConcept": susceptibility_interpretation(status)}
                  for drug, status in drug_status.items()]
    total = len(drug_status)
    resistant = sum(1 for s in drug_status.values() if s == "Resistant")
    acquired = resistant - 1
    rows = "".join(f'<tr><td>{drug}</td><td><strong>{status}</strong></td></tr>'
                   for drug, status in drug_status.items())
    div_text = f"""<div xmlns="{XHTML}">
    <h3>Klebsiella pneumoniae Predicted Susceptibility Panel</h3>
    <p><strong>Sample:</strong> {ctx.fid}</p>
    <p><strong>Method:</strong> Whole Genome Sequencing - Genotypic Prediction</p>
    <p><strong>Panel:</strong> {total} antimicrobials</p>
    <p><strong>Summary:</strong> {resistant} Resistant (including 1 intrinsic), {total - resistant} Susceptible</p>
    <table border="1">
        <tr><th>Antimicrobial</th><th>Predicted Result</th></tr>
        {rows}
    </table>
    </div>"""
    observation = ctx.observation("panel", sanitize_id(f"{ctx.fid}-susceptibility-panel"), {
        "text": {"status": "generated", "div": div_text},
        "status": "final",
        "category": [CATEGORY_LAB, CATEGORY_GE],
        "code": {"coding": [{"system": LOINC, "code": "29576-6", "display": "Bacterial susceptibility panel"}],
                 "text": "Klebsiella pneumoniae Susceptibility Panel (Genotypic Prediction)"},
    })
    observation["interpretation"] = [{
        "coding": [{"system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
                    "code": "N" if acquired == 0 else "A", "display": "Normal" if acquired == 0 else "Abnormal"}],
        "text": ('No acquired resistance detected' if acquired == 0
                 else f'Acquired resistance detected ({acquired} drugs)')}]
    observation["component"] = components
    return observation


def load_cgmlst(path, sample_id):
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().strip().split('\t')
        matches = []
        for line in handle:
            parts = line.strip().split('\t')
            if len(parts) >= 2 and assembly_sample_id(parts[0]) == str(sample_id):
                matches.append(parts)
    if len(matches) > 1:
        raise InputError(f"{path}: {len(matches)} cgMLST rows for sample {sample_id!r}")
    if not matches:
        return None
    alleles = matches[0][1:]
    total = len(alleles)
    called = sum(1 for a in alleles if a.isdigit() and a != '0')
    return {'total_loci': total, 'called_loci': called,
            'percentage': (called / total) * 100 if total > 0 else 0,
            'alleles': {header[i + 1]: a for i, a in enumerate(alleles) if i + 1 < len(header)}}


def build_cgmlst_observation(ctx, info):
    profile = ctx.profile
    div_text = f"""<div xmlns="{XHTML}">
                <strong>Core Genome MLST Analysis</strong><br/>
                Called Loci: {info['called_loci']}/{info['total_loci']} ({info['percentage']:.1f}%)<br/>
                Detailed allele profile included in components.
                </div>"""
    components = [
        {"code": profile.component_code("cgmlst-loci-called"),
         "valueQuantity": {"value": info['called_loci'], "unit": "{count}", "system": UCUM, "code": "1"}},
        {"code": profile.component_code("cgmlst-loci-total"),
         "valueQuantity": {"value": info['total_loci'], "unit": "{count}", "system": UCUM, "code": "1"}},
        {"code": profile.component_code("cgmlst-call-rate"),
         "valueQuantity": {"value": round(info['percentage'], 1), "unit": "%", "system": UCUM, "code": "%"}},
    ]
    for locus, allele in info['alleles'].items():
        value = {"valueInteger": int(allele)} if allele.isdigit() else {"valueString": allele}
        components.append({"code": {"coding": [{"system": CGMLST_NCS, "code": locus, "display": locus}]}, **value})
    observation = ctx.observation("typing", sanitize_id(f"{ctx.fid}-cgmlst"), {
        "status": "final",
        "category": [CATEGORY_LAB],
        "code": profile.observation_code("cgmlst-profile"),
    }, specimen=False)
    observation["component"] = components
    observation["text"] = {"status": "generated", "div": div_text}
    return observation


def plasmid_observation(ctx, obs_id, code, div_text):
    return ctx.observation("typing", obs_id, {
        "text": {"status": "generated", "div": div(html.escape(div_text, quote=False))},
        "status": "final",
        "category": [CATEGORY_LAB, CATEGORY_GE],
        "code": {"coding": [{"system": CS_KP_PLASMID, "code": code, "display": PLASMID_CODE_DISPLAYS[code]}],
                 "text": PLASMID_CODE_DISPLAYS[code]},
    })


def build_plasmid_observations(ctx, plasmid_data, plasmid_members):
    status = plasmid_data.get('status', 'failed')
    plasmids = plasmid_data.get('plasmids', []) if status != 'failed' else []
    confidence = plasmid_data.get('confidence', 'unknown')
    tools = plasmid_data.get('tools', {})
    observations = []

    for plasmid in plasmids:
        cluster_id = plasmid.get('cluster_id') or 'unclustered'
        replicons = ', '.join(plasmid.get('replicon_types') or []) or 'Replicon not typed'
        genes = plasmid.get('amr_genes') or []
        size_bp = plasmid.get('size_bp')
        mobility = plasmid.get('predicted_mobility')

        cluster_text = cluster_id
        if plasmid.get('secondary_cluster_id'):
            cluster_text += f" (secondary {plasmid['secondary_cluster_id']})"
        components = [text_component("Plasmid cluster", cluster_text)]
        if size_bp:
            components.append({"code": {"text": "Plasmid size"}, "valueQuantity": {"value": size_bp, "unit": "bp"}})
        if plasmid.get('num_contigs'):
            components.append({"code": {"text": "Number of contigs"}, "valueInteger": plasmid['num_contigs']})
        components.append(text_component("Relaxase type", plasmid.get('relaxase_types')))
        components.append(text_component("MPF type", plasmid.get('mpf_type')))
        components.append(text_component("oriT type", plasmid.get('orit_types')))
        components.append(text_component("Predicted mobility", mobility))
        if plasmid.get('mash_nearest_neighbor'):
            neighbor = plasmid['mash_nearest_neighbor']
            if plasmid.get('mash_distance') is not None:
                neighbor += f" (Mash distance {plasmid['mash_distance']})"
            components.append(text_component("Nearest reference plasmid", neighbor))
        components.append(text_component("Predicted host range", plasmid.get('predicted_host_range')))
        components.append(text_component("Carried resistance genes", genes))
        components.append(text_component("Carried carbapenemase genes", plasmid.get('carbapenemase_genes')))
        components.append(text_component("Carried ESBL genes", plasmid.get('esbl_genes')))
        components.append(text_component("Analysis confidence", confidence))

        size_text = f"{size_bp / 1000:.1f} kb" if size_bp else "size unknown"
        div_text = (f"Plasmid {cluster_id}: {replicons}, {size_text}, {mobility or 'mobility unknown'}; "
                    f"carries {', '.join(genes) if genes else 'no resistance genes'}. Confidence: {confidence}.")
        observation = plasmid_observation(ctx, sanitize_id(f"{ctx.fid}-plasmid-{cluster_id}"), PLASMID_CODE,
                                          div_text)
        observation["valueCodeableConcept"] = {"text": replicons}
        observation["component"] = [c for c in components if c]
        members = plasmid_members.get(cluster_id, [])
        if members:
            observation["hasMember"] = [{"reference": f"Observation/{obs_id}"} for obs_id in members]
        observations.append(observation)

    if status == 'failed':
        div_text = "Plasmid analysis failed; plasmid content is unknown."
    else:
        div_text = (f"Plasmid analysis (MOB-suite {tools.get('mob_suite', 'unknown')}): "
                    f"{len(plasmids)} plasmid(s) detected. Confidence: {confidence}.")
    summary = plasmid_observation(ctx, sanitize_id(f"{ctx.fid}-plasmid-summary"), PLASMID_ANALYSIS_CODE, div_text)
    if status == 'failed':
        summary["dataAbsentReason"] = {"coding": [{
            "system": "http://terminology.hl7.org/CodeSystem/data-absent-reason", "code": "error",
            "display": "Error"}]}
    else:
        summary["valueInteger"] = len(plasmids)
    summary["component"] = [c for c in (
        text_component("Analysis status", status),
        text_component("Assembly type", plasmid_data.get('assembly_type')),
        text_component("Analysis confidence", confidence),
        text_component("MOB-suite version", tools.get('mob_suite')),
        text_component("Resistance gene reference", tools.get('card_db'))) if c]
    if observations:
        summary["hasMember"] = [{"reference": f"Observation/{obs['id']}"} for obs in observations]
    return [summary] + observations


def build_qc_observation(ctx, qc_data):
    status = qc_data.get('status', 'UNKNOWN')
    warnings = qc_data.get('warnings') or []
    components = [text_component("Sequencing platform", qc_data.get('platform'))]
    for label, key in [("Estimated read coverage", 'estimated_coverage'), ("Minimum read coverage", 'min_coverage')]:
        if qc_data.get(key) is not None:
            components.append({"code": {"text": label}, "valueQuantity": {"value": qc_data[key], "unit": "x"}})
    components.append(text_component("Coverage status", qc_data.get('coverage_status')))
    if qc_data.get('polishing'):
        polishing = f"{qc_data.get('polishing_tool') or 'polishing'} {qc_data['polishing']}"
        if qc_data.get('polishing_model'):
            polishing += f" ({qc_data['polishing_model']})"
        components.append(text_component("Assembly polishing", polishing))
    components.extend(text_component("Quality warning", warning) for warning in warnings)

    div_text = f"Sequencing QC: {status}"
    div_text += f" - {'; '.join(warnings)}. Results for this sample may be unreliable." if warnings else "."
    observation = ctx.observation("typing", sanitize_id(f"{ctx.fid}-sequencing-qc"), {
        "text": {"status": "generated", "div": div(html.escape(div_text, quote=False))},
        "status": "final",
        "category": [CATEGORY_LAB],
        "code": ctx.profile.observation_code("sequencing-qc"),
    })
    observation["valueString"] = status
    observation["component"] = [c for c in components if c]
    if warnings:
        observation["interpretation"] = [{"coding": [{
            "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
            "code": "A", "display": "Abnormal"}]}]
    return observation


def component_value(component):
    if 'valueString' in component:
        return component['valueString']
    if 'valueCodeableConcept' in component:
        return component['valueCodeableConcept'].get('text', '')
    if 'valueQuantity' in component:
        quantity = component['valueQuantity']
        return f"{quantity.get('value')} {quantity.get('unit', '')}".strip()
    if 'valueInteger' in component:
        return str(component['valueInteger'])
    return ''


def labelled(observation):
    return {c.get('code', {}).get('text', ''): component_value(c) for c in observation.get('component', [])}


def typing_summary(records):
    sequence_type = capsule_type = o_type = None
    alleles = {}
    virulence_score = 0
    for kind, obs in records:
        value = obs.get('valueCodeableConcept', {})
        codings = value.get('coding', [])
        if kind == 'mlst':
            sequence_type = value.get('text', '') or (codings[0].get('code', '') if codings else '')
            if sequence_type and not sequence_type.startswith('ST') and sequence_type.isdigit():
                sequence_type = f"ST{sequence_type}"
            for comp in obs.get('component', []):
                text = comp.get('code', {}).get('text', '')
                if 'valueString' in comp and 'MLST allele' in text:
                    alleles[text.replace('MLST allele ', '')] = comp['valueString']
        elif kind == 'o_antigen':
            o_type = value.get('text') or (codings[0].get('code', '') if codings else '')
        elif kind == 'virulence':
            virulence_score = obs.get('valueQuantity', {}).get('value', 0)
        elif kind == 'capsule':
            capsule_type = value.get('text')
    return sequence_type, alleles, capsule_type, virulence_score, o_type


def plasmid_summary_info(records):
    summary = None
    plasmids = []
    gene_locations = {}
    for kind, obs in records:
        if kind == 'plasmid_summary':
            summary = labelled(obs)
        elif kind == 'plasmid':
            plasmids.append({'replicons': obs.get('valueCodeableConcept', {}).get('text', ''), **labelled(obs)})
        elif kind == 'variant':
            gene = location = None
            plasmid_labels = []
            for comp in obs.get('component', []):
                coding = comp.get('code', {}).get('coding', [{}])[0]
                if 'Gene studied' in coding.get('display', '') and 'valueCodeableConcept' in comp:
                    gene = comp['valueCodeableConcept'].get('text', '')
                elif coding.get('code') == GENOMIC_LOCATION_CODE:
                    location = comp.get('valueCodeableConcept', {}).get('text', '')
                elif comp.get('code', {}).get('text') == 'Plasmid':
                    plasmid_labels.append(comp.get('valueString', ''))
            if gene and location:
                gene_locations[gene] = f"plasmid {', '.join(plasmid_labels)}" if plasmid_labels else location.lower()
    return summary, plasmids, gene_locations


def qc_warnings_of(records):
    warnings = []
    for kind, obs in records:
        if kind == 'qc':
            warnings.extend(c.get('valueString', '') for c in obs.get('component', [])
                            if c.get('code', {}).get('text') == 'Quality warning')
    return [w for w in warnings if w]


def plasmid_borne_findings(plasmids):
    findings = []
    for plasmid in plasmids:
        replicons = plasmid.get('replicons', '')
        if not replicons or replicons == 'Replicon not typed':
            replicons = 'untyped'
        mobility = plasmid.get('Predicted mobility', '')
        note = plasmid.get('Plasmid cluster', '').split(' ')[0]
        if plasmid.get('Analysis confidence', '').startswith('predicted'):
            note += "; predicted from short-read assembly"
        for label, kind in [('Carried carbapenemase genes', 'Carbapenemase'), ('Carried ESBL genes', 'ESBL')]:
            if plasmid.get(label):
                carrier = f"{mobility} {replicons}".strip()
                findings.append((kind, f"{kind} {plasmid[label]} on {carrier} plasmid ({note})"))
    return findings


def report_html(ctx, records, resistance_class, conclusion, typing, qc_warnings, gene_location_labels,
                plasmid_summary, plasmids):
    sequence_type, alleles, capsule_type, virulence_score, _o_type = typing
    content = (f"<div xmlns=\"{XHTML}\">\n"
               f"<h1>{REPORT_TITLE}</h1>\n"
               f"<p><strong>Sample ID:</strong> {ctx.fid}</p>\n"
               f"<p><strong>Report Date:</strong> {now().strftime('%Y-%m-%d %H:%M:%S UTC')}</p>\n"
               f"<p><strong>Resistance Classification:</strong> {resistance_class}</p>")
    if qc_warnings:
        content += (f"<p><strong>Quality warning:</strong> {html.escape('; '.join(qc_warnings))}. "
                    f"Interpret these results with caution.</p>")
    if sequence_type:
        content += f"<p><strong>MLST Sequence Type:</strong> {sequence_type}</p>"
    if capsule_type:
        content += f"<p><strong>Capsule Type:</strong> {capsule_type}</p>"
    if virulence_score > 0:
        content += f"<p><strong>Virulence Score:</strong> {virulence_score}</p>"
    content += f"<p><strong>Conclusion:</strong> {conclusion}</p>"
    for kind, obs in records:
        if kind == 'cgmlst' and 'div' in obs.get('text', {}):
            content += "<hr/>" + obs['text']['div']

    resistance_genes, other_genes = [], []
    for kind, obs in records:
        for comp in obs.get('component', []):
            if 'valueCodeableConcept' not in comp:
                continue
            display = comp.get('code', {}).get('coding', [{}])[0].get('display', '')
            text = comp['valueCodeableConcept'].get('text', '')
            if 'Gene studied' in display and text and 'ST' not in text:
                if any(k in text.lower() for k in ['kpc', 'ndm', 'ctx-m', 'tem', 'shv', 'oxa', 'vim', 'imp', 'mcr']):
                    resistance_genes.append(text)
                else:
                    other_genes.append(text)
    for title, genes in (("Detected Resistance Genes", resistance_genes), ("Other Genetic Features", other_genes)):
        if genes:
            content += f"<h2>{title}</h2><ul>"
            for gene in sorted(set(genes)):
                location = f" ({gene_location_labels[gene]})" if gene in gene_location_labels else ""
                content += f"<li>{gene}{location}</li>"
            content += "</ul>"

    if plasmid_summary is not None:
        content += "<h2>Plasmid Content</h2>"
        if plasmid_summary.get('Analysis status') == 'failed':
            content += "<p>Plasmid analysis failed; plasmid content is unknown.</p>"
        elif not plasmids:
            content += "<p>No plasmids detected.</p>"
        else:
            content += ("<table><tr><th>Cluster</th><th>Replicons</th><th>Size</th><th>Mobility</th>"
                        "<th>Resistance genes</th></tr>")
            for plasmid in plasmids:
                cells = [plasmid.get('Plasmid cluster'), plasmid.get('replicons'), plasmid.get('Plasmid size'),
                         plasmid.get('Predicted mobility'), plasmid.get('Carried resistance genes')]
                content += "<tr>" + "".join(f"<td>{html.escape(cell or '-', quote=False)}</td>"
                                            for cell in cells) + "</tr>"
            content += "</table>"
        if plasmid_summary.get('Analysis confidence'):
            content += (f"<p><em>Plasmid assignment confidence: "
                        f"{html.escape(plasmid_summary['Analysis confidence'], quote=False)}.</em></p>")

    if alleles:
        content += "<h2>MLST Allele Profile</h2><ul>"
        for allele, value in sorted(alleles.items()):
            content += f"<li>{allele}: {value}</li>"
        content += "</ul>"
    return content + "</div>"


def build_diagnostic_report(ctx, records):
    gene_obs = [obs for kind, obs in records if kind == 'variant']
    panel = next((obs for kind, obs in records if kind == 'panel'), None)
    resistance_class, description = classify_resistance_profile(gene_obs, panel, log)
    typing = typing_summary(records)
    sequence_type, _alleles, capsule_type, virulence_score, o_type = typing

    parts = [description]
    if sequence_type:
        parts.append(f"MLST {sequence_type}")
    if capsule_type:
        parts.append(f"Capsule type: {capsule_type}")
    if virulence_score > 0:
        parts.append(f"Virulence score: {virulence_score}")
    plasmid_summary, plasmids, gene_location_labels = plasmid_summary_info(records)
    findings = plasmid_borne_findings(plasmids)
    parts.extend(finding for _kind, finding in findings)
    qc_warnings = qc_warnings_of(records)
    if qc_warnings:
        parts.append(f"Quality warning: {'; '.join(qc_warnings)}; interpret these results with caution")
    conclusion = ". ".join(part.rstrip('.') for part in parts) + "."

    display_text = "Other Resistant" if resistance_class == "Resistant" else resistance_class
    coding = ctx.profile.conclusion_coding(resistance_class)
    codes = [{"coding": [coding], "text": display_text} if coding else {"text": display_text}]
    if sequence_type:
        codes.append({"text": f"MLST {sequence_type}"})
    if o_type:
        codes.append({"text": f"O-antigen {o_type}"})
    if capsule_type:
        codes.append({"text": f"Capsule {capsule_type}"})
    for kind in sorted({kind for kind, _f in findings}):
        codes.append({"text": f"Plasmid-borne {kind.lower() if kind == 'Carbapenemase' else kind}"})
    if qc_warnings:
        codes.append({"text": "Sequencing quality warning"})

    content = report_html(ctx, records, resistance_class, conclusion, typing, qc_warnings, gene_location_labels,
                          plasmid_summary, plasmids)
    report = {
        "resourceType": "DiagnosticReport",
        "id": f"{ctx.fid}-genomic-report",
        "text": {"status": "generated",
                 "div": div(xhtml_text(f"KP genomic analysis report for sample {ctx.fid}. {conclusion}"))},
        "status": "final",
        "category": [CATEGORY_LAB, CATEGORY_GE],
        "code": ctx.profile.report_code(),
        "subject": {"reference": f"Patient/{ctx.fid}-patient"},
        "effectiveDateTime": ctx.timestamp,
        "issued": ctx.timestamp,
        "result": [{"reference": f"Observation/{obs['id']}"} for _kind, obs in records if obs.get('id')],
        "specimen": [{"reference": f"Specimen/{ctx.fid}-specimen"}],
        "conclusion": conclusion,
        "conclusionCode": codes,
        "presentedForm": [{"contentType": "text/html", "language": "en-US", "title": REPORT_TITLE,
                           "data": base64.b64encode(content.encode('utf-8')).decode('ascii')}],
    }
    org_id = ctx.performer["reference"].split("/", 1)[1] if ctx.performer else None
    identifier = ctx.profile.identifier("diagnostic-report", f"KP-GEN-{ctx.fid}-{today_compact()}", org_id)
    if identifier:
        report["identifier"] = [identifier]
    if ctx.performer:
        report["performer"] = [dict(ctx.performer)]
    return with_meta(report, ctx.profile.meta("report")), resistance_class


def build_device(profile, version, source_repository):
    device = {
        "resourceType": "Device",
        "id": f"kp-genomic-analysis-{sanitize_id(version) if version else 'unknown'}",
        "status": "active",
        "deviceName": [{"name": "KPtoFHIR", "type": "manufacturer-name"}],
        "type": {"coding": [{"system": SNOMED, "code": "706689003", "display": "Application program software"}],
                 "text": "Bioinformatics analysis pipeline"},
        "version": [{"type": {"text": "software-version"}, "value": str(version)}],
    }
    if source_repository:
        device["property"] = [{
            "type": {"coding": [{"system": CS_PIPELINE_PROPERTY, "code": "source-repository",
                                 "display": "Source repository"}], "text": "Source repository"},
            "valueCode": [{"text": str(source_repository)}]}]
    return with_meta(device, profile.meta("device"))


def build_provenance(ctx, targets, device_id, source_description):
    participant = "http://terminology.hl7.org/CodeSystem/provenance-participant-type"
    agents = [{"type": {"coding": [{"system": participant, "code": "assembler", "display": "Assembler"}]},
               "who": {"reference": f"Device/{device_id}", "display": "KPtoFHIR"}}]
    if ctx.performer:
        agents.append({"type": {"coding": [{"system": participant, "code": "performer", "display": "Performer"}]},
                       "who": dict(ctx.performer)})
    provenance = {
        "resourceType": "Provenance",
        "id": f"{ctx.fid}-provenance",
        "target": targets,
        "recorded": ctx.timestamp,
        "activity": {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v3-DataOperation",
                                 "code": "CREATE", "display": "create"}]},
        "agent": agents,
        "entity": [{"role": "derivation", "what": {"display": source_description}}],
    }
    return with_meta(provenance, ctx.profile.meta("provenance"))


def load_json(path, label, fid):
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError) as exc:
        raise InputError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise InputError(f"{label} {path} is not a JSON object")
    if data.get('sample_id') and fhir_id(data['sample_id']) != fid:
        raise InputError(f"{label} {path} is for sample {data['sample_id']!r}, not {fid!r}")
    return data


def build_records(ctx, typing_data, plasmid_data, qc_data, cgmlst_info):
    plasmids_by_cluster = {p.get('cluster_id'): p for p in (plasmid_data or {}).get('plasmids', [])}
    gene_locations = (plasmid_data or {}).get('gene_locations', {})
    resistance = typing_data.get('resistance', {})
    features = ([('acquired', g) for g in resistance.get('detailed_resistance', [])] +
                [('mutation', m) for m in resistance.get('point_mutations', [])])

    records = []
    plasmid_members = defaultdict(list)
    for idx, (source_type, feature) in enumerate(features):
        built = build_variant_observation(ctx, idx, source_type, feature, typing_data, gene_locations,
                                          plasmids_by_cluster)
        if built is None:
            continue
        observation, clusters = built
        records.append(('variant', observation))
        for cluster_id in clusters:
            plasmid_members[cluster_id].append(observation['id'])
    if not features:
        records.append(('susceptible', build_no_resistance_observation(ctx, typing_data)))

    for kind, builder in (('mlst', build_mlst_observation), ('virulence', build_virulence_observation),
                          ('capsule', build_capsule_observation), ('o_antigen', build_o_antigen_observation),
                          ('panel', build_panel_observation)):
        observation = builder(ctx, typing_data)
        if observation:
            records.append((kind, observation))
    if cgmlst_info:
        records.append(('cgmlst', build_cgmlst_observation(ctx, cgmlst_info)))
    if plasmid_data is not None:
        observations = build_plasmid_observations(ctx, plasmid_data, plasmid_members)
        records.append(('plasmid_summary', observations[0]))
        records.extend(('plasmid', obs) for obs in observations[1:])
    if qc_data is not None:
        records.append(('qc', build_qc_observation(ctx, qc_data)))
    return records


def parse_args(argv=None):
    p = argparse.ArgumentParser(description='Write the PHI-free genomic FHIR Bundle for one K. pneumoniae sample')
    p.add_argument('--typing_json', required=True, help="this sample's <id>_typing.json from kp_parse_kleborate.py")
    p.add_argument('--sample_id', required=True, help='sample identifier (exact)')
    p.add_argument('--output', required=True, help='output Bundle JSON')
    p.add_argument('--plasmid_json', help="this sample's <id>_plasmid.json from parse_plasmids.py")
    p.add_argument('--qc_json', help="this sample's <id>_assembly_qc.json from assembly_qc.py")
    p.add_argument('--cgmlst_file', help='chewBBACA results_alleles.tsv for the run')
    p.add_argument('--fhir_profile', default='generic', choices=['generic', 'satusehat'])
    p.add_argument('--identifier_system_base', default=None,
                   help='base URI for business identifiers in the generic profile')
    p.add_argument('--performer_org_id', default=None,
                   help='FHIR id of the performing laboratory Organization (not PHI)')
    p.add_argument('--performer_org_name', default=None)
    p.add_argument('--pipeline_version', default='unknown')
    p.add_argument('--source_repository', default=SOURCE_REPOSITORY)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    profile = get_profile(args.fhir_profile, args.identifier_system_base)
    fid = fhir_id(args.sample_id)
    performer = None
    if args.performer_org_id:
        performer = {"reference": f"Organization/{args.performer_org_id}"}
        if args.performer_org_name:
            performer["display"] = args.performer_org_name
    ctx = SampleContext(args.sample_id, fid, profile, performer, now_iso())

    try:
        typing_data = load_json(args.typing_json, "typing JSON", fid)
        plasmid_data = load_json(args.plasmid_json, "plasmid JSON", fid) if args.plasmid_json else None
        qc_data = load_json(args.qc_json, "assembly QC JSON", fid) if args.qc_json else None
        cgmlst_info = load_cgmlst(args.cgmlst_file, args.sample_id) if args.cgmlst_file else None
    except (InputError, OSError) as exc:
        sys.exit(f"ERROR: {exc}")
    if args.cgmlst_file and cgmlst_info is None:
        print(f"WARNING: no cgMLST row for {args.sample_id} in {args.cgmlst_file}")

    records = build_records(ctx, typing_data, plasmid_data, qc_data, cgmlst_info)
    observations = [obs for _kind, obs in records]
    report, resistance_class = build_diagnostic_report(ctx, records)
    device = build_device(profile, args.pipeline_version, args.source_repository)
    targets = ([{"reference": f"DiagnosticReport/{report['id']}"}] +
               [{"reference": f"Observation/{obs['id']}"} for obs in observations])
    provenance = build_provenance(ctx, targets, device['id'],
                                  f"Kleborate typing results for {fid} ({os.path.basename(args.typing_json)})")

    resources = point_references_at_entries([report] + observations + [device, provenance])
    bundle = {
        "resourceType": "Bundle",
        "id": bundle_id("genomic", fid),
        "type": "collection",
        "timestamp": ctx.timestamp,
        "entry": [{"fullUrl": full_url(r["resourceType"], r["id"]), "resource": r} for r in resources],
    }
    with open(args.output, 'w', encoding='utf-8') as out:
        json.dump(bundle, out, indent=2, ensure_ascii=False)

    variants = sum(1 for kind, _obs in records if kind == 'variant')
    print(f"Summary for {args.sample_id} (FHIR id {fid}, profile {profile.name}):")
    print(f"  resistance observations : {variants}")
    print(f"  observations in total   : {len(observations)}")
    print(f"  plasmid analysis        : {(plasmid_data or {}).get('status', 'not run')}")
    print(f"  cgMLST                  : {'yes' if cgmlst_info else 'no'}")
    print(f"  classification          : {resistance_class}")


if __name__ == "__main__":
    main()
