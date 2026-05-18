#!/usr/bin/env python3

import json
import argparse
import uuid
import sys
import os  
from datetime import datetime, timezone
from clinical_metadata_parser import load_clinical_metadata, find_matching_sample, get_clinical_value, load_organization_metadata, load_practitioner_metadata
import base64
import re

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

def create_patient_resource(sample_id, clinical_data=None, org_data=None):
    org_data = org_data or {}
    org_id = org_data.get('org_id')
    if not clinical_data:
        debug_print(f"No clinical data found for sample {sample_id}")
        return {
            "resourceType": "Patient",
            "id": f"{sample_id}-patient",
            "meta": {
                "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Patient"]
            },
            "active": True,
            "name": [
                {
                    "use": "official",
                    "family": "Patient",
                    "given": [f"KP-{sample_id}"]
                }
            ],
            "gender": "unknown",
            "identifier": [
                {
                    "use": "usual",
                    "type": {
                        "coding": [
                            {
                                "system": "http://terminology.hl7.org/CodeSystem/v2-0203",
                                "code": "MR",
                                "display": "Medical record number"
                            }
                        ]
                    },
                    "system": "http://sys-ids.kemkes.go.id/mr/100007732",
                    "value": sample_id
                }
            ]
        }
    
    family_name = get_clinical_value(clinical_data, 'family_name')
    given_name = get_clinical_value(clinical_data, 'given_name')
    gender = get_clinical_value(clinical_data, 'gender', 'unknown').lower()
    birth_date = get_clinical_value(clinical_data, 'birth_date')
    nik = get_clinical_value(clinical_data, 'nik')
    address = get_clinical_value(clinical_data, 'address')
    city = get_clinical_value(clinical_data, 'city')
    state = get_clinical_value(clinical_data, 'state')
    province_code = get_clinical_value(clinical_data, 'province_code', '')
    city_code = get_clinical_value(clinical_data, 'city_code', '')
    district_code = get_clinical_value(clinical_data, 'district_code', '')
    village_code = get_clinical_value(clinical_data, 'village_code', '')
    citizenship_status = get_clinical_value(clinical_data, 'citizenship_status', 'WNI')
    lat = get_clinical_value(clinical_data, 'latitude', None)
    lon = get_clinical_value(clinical_data, 'longitude', None)

    if gender in ['laki-laki', 'pria', 'male', 'm']:
        gender = "male"
    elif gender in ['perempuan', 'wanita', 'female', 'f']:
        gender = "female"
    else:
        gender = "unknown"

    geo_extensions = []
    if lat and lon:
        try:
            geo_extensions.append({
                "url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                "extension": [
                    {"url": "latitude",  "valueDecimal": float(lat)},
                    {"url": "longitude", "valueDecimal": float(lon)}
                ]
            })
        except (ValueError, TypeError):
            pass

    return {
        "resourceType": "Patient",
        "id": f"{sample_id}-patient",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Patient"]
        },
        "active": True,
        "name": [
            {
                "use": "official",
                "family": family_name,
                "given": [given_name]
            }
        ],
        "gender": gender,
        "birthDate": birth_date,
        "identifier": [
            {
                "use": "official",
                "system": "https://fhir.kemkes.go.id/id/nik",
                "value": nik
            },
            {
                "use": "usual",
                "type": {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/v2-0203",
                            "code": "MR",
                            "display": "Medical record number"
                        }
                    ]
                },
                "system": f"http://sys-ids.kemkes.go.id/mr/{org_id}",
                "value": sample_id
            }
        ],
        "extension": [
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                "extension": [
                    {"url": "province", "valueCode": province_code},
                    {"url": "city",     "valueCode": city_code},
                    {"url": "district", "valueCode": district_code},
                    {"url": "village",  "valueCode": village_code}
                ]
            },
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/citizenshipStatus",
                "valueCode": citizenship_status
            }
        ],
        "address": [
            {
                "use": "home", 
                "type": "physical", 
                "text": address, 
                "city": city,  
                "state": state, 
                "country": "ID",  
                "extension": [
                    {
                        "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                        "extension": [
                            {"url": "province", "valueCode": province_code},
                            {"url": "city",     "valueCode": city_code},
                            {"url": "district", "valueCode": district_code}
                        ]
                    },
                    *geo_extensions
                ]
            }
        ]
    }

def create_organization_resource(org_data=None):
    org_data = org_data or {}
    org_id        = org_data.get('org_id', 'unknown-org')
    name          = org_data.get('name', 'Unknown Organization')
    alias         = org_data.get('alias', '')
    type_code     = org_data.get('type_code', '')
    type_display  = org_data.get('type_display', '')
    type_text     = org_data.get('type_text', '')
    phone         = org_data.get('phone', '')
    email         = org_data.get('email', '')
    address_line  = org_data.get('address_line', '')
    city          = org_data.get('city', '')
    state         = org_data.get('state', '')
    country       = org_data.get('country', 'ID')
    province_code = org_data.get('province_code', '')
    city_code     = org_data.get('city_code', '')
    district_code = org_data.get('district_code', '')
    lat           = org_data.get('latitude', None)
    lon           = org_data.get('longitude', None)

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})
    if email:
        telecom.append({"system": "email", "value": email, "use": "work"})

    addr_extensions = []
    if province_code or city_code or district_code:
        code_ext = {"url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode", "extension": []}
        if province_code:
            code_ext["extension"].append({"url": "province", "valueCode": province_code})
        if city_code:
            code_ext["extension"].append({"url": "city",     "valueCode": city_code})
        if district_code:
            code_ext["extension"].append({"url": "district", "valueCode": district_code})
        addr_extensions.append(code_ext)
    if lat and lon:
        try:
            addr_extensions.append({
                "url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                "extension": [
                    {"url": "latitude",  "valueDecimal": float(lat)},
                    {"url": "longitude", "valueDecimal": float(lon)}
                ]
            })
        except ValueError:
            pass

    resource = {
        "resourceType": "Organization",
        "id": org_id,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Organization"]
        },
        "identifier": [
            {
                "use": "official",
                "system": "http://sys-ids.kemkes.go.id/organization",
                "value": org_id
            }
        ],
        "active": True,
        "type": [
            {
                "coding": [
                    {
                        "system": "http://terminology.kemkes.go.id/CodeSystem/organization-type",
                        "code": type_code,
                        "display": type_display
                    }
                ],
                "text": type_text
            }
        ],
        "name": name,
        "telecom": telecom,
        "address": [
            {
                "use": "work",
                "type": "physical",
                "line": [address_line] if address_line else [],
                "city": city,
                "state": state,
                "country": country,
                "extension": addr_extensions
            }
        ]
    }
    if alias:
        resource["alias"] = [alias]
    return resource


def create_practitioner_resource(practitioner_data=None):
    practitioner_data = practitioner_data or {}
    pid        = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    nik        = practitioner_data.get('nik', '')
    name       = practitioner_data.get('name', 'Unknown Practitioner')
    phone      = practitioner_data.get('phone', '')
    gender     = practitioner_data.get('gender', 'unknown')
    birth_date = practitioner_data.get('birth_date', '')
    str_kki    = practitioner_data.get('str_kki_number', '')
    qual_start = practitioner_data.get('qualification_period_start', '')

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})

    qualification = []
    if str_kki:
        qual = {
            "code": {
                "coding": [{
                    "system": "https://terminology.kemkes.go.id/v1-0302",
                    "code": "STR-KKI",
                    "display": "Surat Tanda Registrasi Dokter"
                }],
                "text": "Surat Tanda Registrasi Dokter"
            }
        }
        qual["identifier"] = [{"system": "https://fhir.kemkes.go.id/id/str-kki-number", "value": str_kki}]
        if qual_start:
            qual["period"] = {"start": qual_start}
        qualification.append(qual)

    resource = {
        "resourceType": "Practitioner",
        "id": pid,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Practitioner"]
        },
        "active": True,
        "name": [{"use": "official", "text": name}],
        "telecom": telecom,
        "gender": gender
    }
    if nik:
        resource["identifier"] = [{
            "use": "official",
            "system": "https://fhir.kemkes.go.id/id/nik",
            "value": nik
        }]
    if birth_date:
        resource["birthDate"] = birth_date
    if qualification:
        resource["qualification"] = qualification
    return resource


def create_practitioner_role_resource(practitioner_data=None, org_data=None):
    practitioner_data = practitioner_data or {}
    org_data     = org_data or {}
    pid          = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    pname        = practitioner_data.get('name', 'Unknown Practitioner')
    phone        = practitioner_data.get('phone', '')
    role_id      = practitioner_data.get('role_id')
    role_code    = practitioner_data.get('role_code', '')
    role_display = practitioner_data.get('role_display', '')
    role_text    = practitioner_data.get('role_text', '')
    org_id       = org_data.get('org_id', 'unknown-org')
    org_name     = org_data.get('name', 'Unknown Organization')

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})

    return {
        "resourceType": "PractitionerRole",
        "id": role_id,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/PractitionerRole"]
        },
        "active": True,
        "practitioner": {
            "reference": f"Practitioner/{pid}",
            "display": pname
        },
        "organization": {
            "reference": f"Organization/{org_id}",
            "display": org_name
        },
        "code": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": role_code,
                        "display": role_display
                    }
                ],
                "text": role_text
            }
        ],
        "telecom": telecom
    }


def create_specimen_resource(sample_id, clinical_data=None, practitioner_data=None, org_data=None):
    org_data          = org_data or {}
    practitioner_data = practitioner_data or {}

    org_id             = org_data.get('org_id', 'unknown-org')
    practitioner_id    = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    practitioner_name  = practitioner_data.get('name', 'Unknown Practitioner')

    if clinical_data:
        given_name        = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name       = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display   = f"{given_name} {family_name}"
        spec_type_code    = get_clinical_value(clinical_data, 'specimen_type_code', '119303007')
        spec_type_display = get_clinical_value(clinical_data, 'specimen_type_display', 'Microbial isolate specimen')
        method_code       = get_clinical_value(clinical_data, 'specimen_collection_method_code', 'SWA')
        method_display    = get_clinical_value(clinical_data, 'specimen_collection_method_display', 'Swab')
        method_text       = get_clinical_value(clinical_data, 'specimen_collection_method_text', 'Clinical specimen collection')
        qty_value         = get_clinical_value(clinical_data, 'specimen_quantity_value', '1')
        qty_unit          = get_clinical_value(clinical_data, 'specimen_quantity_unit', 'mL')
        collected_date    = get_clinical_value(clinical_data, 'specimen_collected_date', None)
        received_date     = get_clinical_value(clinical_data, 'specimen_received_date', None)
    else:
        patient_display   = f"KP Patient {sample_id}"
        spec_type_code    = '119303007'
        spec_type_display = 'Microbial isolate specimen'
        method_code       = 'SWA'
        method_display    = 'Swab'
        method_text       = 'Clinical specimen collection'
        qty_value         = '1'
        qty_unit          = 'mL'
        collected_date    = None
        received_date     = None

    now          = datetime.now(timezone.utc).isoformat()
    collected_dt = collected_date if collected_date and collected_date != 'Unknown' else now
    received_dt  = received_date  if received_date  and received_date  != 'Unknown' else now

    try:
        qty_float = float(qty_value)
    except (ValueError, TypeError):
        qty_float = 1.0

    return {
        "resourceType": "Specimen",
        "id": f"{sample_id}-specimen",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Specimen"]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/specimen/{org_id}",
                "value": f"SPEC-KP-{sample_id}"
            }
        ],
        "status": "available",
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "receivedTime": received_dt,
        "collection": {
            "collectedDateTime": collected_dt,
            "collector": {
                "reference": f"Practitioner/{practitioner_id}",
                "display": practitioner_name
            },
            "method": {
                "coding": [
                    {
                        "system": "http://terminology.hl7.org/CodeSystem/v2-0488",
                        "code": method_code,
                        "display": method_display
                    }
                ],
                "text": method_text
            },
            "quantity": {
                "value": qty_float,
                "unit": qty_unit,
                "system": "http://unitsofmeasure.org",
                "code": qty_unit
            }
        },
        "type": {
            "coding": [
                {
                    "system": "http://snomed.info/sct",
                    "code": spec_type_code,
                    "display": spec_type_display
                }
            ],
            "text": f"{spec_type_display} for KP testing"
        },
        "note": [
            {
                "text": f"Collected clinical specimen from {patient_display} ({sample_id}) for Klebsiella pneumoniae genomic testing"
            }
        ]
    }


def get_resistance_conclusion_coding(resistance_class):

    coding_map = {
        "XDR": {
            "system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term",
            "code": "SP000681",
            "display": "Extensively drug resistant Klebsiella pneumoniae"
        },
        "CRE": {
            "system": "http://snomed.info/sct",
            "code": "1098201000112108",
            "display": "Carbapenemase-producing Klebsiella pneumoniae (organism)"
        },
        "Multidrug-Resistant ESBL": {
            "system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term",
            "code": "SP000687",
            "display": "Multidrug-resistant ESBL-producing Klebsiella pneumoniae"
        },
        "ESBL": {
            "system": "http://snomed.info/sct",
            "code": "409801009",
            "display": "Extended spectrum beta-lactamase producing Klebsiella pneumoniae (organism)"
        },
        "Multidrug-Resistant AmpC": {
            "system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term",
            "code": "SP000683",
            "display": "Multidrug-resistant AmpC-producing Klebsiella pneumoniae"
        },
        "AmpC": {
            "system": "http://snomed.info/sct",
            "code": "1098101000112102",
            "display": "AmpC beta-lactamase producing Klebsiella pneumoniae (organism)"
        },
        "MDR": {
            "system": "http://snomed.info/sct",
            "code": "714315002",
            "display": "Multidrug-resistant Klebsiella pneumoniae (organism)"
        },
        "Susceptible": {
            "system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term",
            "code": "KP-SO",
            "display": "Klebsiella pneumoniae Sensitif Obat"
        },
        "Resistant": { 
            "system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term",
            "code": "KP-RL",
            "display": "Klebsiella pneumoniae Resisten Lain"
        }
    }
    
    classification_to_key = {
        "Extensively Drug-Resistant (XDR)": "XDR",
        "Carbapenem-Resistant Enterobacteriaceae (CRE)": "CRE",
        "Multidrug-Resistant ESBL": "Multidrug-Resistant ESBL",
        "Extended-Spectrum Beta-Lactamase (ESBL)": "ESBL",
        "Multidrug-Resistant AmpC": "Multidrug-Resistant AmpC",
        "AmpC Beta-Lactamase (AmpC)": "AmpC",
        "Multidrug-Resistant (MDR)": "MDR",
        "Resistant": "Resistant",
        "Susceptible": "Susceptible"
    }
    
    if resistance_class in classification_to_key:
        key = classification_to_key[resistance_class]
        return coding_map.get(key)
    
    if resistance_class in coding_map:
        return coding_map[resistance_class]
    
    priority_order = [
        "Multidrug-Resistant ESBL",
        "Multidrug-Resistant AmpC",
        "XDR",
        "CRE", 
        "ESBL",
        "AmpC",
        "MDR"
    ]
    
    for key in priority_order:
        if key in resistance_class:
            return coding_map[key]
    
    if "Susceptible" in resistance_class:
        return coding_map["Susceptible"]
    if "Resistant" in resistance_class:
        return coding_map["Resistant"]
    
    return None

def classify_resistance_profile(observations):

    carbapenemases = []
    esbl_genes = []
    ampc_genes = [] 
    colistin_resistance = []
    other_resistance = []
    
    affected_drug_classes = set()
    specific_drugs_resistant = set()
    
    drug_class_normalization = {
        'beta-lactam': 'beta-lactam',
        'beta_lactam': 'beta-lactam',
        'betalactam': 'beta-lactam',
        'aminoglycoside': 'aminoglycoside',
        'fluoroquinolone': 'fluoroquinolone',
        'quinolone': 'fluoroquinolone',
        'macrolide': 'macrolide',
        'sulfonamide': 'sulfonamide',
        'trimethoprim': 'trimethoprim',
        'tetracycline': 'tetracycline',
        'chloramphenicol': 'chloramphenicol',
        'phenicol': 'chloramphenicol',
        'colistin': 'colistin',
        'polymyxin': 'colistin',
        'carbapenem': 'carbapenem',
        'fosfomycin': 'fosfomycin',
        'rifampicin': 'rifampicin',
        'bla': 'beta-lactam',
        'agly': 'aminoglycoside',
        'flq': 'fluoroquinolone',
        'mls': 'macrolide',
        'mlsb': 'macrolide',
        'sul': 'sulfonamide',
        'tmt': 'trimethoprim',
        'tet': 'tetracycline',
        'col': 'colistin',
        'fos': 'fosfomycin',
        'rif': 'rifampicin',
        'dfr': 'trimethoprim',
    }
    
    def normalize_drug_class(drug_class):
        if not drug_class:
            return None
        dc_lower = drug_class.lower().strip()
        
        if dc_lower in drug_class_normalization:
            return drug_class_normalization[dc_lower]
        
        dc_normalized = dc_lower.replace('-', '_').replace(' ', '_')
        if dc_normalized in drug_class_normalization:
            return drug_class_normalization[dc_normalized]
        
        for key, normalized in drug_class_normalization.items():
            if len(key) >= 3 and (key in dc_lower or dc_lower in key):
                return normalized
        
        return None
    
    narrow_spectrum_shv = [
        'shv-1', 'shv-11', 'shv-26', 'shv-28', 'shv-33', 'shv-36', 'shv-38',
        'shv-41', 'shv-42', 'shv-60', 'shv-62', 'shv-65', 'shv-75', 'shv-76',
        'shv-77', 'shv-79', 'shv-99', 'shv-100', 'shv-108', 'shv-119', 'shv-155',
        'shv-164', 'shv-165', 'shv-166', 'shv-167', 'shv-168', 'shv-169', 'shv-170',
        'shv-171', 'shv-172', 'shv-173', 'shv-174', 'shv-175', 'shv-176', 'shv-177',
        'shv-178', 'shv-179', 'shv-180', 'shv-181', 'shv-182', 'shv-183', 'shv-184',
        'shv-185', 'shv-186', 'shv-187', 'shv-188', 'shv-189', 'shv-190', 'shv-191',
        'shv-192', 'shv-193', 'shv-194', 'shv-195', 'shv-196', 'shv-197', 'shv-198',
        'shv-199', 'shv-200', 'shv-206', 'shv-207'
    ]
    
    esbl_shv_variants = [
        'shv-2', 'shv-3', 'shv-4', 'shv-5', 'shv-6', 'shv-7', 'shv-8', 'shv-9',
        'shv-12', 'shv-13', 'shv-14', 'shv-15', 'shv-16', 'shv-17', 'shv-18',
        'shv-24', 'shv-25', 'shv-27', 'shv-29', 'shv-30', 'shv-31', 'shv-32',
        'shv-34', 'shv-35', 'shv-37', 'shv-39', 'shv-40', 'shv-43', 'shv-44',
        'shv-45', 'shv-46', 'shv-48', 'shv-49', 'shv-52', 'shv-55', 'shv-56',
        'shv-57', 'shv-59', 'shv-61', 'shv-63', 'shv-64', 'shv-66', 'shv-67',
        'shv-68', 'shv-69', 'shv-70', 'shv-71', 'shv-72', 'shv-73', 'shv-74',
        'shv-78', 'shv-80', 'shv-81', 'shv-82', 'shv-83', 'shv-84', 'shv-85',
        'shv-86', 'shv-87', 'shv-88', 'shv-89', 'shv-90', 'shv-91', 'shv-92',
        'shv-93', 'shv-94', 'shv-95', 'shv-96', 'shv-97', 'shv-98', 'shv-101',
        'shv-102', 'shv-103', 'shv-104', 'shv-105', 'shv-106', 'shv-107',
        'shv-109', 'shv-110', 'shv-111', 'shv-112', 'shv-113', 'shv-114',
        'shv-115', 'shv-116', 'shv-117', 'shv-118', 'shv-120', 'shv-121',
        'shv-122', 'shv-123', 'shv-124', 'shv-125', 'shv-126', 'shv-127',
        'shv-128', 'shv-129', 'shv-130', 'shv-131', 'shv-132', 'shv-133',
        'shv-134', 'shv-135', 'shv-136', 'shv-137', 'shv-138', 'shv-139',
        'shv-140', 'shv-141', 'shv-142', 'shv-143', 'shv-144', 'shv-145',
        'shv-146', 'shv-147', 'shv-148', 'shv-149', 'shv-150', 'shv-151',
        'shv-152', 'shv-153', 'shv-154', 'shv-156', 'shv-157', 'shv-158',
        'shv-159', 'shv-160', 'shv-161', 'shv-162', 'shv-163'
    ]
    
    esbl_tem_variants = [
        'tem-3', 'tem-4', 'tem-5', 'tem-6', 'tem-7', 'tem-8', 'tem-9', 'tem-10',
        'tem-12', 'tem-15', 'tem-16', 'tem-17', 'tem-19', 'tem-21', 'tem-22',
        'tem-24', 'tem-25', 'tem-26', 'tem-27', 'tem-28', 'tem-29', 'tem-42',
        'tem-43', 'tem-46', 'tem-47', 'tem-49', 'tem-50', 'tem-52', 'tem-53',
        'tem-60', 'tem-63', 'tem-64', 'tem-70', 'tem-71', 'tem-87', 'tem-88',
        'tem-89', 'tem-90', 'tem-92', 'tem-109', 'tem-116', 'tem-121', 'tem-126',
        'tem-127', 'tem-128', 'tem-129', 'tem-131', 'tem-132', 'tem-134',
        'tem-135', 'tem-136', 'tem-137', 'tem-138', 'tem-139', 'tem-140',
        'tem-141', 'tem-142', 'tem-143', 'tem-144', 'tem-146', 'tem-147',
        'tem-148', 'tem-149', 'tem-150', 'tem-151', 'tem-152', 'tem-153',
        'tem-154', 'tem-155', 'tem-156', 'tem-157', 'tem-158', 'tem-159',
        'tem-160', 'tem-161', 'tem-162', 'tem-163', 'tem-164', 'tem-167',
        'tem-168', 'tem-169', 'tem-170', 'tem-171', 'tem-172', 'tem-176',
        'tem-178', 'tem-179', 'tem-180', 'tem-181', 'tem-182', 'tem-183',
        'tem-184', 'tem-185', 'tem-186', 'tem-187', 'tem-188', 'tem-189',
        'tem-190', 'tem-191', 'tem-192', 'tem-193', 'tem-194', 'tem-195',
        'tem-196', 'tem-197', 'tem-198', 'tem-199', 'tem-200'
    ]
    
    ampc_gene_patterns = [
        'dha', 'cmy', 'mox', 'fox', 'acc', 'act', 'mir', 'lat', 'bil', 'cfea'
    ]
    
    def normalize_gene_name(gene_name):
        if not gene_name:
            return ''
        normalized = gene_name.lower()
        normalized = normalized.replace('^', '').replace('*', '').replace('.', '-').strip()
        normalized = re.sub(r'-v\d+$', '', normalized)
        normalized = re.sub(r'\.v\d+$', '', normalized)
        return normalized
    
    def is_narrow_spectrum_shv(gene_lower):
        normalized = normalize_gene_name(gene_lower)
        for ns_shv in narrow_spectrum_shv:
            if normalized == ns_shv or normalized == ns_shv.replace('-', ''):
                return True
        return False
    
    def is_esbl_gene(gene_lower):
        normalized = normalize_gene_name(gene_lower)
        
        if is_narrow_spectrum_shv(normalized):
            return False
        
        if 'ctx-m' in normalized or 'ctxm' in normalized:
            return True
        
        for esbl_shv in esbl_shv_variants:
            if normalized == esbl_shv or normalized == esbl_shv.replace('-', ''):
                return True
        
        for esbl_tem in esbl_tem_variants:
            if normalized == esbl_tem or normalized == esbl_tem.replace('-', ''):
                return True
        
        return False
    
    for obs in observations:
        div_text = obs.get('text', {}).get('div', '')
        is_intrinsic = '(chromosomal)' in div_text or 'intrinsic' in div_text
        
        code_coding = obs.get('code', {}).get('coding', [{}])[0]
        code_val = code_coding.get('code', '')
        
        if code_val in ['612-2', 'SP000678', 'SP000680', 'SP000682']:
            continue
        
        if code_val == '29576-6':
            components = obs.get('component', [])
            for comp in components:
                drug_text = comp.get('code', {}).get('text', '')
                drug_name = drug_text.split(' - ')[0] if ' - ' in drug_text else drug_text
                
                val_code = comp.get('valueCodeableConcept', {}).get('coding', [{}])[0].get('code')
                if val_code == 'LA6676-6':  
                    specific_drugs_resistant.add(drug_name)
            continue

        components = obs.get('component', [])
        gene = None
        drug_class = None
        
        for component in components:
            if 'valueCodeableConcept' in component:
                code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '')
                value_text = component['valueCodeableConcept'].get('text', '')
                
                if 'Gene studied' in code_display:
                    gene = value_text
                elif 'drug efficacy' in code_display:
                    drug_class = value_text
        
        if gene:
            gene_display = gene.upper().replace('^', '').replace('*', '').strip()
            gene_normalized = normalize_gene_name(gene)
            
            is_carb = any(carb in gene_normalized for carb in ['kpc', 'ndm', 'oxa-48', 'oxa-181', 'oxa-232', 'vim', 'imp', 'ges'])
            
            is_esbl = is_esbl_gene(gene_normalized)
            
            is_ampc = any(ampc in gene_normalized for ampc in ampc_gene_patterns)
            
            is_col = 'mcr' in gene_normalized or (drug_class and 'colistin' in str(drug_class).lower())
            
            is_ns_shv = is_narrow_spectrum_shv(gene_normalized)
            
            if is_ns_shv:
                debug_print(f"Skipping narrow-spectrum SHV: {gene} (intrinsic)")
                continue
            
            if is_intrinsic and not (is_carb or is_esbl or is_ampc or is_col):
                debug_print(f"Skipping chromosomal gene: {gene}")
                continue

            normalized_class = normalize_drug_class(drug_class)
            if normalized_class:
                affected_drug_classes.add(normalized_class)
                debug_print(f"Gene {gene} -> drug class {drug_class} -> normalized: {normalized_class}")

            if is_carb:
                carbapenemases.append(gene_display)
                affected_drug_classes.add('carbapenem')
                affected_drug_classes.add('beta-lactam')
            elif is_esbl:
                esbl_genes.append(gene_display)
                affected_drug_classes.add('beta-lactam')
            elif is_ampc:
                ampc_genes.append(gene_display)
                affected_drug_classes.add('beta-lactam')
            elif is_col:
                colistin_resistance.append(gene_display)
                affected_drug_classes.add('colistin')
            else:
                other_resistance.append(gene_display)
    
    affected_drugs_text = ""
    if specific_drugs_resistant:
        acquired_resistant = [d for d in sorted(specific_drugs_resistant) if d != "Ampicillin"]
        if acquired_resistant:
            affected_drugs_text = f" Predicted resistance to: {', '.join(acquired_resistant)}."

    unique_class_count = len(affected_drug_classes)
    debug_print(f"Total affected drug classes: {unique_class_count} - {affected_drug_classes}")
    
    all_acquired_genes = list(set(carbapenemases + esbl_genes + ampc_genes + colistin_resistance + other_resistance))
    resistance_genes_text = f"Resistance genes: {', '.join(sorted(all_acquired_genes))}." if all_acquired_genes else ""

    if carbapenemases:
        if colistin_resistance:
            return "Extensively Drug-Resistant (XDR)", (
                f"Carbapenem and colistin resistant K. pneumoniae. "
                f"Carbapenemases: {', '.join(set(carbapenemases))}. "
                f"Colistin resistance: {', '.join(set(colistin_resistance))}."
                f"{affected_drugs_text}"
            )
        else:
            return "Carbapenem-Resistant Enterobacteriaceae (CRE)", (
                f"Carbapenem-resistant K. pneumoniae (CP-CRE). "
                f"Carbapenemases detected: {', '.join(set(carbapenemases))}."
                f"{affected_drugs_text}"
            )
    
    elif esbl_genes:
        if unique_class_count >= 3:
            return "Multidrug-Resistant ESBL", (
                f"ESBL-producing K. pneumoniae with additional resistance (MDR). "
                f"ESBL genes: {', '.join(set(esbl_genes))}. "
                f"Resistance affects {unique_class_count} drug classes. "
                f"{resistance_genes_text}"
                f"{affected_drugs_text}"
            )
        else:
            return "Extended-Spectrum Beta-Lactamase (ESBL)", (
                f"ESBL-producing K. pneumoniae. "
                f"ESBL genes: {', '.join(set(esbl_genes))}."
                f"{affected_drugs_text}"
            )
    
    elif ampc_genes:
        if unique_class_count >= 3:
            return "Multidrug-Resistant AmpC", (
                f"AmpC beta-lactamase producing K. pneumoniae with additional resistance (MDR). "
                f"AmpC genes: {', '.join(set(ampc_genes))}. "
                f"Resistance affects {unique_class_count} drug classes. "
                f"{resistance_genes_text}"
                f"{affected_drugs_text}"
            )
        else:
            return "AmpC Beta-Lactamase (AmpC)", (
                f"AmpC beta-lactamase producing K. pneumoniae. "
                f"AmpC genes: {', '.join(set(ampc_genes))}."
                f"{affected_drugs_text}"
            )
    
    elif unique_class_count >= 3:
        return "Multidrug-Resistant (MDR)", (
            f"Multiple resistance genes detected affecting {unique_class_count} drug classes. "
            f"{resistance_genes_text}"
            f"{affected_drugs_text}"
        )
    
    elif other_resistance or unique_class_count > 0:
        return "Resistant", (
            f"Resistance genes detected: {', '.join(set(other_resistance)) if other_resistance else 'See panel'}. "
            f"{affected_drugs_text}"
        )
    
    else:
        return "Susceptible", "No acquired resistance genes detected. Intrinsic ampicillin resistance only."

def extract_mlst_info(observations):
    sequence_type = None
    clonal_complex = None
    alleles = {}
    capsule_type = None
    virulence_score = 0
    o_type = None
    
    for obs in observations:
        code_coding = obs.get('code', {}).get('coding', [{}])[0]
        code_val = code_coding.get('code', '')
        
        if code_val == '612-2':
            val_cc = obs.get('valueCodeableConcept', {})
            val_codings = val_cc.get('coding', [])
            val_system = val_codings[0].get('system', '') if val_codings else ''
            val_text = val_cc.get('text', '')
            val_code = val_codings[0].get('code', '') if val_codings else ''
            
            if val_system == 'http://kaptive.holtlab.net/o-antigen':
                o_type = val_text or val_code
            else:
                sequence_type = val_text
                if not sequence_type and val_code:
                    sequence_type = val_code
                if sequence_type and not sequence_type.startswith('ST') and sequence_type.isdigit():
                    sequence_type = f"ST{sequence_type}"

        elif code_val == 'SP000680':
            val_qty = obs.get('valueQuantity', {})
            virulence_score = val_qty.get('value', 0)
            
        elif code_val == 'SP000678':
            val_cc = obs.get('valueCodeableConcept', {})
            capsule_type = val_cc.get('text')

        components = obs.get('component', [])
        for component in components:
            if 'valueCodeableConcept' in component:
                code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '')
                value_text = component['valueCodeableConcept'].get('text', '')
                
                if 'bacterial subtyping' in code_display.lower() and 'ST' in value_text:
                    sequence_type = value_text
            elif 'valueString' in component:
                code_text = component.get('code', {}).get('text', '')
                value = component.get('valueString')
                
                if 'MLST allele' in code_text:
                    allele_name = code_text.replace('MLST allele ', '')
                    alleles[allele_name] = value
                elif 'Clonal Complex' in code_text:
                    clonal_complex = value
                elif 'K Locus' in code_text or 'Capsule' in code_text:
                    capsule_type = value
    
    return sequence_type, clonal_complex, alleles, capsule_type, virulence_score, o_type

def create_cgmlst_report(sample_id, observations, clinical_data=None, org_data=None, practitioner_data=None):
    org_data          = org_data or {}
    practitioner_data = practitioner_data or {}
    org_id            = org_data.get('org_id', 'unknown-org')
    org_name          = org_data.get('name', 'Unknown Organization')
    practitioner_id   = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    practitioner_name = practitioner_data.get('name', 'Unknown Practitioner')
    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"KP Patient {sample_id}"
    
    report_id = f"{sample_id}-cgmlst-report"
    current_time = datetime.now(timezone.utc).isoformat()
    
    html_content = f"""<div xmlns="http://www.w3.org/1999/xhtml">
<h1>Klebsiella pneumoniae cgMLST Analysis Report</h1>
<p><strong>Patient:</strong> {patient_display}</p>
<p><strong>Sample ID:</strong> {sample_id}</p>
<p><strong>Report Date:</strong> {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}</p>
"""
    
    for obs in observations:
        if 'text' in obs and 'div' in obs['text']:
             html_content += obs['text']['div']
             
    html_content += "</div>"
    
    html_base64 = base64.b64encode(html_content.encode('utf-8')).decode('utf-8')
    
    return {
        "resourceType": "DiagnosticReport",
        "id": report_id,
        "meta": {
            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report"],
            "tag": [
                {
                    "system": "http://terminology.kemkes.go.id/sp",
                    "code": "genomics",
                    "display": "Genomics"
                }
            ]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/diagnostic-report/{org_id}",
                "value": f"KP-CGMLST-{sample_id}-{datetime.now().strftime('%Y%m%d')}"
            }
        ],
        "basedOn": [
            {
                "reference": f"ServiceRequest/{sample_id}-service-request"
            }
        ],
        "status": "final",
        "category": [
            {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                    "code": "laboratory",
                    "display": "Laboratory"
                }]
            },
            {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                    "code": "GE",
                    "display": "Genetics"
                }]
            }
        ],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "SP000682",
                "display": "Core genome MLST [Type] in Isolate by Sequencing"
            }],
            "text": "Klebsiella pneumoniae cgMLST"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "KP Testing Encounter"
        },
        "effectiveDateTime": current_time,
        "issued": current_time,
        "performer": [
            {
                "reference": f"Organization/{org_id}",
                "display": org_name
            },
            {
                "reference": f"Practitioner/{practitioner_id}",
                "display": practitioner_name
            }
        ],
        "result": [{"reference": f"Observation/{obs['id']}"} for obs in observations if obs.get('id')],
        "specimen": [{
            "reference": f"Specimen/{sample_id}-specimen",
            "display": f"Clinical specimen from {patient_display}"
        }],
        "presentedForm": [
            {
                "contentType": "text/html",
                "language": "en-US", 
                "title": "Klebsiella pneumoniae cgMLST Analysis Report",
                "data": html_base64
            }
        ]
    }

def create_diagnostic_report(sample_id, observations, clinical_data=None, org_data=None, practitioner_data=None):
    org_data          = org_data or {}
    practitioner_data = practitioner_data or {}
    org_id            = org_data.get('org_id', 'unknown-org')
    org_name          = org_data.get('name', 'Unknown Organization')
    practitioner_id   = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    practitioner_name = practitioner_data.get('name', 'Unknown Practitioner')
    resistance_class, resistance_description = classify_resistance_profile(observations)
    sequence_type, clonal_complex, alleles, capsule_type, virulence_score, o_type = extract_mlst_info(observations)
    
    conclusion_parts = [resistance_description]
    if sequence_type:
        conclusion_parts.append(f"MLST {sequence_type}")
    if clonal_complex:
        conclusion_parts.append(f"belonging to {clonal_complex}")
    if capsule_type:
        conclusion_parts.append(f"Capsule type: {capsule_type}")
    if virulence_score > 0:
        conclusion_parts.append(f"Virulence score: {virulence_score}")
    
    conclusion = ". ".join(conclusion_parts) + "."
    
    conclusion_codes = []
    resistance_coding = get_resistance_conclusion_coding(resistance_class)
    
    display_text = resistance_class
    if resistance_class == "Resistant":
        display_text = "Other Resistant"
    
    if resistance_coding:
        conclusion_codes.append({
            "coding": [resistance_coding],
            "text": display_text
        })
    else:
        conclusion_codes.append({
            "text": display_text
        })
    
    if sequence_type:
        conclusion_codes.append({
            "text": f"MLST {sequence_type}"
        })
    
    if o_type:
        conclusion_codes.append({
            "text": f"O-antigen {o_type}"
        })
    
    if capsule_type:
        conclusion_codes.append({
            "text": f"Capsule {capsule_type}"
        })
    
    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"KP Patient {sample_id}"
    
    report_id = f"{sample_id}-genomic-report"
    current_time = datetime.now(timezone.utc).isoformat()
    
    html_content = f"""<div xmlns="http://www.w3.org/1999/xhtml">
<h1>Klebsiella pneumoniae Genomic Analysis Report</h1>
<p><strong>Patient:</strong> {patient_display}</p>
<p><strong>Sample ID:</strong> {sample_id}</p>
<p><strong>Report Date:</strong> {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}</p>
<p><strong>Resistance Classification:</strong> {resistance_class}</p>"""

    if sequence_type:
        html_content += f"<p><strong>MLST Sequence Type:</strong> {sequence_type}</p>"
    if clonal_complex:
        html_content += f"<p><strong>Clonal Complex:</strong> {clonal_complex}</p>"
    if capsule_type:
        html_content += f"<p><strong>Capsule Type:</strong> {capsule_type}</p>"
    if virulence_score > 0:
        html_content += f"<p><strong>Virulence Score:</strong> {virulence_score}</p>"
    
    html_content += f"<p><strong>Conclusion:</strong> {conclusion}</p>"
    
    for obs in observations:
        if 'div' in obs.get('text', {}):
             code_coding = obs.get('code', {}).get('coding', [{}])[0]
             if code_coding.get('code') == 'SP000682':
                 html_content += "<hr/>"
                 html_content += obs['text']['div']
    
    if observations:
        resistance_genes = []
        other_genes = []
        
        for obs in observations:
            components = obs.get('component', [])
            for component in components:
                if 'valueCodeableConcept' in component:
                    code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '')
                    value_text = component['valueCodeableConcept'].get('text', '')
                    
                    if 'Gene studied' in code_display and value_text and 'ST' not in value_text:
                        if any(keyword in value_text.lower() for keyword in ['kpc', 'ndm', 'ctx-m', 'tem', 'shv', 'oxa', 'vim', 'imp', 'mcr']):
                            resistance_genes.append(value_text)
                        else:
                            other_genes.append(value_text)
        
        if resistance_genes:
            html_content += "<h2>Detected Resistance Genes</h2><ul>"
            for gene in sorted(set(resistance_genes)):
                html_content += f"<li>{gene}</li>"
            html_content += "</ul>"
        
        if other_genes:
            html_content += "<h2>Other Genetic Features</h2><ul>"
            for gene in sorted(set(other_genes)):
                html_content += f"<li>{gene}</li>"
            html_content += "</ul>"
    
    if alleles:
        html_content += "<h2>MLST Allele Profile</h2><ul>"
        for allele, value in sorted(alleles.items()):
            html_content += f"<li>{allele}: {value}</li>"
        html_content += "</ul>"
    
    html_content += "</div>"
    
    html_base64 = base64.b64encode(html_content.encode('utf-8')).decode('utf-8')
    
    return {
        "resourceType": "DiagnosticReport",
        "id": report_id,
        "meta": {
            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report"],
            "tag": [
                {
                    "system": "http://terminology.kemkes.go.id/sp",
                    "code": "genomics",
                    "display": "Genomics"
                }
            ]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/diagnostic-report/{org_id}",
                "value": f"KP-GEN-{sample_id}-{datetime.now().strftime('%Y%m%d')}"
            }
        ],
        "basedOn": [
            {
                "reference": f"ServiceRequest/{sample_id}-service-request"
            }
        ],
        "status": "final",
        "category": [
            {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                    "code": "laboratory",
                    "display": "Laboratory"
                }]
            },
            {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                    "code": "GE",
                    "display": "Genetics"
                }]
            }
        ],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "81247-9",
                "display": "Master HL7 genetic variant reporting panel"
            }],
            "text": "Klebsiella pneumoniae Genomic Analysis Report"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "KP Testing Encounter"
        },
        "effectiveDateTime": current_time,
        "issued": current_time,
        "performer": [
            {
                "reference": f"Organization/{org_id}",
                "display": org_name
            },
            {
                "reference": f"Practitioner/{practitioner_id}",
                "display": practitioner_name
            }
        ],
        "result": [{"reference": f"Observation/{obs['id']}"} for obs in observations if obs.get('id')],
        "specimen": [{
            "reference": f"Specimen/{sample_id}-specimen",
            "display": f"Clinical specimen from {patient_display}"
        }],
        "conclusion": conclusion,
        "conclusionCode": conclusion_codes,
        "presentedForm": [
            {
                "contentType": "text/html",
                "language": "en-US", 
                "title": "Klebsiella pneumoniae Genomic Analysis Report",
                "data": html_base64
            }
        ]
    }

def create_service_request_resource(sample_id, clinical_data=None, practitioner_data=None, org_data=None):
    org_data          = org_data or {}
    practitioner_data = practitioner_data or {}
    org_id            = org_data.get('org_id', 'unknown-org')
    practitioner_id   = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    practitioner_name = practitioner_data.get('name', 'Unknown Practitioner')
    role_id           = practitioner_data.get('role_id', 'KP-Lab-Tech-Role')
    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"KP Patient {sample_id}"

    return {
        "resourceType": "ServiceRequest",
        "id": f"{sample_id}-service-request",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/ServiceRequest"]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/servicerequest/{org_id}",
                "value": f"SR-KP-{sample_id}"
            }
        ],
        "status": "active",
        "intent": "original-order",
        "priority": "routine",
        "category": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": "108252007",
                        "display": "Laboratory procedure"
                    }
                ]
            }
        ],
        "code": {
            "coding": [
                {
                    "system": "http://loinc.org",
                    "code": "69548-6",
                    "display": "Genetic variant assessment"
                }
            ],
            "text": "Klebsiella pneumoniae Genomic Analysis"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "KP Testing Encounter"
        },
        "occurrenceDateTime": datetime.now(timezone.utc).isoformat(),
        "requester": {
            "reference": f"Practitioner/{practitioner_id}",
            "display": practitioner_name
        },
        "performer": [
            {
                "reference": f"PractitionerRole/{role_id}",
                "display": "Laboratory Technician"
            }
        ],
        "reasonCode": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": "56415008",
                        "display": "Klebsiella pneumoniae"
                    }
                ],
                "text": "Antimicrobial resistance screening for Klebsiella pneumoniae"
            }
        ]
    }

def main():
    parser = argparse.ArgumentParser(description='Merge KP genomic data with clinical metadata into comprehensive FHIR bundle')
    parser.add_argument('--input', required=True, help='Path to input FHIR bundle')
    parser.add_argument('--output', required=True, help='Path to output merged FHIR bundle')
    parser.add_argument('--patient_metadata',      help='Path to patient clinical metadata CSV/Excel file')
    parser.add_argument('--organization_metadata', help='Path to organization metadata CSV/Excel file')
    parser.add_argument('--practitioner_metadata', help='Path to practitioner metadata CSV/Excel file')
    args = parser.parse_args()

    clinical_data = {}
    if args.patient_metadata and os.path.exists(args.patient_metadata):
        clinical_data = load_clinical_metadata(args.patient_metadata)
    else:
        debug_print(f"Patient metadata file not found or not provided: {args.patient_metadata}")

    org_data = {}
    if args.organization_metadata and os.path.exists(args.organization_metadata):
        org_data = load_organization_metadata(args.organization_metadata)
    else:
        debug_print(f"Organization metadata file not found or not provided: {args.organization_metadata}")

    practitioner_data = {}
    if args.practitioner_metadata and os.path.exists(args.practitioner_metadata):
        practitioner_data = load_practitioner_metadata(args.practitioner_metadata)
    else:
        debug_print(f"Practitioner metadata file not found or not provided: {args.practitioner_metadata}")

    try:
        with open(args.input, 'r') as f:
            fhir_bundle = json.load(f)

        sample_ids = set()
        all_observations = []
        
        for entry in fhir_bundle.get('entry', []):
            resource = entry.get('resource', {})
            if resource.get('resourceType') == 'Observation':
                all_observations.append(resource)
                subject_ref = resource.get('subject', {}).get('reference', '')
                
                if subject_ref.startswith('Patient/'):
                    sample_id = subject_ref.replace('Patient/', '').replace('-patient', '')
                    sample_ids.add(sample_id)

        if not sample_ids:
            filename = os.path.basename(args.input)
            filename_sample_id = filename.replace('.fhir.json', '').replace('_ont', '').replace('_illumina', '')
            sample_ids.add(filename_sample_id)

        matched_samples = {}
        for sample_id in sample_ids:
            sample_clinical_data = find_matching_sample(sample_id, clinical_data)
            matched_samples[sample_id] = sample_clinical_data
            if sample_clinical_data:
                debug_print(f"Found clinical data for sample: {sample_id}")
            else:
                debug_print(f"No clinical data found for sample: {sample_id}")

        merged_bundle = {
            "resourceType": "Bundle",
            "id": str(uuid.uuid4()),
            "meta": {
                "lastUpdated": datetime.now(timezone.utc).isoformat(),
                "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Bundle"]
            },
            "type": "transaction", 
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "entry": []
        }

        org_resource = create_organization_resource(org_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": org_resource,
            "request": {
                "method": "PUT",
                "url": f"Organization/{org_resource['id']}"
            }
        })

        practitioner_resource = create_practitioner_resource(practitioner_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": practitioner_resource,
            "request": {
                "method": "PUT",
                "url": f"Practitioner/{practitioner_resource['id']}"
            }
        })

        role_resource = create_practitioner_role_resource(practitioner_data, org_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": role_resource,
            "request": {
                "method": "PUT",
                "url": f"PractitionerRole/{role_resource['id']}"
            }
        })

        for sample_id, sample_clinical_data in matched_samples.items():
            debug_print(f"Adding patient resources for sample: {sample_id}")

            patient_resource = create_patient_resource(sample_id, sample_clinical_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": patient_resource,
                "request": {
                    "method": "PUT",
                    "url": f"Patient/{patient_resource['id']}"
                }
            })

            specimen_resource = create_specimen_resource(sample_id, sample_clinical_data, practitioner_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": specimen_resource,
                "request": {
                    "method": "PUT",
                    "url": f"Specimen/{specimen_resource['id']}"
                }
            })

            service_request_resource = create_service_request_resource(sample_id, sample_clinical_data, practitioner_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": service_request_resource,
                "request": {
                    "method": "PUT",
                    "url": f"ServiceRequest/{service_request_resource['id']}"
                }
            })

        observations_by_sample = {}
        for obs in all_observations:
            subject_ref = obs.get('subject', {}).get('reference', '')
            if subject_ref.startswith('Patient/'):
                sample_id = subject_ref.replace('Patient/', '').replace('-patient', '')
                if sample_id not in observations_by_sample:
                    observations_by_sample[sample_id] = []
                observations_by_sample[sample_id].append(obs)

        for sample_id, sample_observations in observations_by_sample.items():
            sample_clinical_data = matched_samples.get(sample_id)
            
            if sample_observations:
                diagnostic_report = create_diagnostic_report(
                    sample_id, 
                    sample_observations, 
                    sample_clinical_data,
                    org_data,
                    practitioner_data
                )
                
                merged_bundle['entry'].append({
                    "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                    "resource": diagnostic_report,
                    "request": {
                        "method": "PUT",
                        "url": f"DiagnosticReport/{diagnostic_report['id']}"
                    }
                })

        for entry in fhir_bundle.get('entry', []):
            resource = entry.get('resource', {})
            resource_type = resource.get('resourceType')
            resource_id = resource.get('id')
            
            entry_with_request = {
                "fullUrl": entry.get('fullUrl', f"urn:uuid:{str(uuid.uuid4())}"),
                "resource": resource,
                "request": {
                    "method": "PUT",
                    "url": f"{resource_type}/{resource_id}" if resource_id else f"{resource_type}"
                }
            }
            merged_bundle['entry'].append(entry_with_request)

        with open(args.output, 'w') as f:
            json.dump(merged_bundle, f, indent=2)

    except Exception as e:
        import traceback
        sys.exit(1)

if __name__ == "__main__":
    main()