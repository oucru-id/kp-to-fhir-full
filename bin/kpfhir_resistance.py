
import re

from kpfhir_terminology import COLISTIN_MUTATION_GENES, DRUG_CLASS_ALIASES, SNOMED

PANEL_DRUGS = [
    "Ampicillin", "Amoxicillin-Clavulanate", "Piperacillin-Tazobactam", "Cefotaxime", "Ceftazidime",
    "Cefepime", "Meropenem", "Ciprofloxacin", "Gentamicin", "Amikacin", "Trimethoprim-Sulfamethoxazole",
    "Colistin",
]

AMINOGLYCOSIDE_PROFILES = {
    'aac(3)-i': {'Gentamicin': True, 'Amikacin': False},
    'aac(3)-ia': {'Gentamicin': True, 'Amikacin': False},
    'aac(3)-ib': {'Gentamicin': True, 'Amikacin': False},
    'aac(3)-ii': {'Gentamicin': True, 'Amikacin': False},
    'aac(3)-iia': {'Gentamicin': True, 'Amikacin': False},
    'aac(3)-iib': {'Gentamicin': True, 'Amikacin': False},
    'aac(3)-iii': {'Gentamicin': True, 'Amikacin': False},
    'aac(3)-iv': {'Gentamicin': True, 'Amikacin': True},
    'aac(3)-vi': {'Gentamicin': True, 'Amikacin': False},
    "aac(6')-i": {'Gentamicin': False, 'Amikacin': True},
    "aac(6')-ia": {'Gentamicin': False, 'Amikacin': True},
    "aac(6')-ib": {'Gentamicin': False, 'Amikacin': True},
    "aac(6')-ib-cr": {'Gentamicin': False, 'Amikacin': True},
    "aac(6')-ii": {'Gentamicin': True, 'Amikacin': True},
    "aac(6')-30": {'Gentamicin': False, 'Amikacin': True},
    "ant(2'')-i": {'Gentamicin': True, 'Amikacin': False},
    "ant(2'')-ia": {'Gentamicin': True, 'Amikacin': False},
    "ant(3'')-i": {'Gentamicin': False, 'Amikacin': False},
    "ant(3'')-ia": {'Gentamicin': False, 'Amikacin': False},
    "ant(4')-i": {'Gentamicin': False, 'Amikacin': True},
    "ant(4')-ia": {'Gentamicin': False, 'Amikacin': True},
    "ant(4')-ii": {'Gentamicin': False, 'Amikacin': True},
    "ant(4')-iia": {'Gentamicin': False, 'Amikacin': True},
    "aph(3')-i": {'Gentamicin': False, 'Amikacin': False},
    "aph(3')-ia": {'Gentamicin': False, 'Amikacin': False},
    "aph(3')-ii": {'Gentamicin': False, 'Amikacin': False},
    "aph(3')-iia": {'Gentamicin': False, 'Amikacin': False},
    "aph(3')-vi": {'Gentamicin': False, 'Amikacin': True},
    "aph(3')-via": {'Gentamicin': False, 'Amikacin': True},
    "aph(3'')-i": {'Gentamicin': False, 'Amikacin': False},
    "aph(3'')-ib": {'Gentamicin': False, 'Amikacin': False},
    'aph(6)-i': {'Gentamicin': False, 'Amikacin': False},
    'aph(6)-ia': {'Gentamicin': False, 'Amikacin': False},
    'aph(6)-id': {'Gentamicin': False, 'Amikacin': False},
    'arma': {'Gentamicin': True, 'Amikacin': True},
    'rmta': {'Gentamicin': True, 'Amikacin': True},
    'rmtb': {'Gentamicin': True, 'Amikacin': True},
    'rmtc': {'Gentamicin': True, 'Amikacin': True},
    'rmtd': {'Gentamicin': True, 'Amikacin': True},
    'rmte': {'Gentamicin': True, 'Amikacin': True},
    'rmtf': {'Gentamicin': True, 'Amikacin': True},
    'rmtg': {'Gentamicin': True, 'Amikacin': True},
    'rmth': {'Gentamicin': True, 'Amikacin': True},
    'npma': {'Gentamicin': True, 'Amikacin': True},
}


def normalise_aminoglycoside_name(name):
    for char in "()'^*?":
        name = name.replace(char, '')
    return name.lower().replace('_', '-').replace(' ', '-')


_AMINOGLYCOSIDE_LOOKUP = sorted(
    ((normalise_aminoglycoside_name(name), profile) for name, profile in AMINOGLYCOSIDE_PROFILES.items()),
    key=lambda item: len(item[0]), reverse=True)


def predict_panel(typing_data, log=None):
    drug_status = {drug: "Susceptible" for drug in PANEL_DRUGS}
    drug_status["Ampicillin"] = "Resistant"

    resistance = typing_data.get('resistance', {})
    has_carb = resistance.get('carbapenemase', 'None') != 'None'
    has_esbl = resistance.get('esbl', 'None') != 'None'

    has_ompk35 = has_ompk36 = False
    has_flq_gene = has_flq_mutation = False
    has_col_gene = has_col_mutation = False
    has_trimethoprim = has_sulfonamide = False

    for gene in resistance.get('detailed_resistance', []):
        gene_name = gene.get('gene', '').lower().replace('^', '').strip()
        drug_class = gene.get('drug_class', '').lower()
        drug_class = DRUG_CLASS_ALIASES.get(drug_class, drug_class)

        if 'aminoglycoside' in drug_class:
            lookup = normalise_aminoglycoside_name(gene_name)
            profile = next((p for name, p in _AMINOGLYCOSIDE_LOOKUP if lookup.startswith(name)), None)
            if profile is not None:
                if profile.get('Gentamicin'):
                    drug_status["Gentamicin"] = "Resistant"
                if profile.get('Amikacin'):
                    drug_status["Amikacin"] = "Resistant"
            elif 'aac' in lookup or 'ant' in lookup or 'aph' in lookup:
                drug_status["Gentamicin"] = "Resistant"
                if log:
                    log(f"Unknown aminoglycoside gene {gene_name} - marking Gentamicin resistant (conservative)")

        if drug_class == 'fluoroquinolone' or 'qnr' in gene_name:
            has_flq_gene = True
        if "aac(6')-ib-cr" in gene_name or 'aac(6)-ib-cr' in gene_name:
            has_flq_gene = True
            drug_status["Amikacin"] = "Resistant"
        if 'trimethoprim' in drug_class or gene_name.startswith('dfr'):
            has_trimethoprim = True
        if 'sulfonamide' in drug_class or gene_name.startswith('sul'):
            has_sulfonamide = True
        if 'colistin' in drug_class or gene_name.startswith('mcr'):
            has_col_gene = True

    for mut in resistance.get('point_mutations', []):
        gene = mut.get('gene', '')
        if gene in ['GyrA', 'ParC']:
            has_flq_mutation = True
        if gene == 'OmpK35':
            has_ompk35 = True
        if gene == 'OmpK36':
            has_ompk36 = True
        if gene in COLISTIN_MUTATION_GENES:
            has_col_mutation = True

    if has_carb:
        for drug in ("Meropenem", "Amoxicillin-Clavulanate", "Piperacillin-Tazobactam", "Cefotaxime",
                     "Ceftazidime", "Cefepime"):
            drug_status[drug] = "Resistant"
    if has_esbl:
        for drug in ("Cefotaxime", "Ceftazidime", "Cefepime", "Piperacillin-Tazobactam"):
            drug_status[drug] = "Resistant"
    if has_ompk35 and has_ompk36:
        if has_esbl:
            drug_status["Meropenem"] = "Resistant"
            drug_status["Piperacillin-Tazobactam"] = "Resistant"
    elif has_ompk36 and has_esbl:
        drug_status["Piperacillin-Tazobactam"] = "Resistant"
    if has_flq_gene or has_flq_mutation:
        drug_status["Ciprofloxacin"] = "Resistant"
    if has_col_gene or has_col_mutation:
        drug_status["Colistin"] = "Resistant"
    if has_trimethoprim or has_sulfonamide:
        drug_status["Trimethoprim-Sulfamethoxazole"] = "Resistant"
    return drug_status


DRUG_CLASS_NORMALIZATION = {
    'beta-lactam': 'beta-lactam', 'beta_lactam': 'beta-lactam', 'betalactam': 'beta-lactam',
    'aminoglycoside': 'aminoglycoside', 'fluoroquinolone': 'fluoroquinolone', 'quinolone': 'fluoroquinolone',
    'macrolide': 'macrolide', 'sulfonamide': 'sulfonamide', 'trimethoprim': 'trimethoprim',
    'tetracycline': 'tetracycline', 'chloramphenicol': 'chloramphenicol', 'phenicol': 'chloramphenicol',
    'colistin': 'colistin', 'polymyxin': 'colistin', 'carbapenem': 'carbapenem', 'fosfomycin': 'fosfomycin',
    'rifampicin': 'rifampicin', 'bla': 'beta-lactam', 'agly': 'aminoglycoside', 'flq': 'fluoroquinolone',
    'mls': 'macrolide', 'mlsb': 'macrolide', 'sul': 'sulfonamide', 'tmt': 'trimethoprim', 'tet': 'tetracycline',
    'col': 'colistin', 'fos': 'fosfomycin', 'rif': 'rifampicin', 'dfr': 'trimethoprim',
}

NARROW_SPECTRUM_SHV = {f'shv-{n}' for n in (
    1, 11, 26, 28, 33, 36, 38, 41, 42, 60, 62, 65, 75, 76, 77, 79, 99, 100, 108, 119, 155,
    *range(164, 201), 206, 207)}

ESBL_SHV = {f'shv-{n}' for n in (
    *range(2, 10), 12, *range(13, 19), 24, 25, 27, *range(29, 33), 34, 35, 37, 39, 40, 43, 44, 45, 46, 48, 49,
    52, 55, 56, 57, 59, 61, 63, 64, *range(66, 75), 78, *range(80, 99), 101, *range(102, 108),
    *range(109, 119), *range(120, 155), *range(156, 164))}

ESBL_TEM = {f'tem-{n}' for n in (
    *range(3, 11), 12, 15, 16, 17, 19, 21, 22, 24, 25, 26, 27, 28, 29, 42, 43, 46, 47, 49, 50, 52, 53, 60,
    63, 64, 70, 71, 87, 88, 89, 90, 92, 109, 116, 121, 126, 127, 128, 129, 131, 132, 134, *range(135, 145),
    146, *range(147, 165), *range(167, 173), 176, *range(178, 201))}

AMPC_PATTERNS = ['dha', 'cmy', 'mox', 'fox', 'acc', 'act', 'mir', 'lat', 'bil', 'cfea']
CARBAPENEMASE_PATTERNS = ['kpc', 'ndm', 'oxa-48', 'oxa-181', 'oxa-232', 'vim', 'imp', 'ges']


def normalize_drug_class(drug_class):
    if not drug_class:
        return None
    lower = drug_class.lower().strip()
    if lower in DRUG_CLASS_NORMALIZATION:
        return DRUG_CLASS_NORMALIZATION[lower]
    underscored = lower.replace('-', '_').replace(' ', '_')
    if underscored in DRUG_CLASS_NORMALIZATION:
        return DRUG_CLASS_NORMALIZATION[underscored]
    for key, normalized in DRUG_CLASS_NORMALIZATION.items():
        if len(key) >= 3 and (key in lower or lower in key):
            return normalized
    return None


def normalize_gene_name(gene_name):
    if not gene_name:
        return ''
    normalized = gene_name.lower().replace('^', '').replace('*', '').replace('.', '-').strip()
    normalized = re.sub(r'-v\d+$', '', normalized)
    return re.sub(r'\.v\d+$', '', normalized)


def _in_family(normalized, family):
    return normalized in family or normalized in {g.replace('-', '') for g in family}


def is_narrow_spectrum_shv(gene):
    return _in_family(normalize_gene_name(gene), NARROW_SPECTRUM_SHV)


def is_esbl_gene(gene):
    normalized = normalize_gene_name(gene)
    if is_narrow_spectrum_shv(normalized):
        return False
    if 'ctx-m' in normalized or 'ctxm' in normalized:
        return True
    return _in_family(normalized, ESBL_SHV) or _in_family(normalized, ESBL_TEM)


def _joined(values):
    return ', '.join(dict.fromkeys(values))


def classify_resistance_profile(gene_observations, panel_observation, log=None):
    carbapenemases, esbl_genes, ampc_genes, colistin, other = [], [], [], [], []
    classes = set()
    resistant_drugs = set()

    if panel_observation:
        for comp in panel_observation.get('component', []):
            drug_text = comp.get('code', {}).get('text', '')
            drug_name = drug_text.split(' - ')[0] if ' - ' in drug_text else drug_text
            if comp.get('valueCodeableConcept', {}).get('coding', [{}])[0].get('code') == 'LA6676-6':
                resistant_drugs.add(drug_name)

    for obs in gene_observations:
        div_text = obs.get('text', {}).get('div', '')
        is_intrinsic = '(chromosomal)' in div_text or 'intrinsic' in div_text
        gene = drug_class = None
        for component in obs.get('component', []):
            if 'valueCodeableConcept' not in component:
                continue
            code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '')
            value_text = component['valueCodeableConcept'].get('text', '')
            if 'Gene studied' in code_display:
                gene = value_text
            elif 'drug efficacy' in code_display:
                drug_class = value_text
        if not gene:
            continue

        gene_display = gene.upper().replace('^', '').replace('*', '').strip()
        gene_normalized = normalize_gene_name(gene)
        is_carb = any(c in gene_normalized for c in CARBAPENEMASE_PATTERNS)
        is_esbl = is_esbl_gene(gene_normalized)
        is_ampc = any(a in gene_normalized for a in AMPC_PATTERNS)
        is_col = 'mcr' in gene_normalized or bool(drug_class and 'colistin' in str(drug_class).lower())

        if is_narrow_spectrum_shv(gene_normalized):
            continue
        if is_intrinsic and not (is_carb or is_esbl or is_ampc or is_col):
            continue

        normalized_class = normalize_drug_class(drug_class)
        if normalized_class:
            classes.add(normalized_class)
        if is_carb:
            carbapenemases.append(gene_display)
            classes.update({'carbapenem', 'beta-lactam'})
        elif is_esbl:
            esbl_genes.append(gene_display)
            classes.add('beta-lactam')
        elif is_ampc:
            ampc_genes.append(gene_display)
            classes.add('beta-lactam')
        elif is_col:
            colistin.append(gene_display)
            classes.add('colistin')
        else:
            other.append(gene_display)

    if log:
        log(f"Total affected drug classes: {len(classes)} - {sorted(classes)}")

    drugs_text = ""
    acquired = [d for d in sorted(resistant_drugs) if d != "Ampicillin"]
    if acquired:
        drugs_text = f" Predicted resistance to: {', '.join(acquired)}."

    n_classes = len(classes)
    all_genes = sorted(set(carbapenemases + esbl_genes + ampc_genes + colistin + other))
    genes_text = f"Resistance genes: {', '.join(all_genes)}." if all_genes else ""

    if carbapenemases:
        if colistin:
            return "Extensively Drug-Resistant (XDR)", (
                f"Carbapenem and colistin resistant K. pneumoniae. "
                f"Carbapenemases: {_joined(carbapenemases)}. "
                f"Colistin resistance: {_joined(colistin)}."
                f"{drugs_text}")
        return "Carbapenem-Resistant Enterobacteriaceae (CRE)", (
            f"Carbapenem-resistant K. pneumoniae (CP-CRE). "
            f"Carbapenemases detected: {_joined(carbapenemases)}."
            f"{drugs_text}")
    if esbl_genes:
        if n_classes >= 3:
            return "Multidrug-Resistant ESBL", (
                f"ESBL-producing K. pneumoniae with additional resistance (MDR). "
                f"ESBL genes: {_joined(esbl_genes)}. "
                f"Resistance affects {n_classes} drug classes. "
                f"{genes_text}"
                f"{drugs_text}")
        return "Extended-Spectrum Beta-Lactamase (ESBL)", (
            f"ESBL-producing K. pneumoniae. "
            f"ESBL genes: {_joined(esbl_genes)}."
            f"{drugs_text}")
    if ampc_genes:
        if n_classes >= 3:
            return "Multidrug-Resistant AmpC", (
                f"AmpC beta-lactamase producing K. pneumoniae with additional resistance (MDR). "
                f"AmpC genes: {_joined(ampc_genes)}. "
                f"Resistance affects {n_classes} drug classes. "
                f"{genes_text}"
                f"{drugs_text}")
        return "AmpC Beta-Lactamase (AmpC)", (
            f"AmpC beta-lactamase producing K. pneumoniae. "
            f"AmpC genes: {_joined(ampc_genes)}."
            f"{drugs_text}")
    if n_classes >= 3:
        return "Multidrug-Resistant (MDR)", (
            f"Multiple resistance genes detected affecting {n_classes} drug classes. "
            f"{genes_text}"
            f"{drugs_text}")
    if other or n_classes > 0:
        return "Resistant", (
            f"Resistance genes detected: {_joined(other) if other else 'See panel'}. "
            f"{drugs_text}")
    return "Susceptible", "No acquired resistance genes detected. Intrinsic ampicillin resistance only."


CLASSIFICATION_KEYS = {
    "Extensively Drug-Resistant (XDR)": "XDR",
    "Carbapenem-Resistant Enterobacteriaceae (CRE)": "CRE",
    "Multidrug-Resistant ESBL": "Multidrug-Resistant ESBL",
    "Extended-Spectrum Beta-Lactamase (ESBL)": "ESBL",
    "Multidrug-Resistant AmpC": "Multidrug-Resistant AmpC",
    "AmpC Beta-Lactamase (AmpC)": "AmpC",
    "Multidrug-Resistant (MDR)": "MDR",
    "Resistant": "Resistant",
    "Susceptible": "Susceptible",
}

SNOMED_CONCLUSIONS = {
    "CRE": {"system": SNOMED, "code": "1098201000112108",
            "display": "Carbapenemase-producing Klebsiella pneumoniae (organism)"},
    "ESBL": {"system": SNOMED, "code": "409801009",
             "display": "Extended spectrum beta-lactamase producing Klebsiella pneumoniae (organism)"},
    "AmpC": {"system": SNOMED, "code": "1098101000112102",
             "display": "AmpC beta-lactamase producing Klebsiella pneumoniae (organism)"},
    "MDR": {"system": SNOMED, "code": "714315002",
            "display": "Multidrug-resistant Klebsiella pneumoniae (organism)"},
}

KEMKES_CLINICAL_TERM = "http://terminology.kemkes.go.id/CodeSystem/clinical-term"
KEMKES_CONCLUSIONS = {
    "XDR": {"system": KEMKES_CLINICAL_TERM, "code": "SP000681",
            "display": "Extensively drug resistant Klebsiella pneumoniae"},
    "Multidrug-Resistant ESBL": {"system": KEMKES_CLINICAL_TERM, "code": "SP000687",
                                 "display": "Multidrug-resistant ESBL-producing Klebsiella pneumoniae"},
    "Multidrug-Resistant AmpC": {"system": KEMKES_CLINICAL_TERM, "code": "SP000683",
                                 "display": "Multidrug-resistant AmpC-producing Klebsiella pneumoniae"},
    "Susceptible": {"system": KEMKES_CLINICAL_TERM, "code": "KP-SO",
                    "display": "Klebsiella pneumoniae Sensitif Obat"},
    "Resistant": {"system": KEMKES_CLINICAL_TERM, "code": "KP-RL",
                  "display": "Klebsiella pneumoniae Resisten Lain"},
}

FALLBACK_ORDER = ["Multidrug-Resistant ESBL", "Multidrug-Resistant AmpC", "XDR", "CRE", "ESBL", "AmpC", "MDR"]


def conclusion_key(resistance_class):
    if resistance_class in CLASSIFICATION_KEYS:
        return CLASSIFICATION_KEYS[resistance_class]
    known = set(CLASSIFICATION_KEYS.values())
    if resistance_class in known:
        return resistance_class
    for key in FALLBACK_ORDER:
        if key in resistance_class:
            return key
    if "Susceptible" in resistance_class:
        return "Susceptible"
    if "Resistant" in resistance_class:
        return "Resistant"
    return None


def conclusion_coding(resistance_class, include_kemkes):
    key = conclusion_key(resistance_class)
    if key is None:
        return None
    if key in SNOMED_CONCLUSIONS:
        return dict(SNOMED_CONCLUSIONS[key])
    if include_kemkes and key in KEMKES_CONCLUSIONS:
        return dict(KEMKES_CONCLUSIONS[key])
    return None
