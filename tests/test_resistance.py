
import pytest

from kpfhir_resistance import classify_resistance_profile, conclusion_coding, predict_panel


def gene_obs(gene, drug_class, intrinsic=False):
    div = f"Resistance gene {gene}" + (" (intrinsic chromosomal)" if intrinsic else "")
    return {"text": {"div": div}, "component": [
        {"code": {"coding": [{"display": "Gene studied [ID]"}]}, "valueCodeableConcept": {"text": gene}},
        {"code": {"coding": [{"display": "Genetic variation's effect on drug efficacy"}]},
         "valueCodeableConcept": {"text": drug_class}}]}


CASES = [
    ([], "Susceptible"),
    ([("SHV-11", "beta-lactam", True)], "Susceptible"),
    ([("SHV-1", "beta-lactam", False)], "Susceptible"),
    ([("KPC-2", "beta-lactam", False)], "Carbapenem-Resistant Enterobacteriaceae (CRE)"),
    ([("NDM-1", "beta-lactam", False), ("mcr-1.1", "col", False)], "Extensively Drug-Resistant (XDR)"),
    ([("KPC-3", "beta-lactam", False), ("MgrB -62%", "colistin", False)], "Extensively Drug-Resistant (XDR)"),
    ([("CTX-M-15", "beta-lactam", False)], "Extended-Spectrum Beta-Lactamase (ESBL)"),
    ([("SHV-12", "beta-lactam", False)], "Extended-Spectrum Beta-Lactamase (ESBL)"),
    ([("CTX-M-15", "beta-lactam", False), ("sul2", "sul", False), ("tet(A).v1", "tet", False)],
     "Multidrug-Resistant ESBL"),
    ([("DHA-1", "beta-lactam", False)], "AmpC Beta-Lactamase (AmpC)"),
    ([("DHA-1", "beta-lactam", False), ("sul1", "sul", False), ("dfrA12", "tmt", False)],
     "Multidrug-Resistant AmpC"),
    ([("sul1", "sul", False), ("dfrA12", "tmt", False), ("arr-3", "rif", False)], "Multidrug-Resistant (MDR)"),
    ([("sul1", "sul", False)], "Resistant"),
    ([("OXA-1", "beta-lactam", True)], "Susceptible"),
]


@pytest.mark.parametrize("genes,expected", CASES)
def test_classification_table(genes, expected):
    resistance_class, description = classify_resistance_profile([gene_obs(*g) for g in genes], None)
    assert resistance_class == expected
    assert description


def test_conclusion_text_is_deterministic():
    observations = [gene_obs(g, "beta-lactam") for g in ("KPC-2", "NDM-1", "OXA-48", "VIM-1")]
    texts = {classify_resistance_profile(observations, None)[1] for _ in range(20)}
    assert texts == {"Carbapenem-resistant K. pneumoniae (CP-CRE). "
                     "Carbapenemases detected: KPC-2, NDM-1, OXA-48, VIM-1."}


@pytest.mark.parametrize("resistance_class,generic,satusehat", [
    ("Carbapenem-Resistant Enterobacteriaceae (CRE)", "1098201000112108", "1098201000112108"),
    ("Extensively Drug-Resistant (XDR)", None, "SP000681"),
    ("Multidrug-Resistant ESBL", None, "SP000687"),
    ("Extended-Spectrum Beta-Lactamase (ESBL)", "409801009", "409801009"),
    ("Multidrug-Resistant AmpC", None, "SP000683"),
    ("AmpC Beta-Lactamase (AmpC)", "1098101000112102", "1098101000112102"),
    ("Multidrug-Resistant (MDR)", "714315002", "714315002"),
    ("Resistant", None, "KP-RL"),
    ("Susceptible", None, "KP-SO"),
])
def test_conclusion_codings_per_profile(resistance_class, generic, satusehat):
    assert (conclusion_coding(resistance_class, include_kemkes=False) or {}).get("code") == generic
    assert (conclusion_coding(resistance_class, include_kemkes=True) or {}).get("code") == satusehat


def typing(genes=(), mutations=(), esbl='None', carbapenemase='None'):
    return {"resistance": {"esbl": esbl, "carbapenemase": carbapenemase,
                           "detailed_resistance": [{"gene": g, "drug_class": c} for g, c in genes],
                           "point_mutations": [{"gene": g, "mutation": m} for g, m in mutations]}}


def resistant(panel):
    return sorted(d for d, s in panel.items() if s == "Resistant")


@pytest.mark.parametrize("data,expected", [
    (typing(), ["Ampicillin"]),
    (typing([("rmtB.v1", "agly")]), ["Amikacin", "Ampicillin", "Gentamicin"]),
    (typing([("aac(3)-IIa.v1^", "agly")]), ["Ampicillin", "Gentamicin"]),
    (typing([("aph(3')-Ia.v1", "agly")]), ["Ampicillin"]),
    (typing([("aac(6')-Ib-cr.v2", "flq")]), ["Amikacin", "Ampicillin", "Ciprofloxacin"]),
    (typing(mutations=[("GyrA", "p.Ser83Ile")]), ["Ampicillin", "Ciprofloxacin"]),
    (typing(mutations=[("MgrB", "-62%")]), ["Ampicillin", "Colistin"]),
    (typing([("sul1", "sul")]), ["Ampicillin", "Trimethoprim-Sulfamethoxazole"]),
    (typing(esbl="CTX-M-15"), ["Ampicillin", "Cefepime", "Cefotaxime", "Ceftazidime", "Piperacillin-Tazobactam"]),
    (typing(mutations=[("OmpK35", "-25%"), ("OmpK36", "GD")], esbl="CTX-M-15"),
     ["Ampicillin", "Cefepime", "Cefotaxime", "Ceftazidime", "Meropenem", "Piperacillin-Tazobactam"]),
    (typing(carbapenemase="KPC-2"),
     ["Amoxicillin-Clavulanate", "Ampicillin", "Cefepime", "Cefotaxime", "Ceftazidime", "Meropenem",
      "Piperacillin-Tazobactam"]),
])
def test_panel_prediction_table(data, expected):
    assert resistant(predict_panel(data)) == expected
