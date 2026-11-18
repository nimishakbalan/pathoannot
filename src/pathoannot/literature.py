"""Literature evidence from Europe PMC: is the species mostly discussed as a disease agent,
a commensal / microbiome member, or an environmental organism?"""
import time

import requests

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
DISEASE_TERMS = ('infection OR infections OR pathogen OR pathogenic OR disease OR outbreak OR '
                 '"case report" OR virulence OR sepsis OR pneumonia OR meningitis OR bacteremia')
COMMENSAL_TERMS = ('commensal OR "normal flora" OR probiotic OR "healthy volunteers" OR '
                   '"healthy individuals" OR "healthy subjects"')
ENV_TERMS = ('soil OR sediment OR rhizosphere OR marine OR freshwater OR wastewater OR '
             '"environmental isolate" OR "environmental samples"')
MIN_PAPERS = 20


def europepmc_count(query):
    r = requests.get(EPMC, params={"query": query, "format": "json", "pageSize": 1,
                                   "resultType": "lite"}, timeout=60)
    r.raise_for_status()
    return int(r.json().get("hitCount", 0))


def lit_summary(species):
    name = species.replace('"', "").replace("[", "").replace("]", "").strip()
    base = f'TITLE_ABS:"{name}"'
    total = europepmc_count(base)
    time.sleep(0.15)
    if total == 0:
        return dict(lit_total=0, lit_disease_frac=None, lit_commensal_frac=None, lit_environment_frac=None)
    out = {}
    for key, terms in (("disease", DISEASE_TERMS), ("commensal", COMMENSAL_TERMS), ("environment", ENV_TERMS)):
        n = europepmc_count(f"{base} AND TITLE_ABS:({terms})")
        out[f"lit_{key}_frac"] = round(n / total, 3)
        time.sleep(0.15)
    return dict(lit_total=total, **out)


def lit_label(s, min_papers=MIN_PAPERS):
    if not s or s.get("lit_total", 0) < min_papers:
        return "Insufficient literature"
    d, c, e = s["lit_disease_frac"], s["lit_commensal_frac"], s["lit_environment_frac"]
    if d >= 0.5:
        return "Disease-dominated literature"
    if e >= 0.4 and d < 0.25:
        return "Environmental-dominated literature"
    if c >= 0.15 and d < 0.35:
        return "Commensal/microbiome-dominated literature"
    return "Mixed literature"
