"""Shared Pfam prioritisation rules for Phorific hotspot annotation."""

from __future__ import annotations


def excluded_reason(hit: dict[str, str]) -> str:
    """Return why a Pfam hit is uninformative for accessory-biology review."""
    text = f"{hit.get('name', '')} {hit.get('description', '')}".lower()
    groups = {
        "integration_or_transposition": (
            "integrase", "recombinase", "resolvase", "transpos", "tnp", "rve",
            "is605", "isl3", "insertion sequence", "retrotrans", "mule",
        ),
        "generic_regulation_or_dna_binding": (
            "helix-turn-helix", "hth", "trans_reg", "response_reg", "laci",
            "gntr", "deorc", "sigma54_activat", "cro/c1", "winged helix",
            "hatpase", "hiska",
        ),
        "phage_structure_or_lysis": (
            "capsid", "phage_tail", "prophage_tail", "sipho_tail", "terminase",
            "portal", "tape_measure", "holin", "phage_gp", "phage_g",
        ),
    }
    for reason, markers in groups.items():
        if any(marker in text for marker in markers):
            return reason
    return ""


def interesting_category(hit: dict[str, str]) -> str:
    """Prioritise residual domains with plausible accessory biological value."""
    accession = hit.get("accession", "").split(".")[0]
    text = f"{hit.get('name', '')} {hit.get('description', '')}".lower()
    explicit = {
        "PF07275": "anti-restriction",
        "PF12648": "conjugation/secretion",
        "PF02604": "toxin-antitoxin",
        "PF02452": "toxin-antitoxin",
        "PF05016": "toxin-antitoxin",
        "PF15781": "toxin-antitoxin",
        "PF02661": "toxin-antitoxin",
        "PF20613": "toxin-antitoxin",
        "PF13274": "defence-like",
        "PF07282": "defence-like",
        "PF07751": "defence-like",
        "PF04326": "defence-like",
        "PF18813": "defence-like",
        "PF09607": "defence-like",
        "PF02534": "conjugation/secretion",
        "PF17802": "adhesion/cell-surface",
        "PF04203": "adhesion/cell-surface",
        "PF16555": "adhesion/cell-surface",
        "PF13614": "ATPase/helicase/nuclease",
        "PF00004": "ATPase/helicase/nuclease",
        "PF07728": "ATPase/helicase/nuclease",
    }
    if accession in explicit:
        return explicit[accession]
    rules = (
        ("anti-restriction", ("arda", "ardb", "anti-restriction", "antirestriction")),
        ("toxin-antitoxin", (
            "toxin", "antitoxin", "higb", "hipa", "hipb", "mazef", "mazf",
            "maze", "rele", "relb", "vapb", "vapc", "pemk", "yefm", "phd",
            "pare", "doc", "fic",
        )),
        ("defence-like", (
            "abortive", "abi_", "crispr", "cas", "restriction", "cbass", "gabija",
            "thoeris", "zorya", "brex", "brx", "dnd", "panacea", "pbecr",
            "schlafen", "alba", "brkdbd",
        )),
        ("conjugation/secretion", (
            "tcpe", "tcp", "t4ss", "t7ss", "type_iv", "type_vii", "secretion",
            "conjug", "relaxase", "coupling", "virb", "traa", "trab", "mobl",
        )),
        ("adhesion/cell-surface", (
            "adhes", "pilin", "sortase", "spaa", "lpxtg", "collagen",
        )),
        ("ATPase/helicase/nuclease", (
            "nuclease", "endonuclease", "ribonuclease", "rnase", "helicase",
            "aaa", "atpase", "hepn", "hnh", "giy-yig", "pin_",
        )),
    )
    for category, markers in rules:
        if any(marker in text for marker in markers):
            return category
    return ""
