"""
==============================================================================
PROJET : Analyse du réseau d'interaction (interactome) des cytokines
         impliquées dans la réponse inflammatoire
==============================================================================

Ce script récupère automatiquement :
1. Les annotations biologiques de chaque cytokine depuis UniProtKB
2. Le réseau d'interactions protéine-protéine depuis STRING-db

Aucune donnée patient n'est nécessaire : toutes les données proviennent
de bases de connaissances publiques et libres d'accès.

Auteur : [Ton nom]
==============================================================================
"""

import requests
import pandas as pd
import time
import os

# ==============================================================================
# ÉTAPE 0 — LISTE DES CYTOKINES ÉTUDIÉES
# ==============================================================================
# Ce panel de 20 cytokines a été choisi car il représente un ensemble
# cohérent et largement étudié dans la littérature sur l'inflammation
# (voir par exemple Fjell et al. 2013, PLOS ONE, panel de 39 cytokines
# dans le sepsis ; et études sur l'asthme/COPD consultées sur ImmPort).
#
# Le panel couvre plusieurs familles fonctionnelles, ce qui rend le
# réseau et l'analyse de communautés plus intéressants :
#
#   - Pro-inflammatoires classiques : IL1B, IL6, TNF, IL8 (CXCL8)
#   - Anti-inflammatoires / régulatrices : IL10, TGFB1, IL1RN
#   - Réponse Th1 : IFNG, IL12A, IL12B, IL2
#   - Réponse Th2 (allergie/asthme) : IL4, IL5, IL13, IL9
#   - Réponse Th17 : IL17A, IL22, IL23A
#   - Chimiokines (recrutement cellulaire) : CCL2, CCL11, CXCL10
#
# Tu peux librement ajouter/retirer des gènes ici selon ton besoin.
# ==============================================================================

CYTOKINES = [
    "IL1B", "IL6", "TNF", "CXCL8",       # pro-inflammatoires
    "IL10", "TGFB1", "IL1RN",            # anti-inflammatoires / régulatrices
    "IFNG", "IL12A", "IL12B", "IL2",     # Th1
    "IL4", "IL5", "IL13", "IL9",         # Th2
    "IL17A", "IL22", "IL23A",            # Th17
    "CCL2", "CCL11", "CXCL10",           # chimiokines
]

ORGANISM = "human"
SPECIES_TAXID = 9606  # Homo sapiens, identifiant utilisé par STRING

OUTPUT_DIR = "data/raw"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ==============================================================================
# ÉTAPE 1 — RÉCUPÉRATION DES ANNOTATIONS UNIPROT
# ==============================================================================

def get_uniprot_info(gene_name, organism="human"):
    """
    Interroge l'API REST d'UniProtKB pour récupérer les informations
    d'une protéine à partir de son nom de gène.

    Retourne un dictionnaire avec : accession, nom de la protéine,
    fonction biologique, termes Gene Ontology (GO).
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_name} AND organism_name:{organism} AND reviewed:true",
        "format": "json",
        "fields": "accession,id,protein_name,gene_names,cc_function,go_id,go_p,go_f,go_c",
        "size": 1,
    }

    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    data = response.json()

    if not data.get("results"):
        print(f"  [!] Aucun résultat UniProt trouvé pour {gene_name}")
        return None

    entry = data["results"][0]

    # Extraction propre des champs utiles
    accession = entry.get("primaryAccession", "")
    protein_name = (
        entry.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", "")
    )

    function_text = ""
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                function_text = texts[0].get("value", "")
            break

    go_terms = [
        xref.get("properties", [{}])[0].get("value", "")
        for xref in entry.get("uniProtKBCrossReferences", [])
        if xref.get("database") == "GO"
    ]

    return {
        "gene_name": gene_name,
        "uniprot_accession": accession,
        "protein_name": protein_name,
        "function": function_text,
        "go_terms": "; ".join(go_terms[:10]),  # limite à 10 termes pour lisibilité
    }


def collect_all_uniprot_data(gene_list, organism="human"):
    """Boucle sur la liste de gènes et collecte les infos UniProt pour chacun."""
    records = []
    for gene in gene_list:
        print(f"[UniProt] Récupération : {gene}")
        info = get_uniprot_info(gene, organism)
        if info:
            records.append(info)
        time.sleep(0.5)  # pause polie pour ne pas surcharger l'API
    return pd.DataFrame(records)


# ==============================================================================
# ÉTAPE 2 — RÉCUPÉRATION DU RÉSEAU D'INTERACTIONS STRING
# ==============================================================================

def get_string_network(gene_list, species=9606, required_score=700):
    """
    Interroge l'API STRING pour récupérer le réseau d'interactions
    fonctionnelles entre une liste de gènes/protéines.

    required_score : score de confiance minimum (0-1000).
        400 = confiance moyenne (valeur recommandée par défaut)
        700 = confiance élevée
        900 = confiance très élevée
    """
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "network"

    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "required_score": required_score,
        "network_type": "functional",
        "caller_identity": "portfolio_project_cytokines",
    }

    url = f"{string_api_url}/{output_format}/{method}"
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()

    # Le résultat TSV est parsé directement en DataFrame
    from io import StringIO
    df = pd.read_csv(StringIO(response.text), sep="\t")
    return df


def get_string_functional_enrichment(gene_list, species=9606):
    """
    Interroge l'API STRING pour récupérer l'enrichissement fonctionnel
    (GO, KEGG, etc.) de la liste de gènes fournie.
    Utile pour interpréter les communautés détectées plus tard.
    """
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "enrichment"

    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "caller_identity": "portfolio_project_cytokines",
    }

    url = f"{string_api_url}/{output_format}/{method}"
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()

    from io import StringIO
    df = pd.read_csv(StringIO(response.text), sep="\t")
    return df


# ==============================================================================
# EXÉCUTION PRINCIPALE
# ==============================================================================

if __name__ == "__main__":

    print("=" * 70)
    print(f"Collecte de données pour {len(CYTOKINES)} cytokines")
    print("=" * 70)

    # --- 1. UniProt ---
    print("\n[1/3] Récupération des annotations UniProt...")
    uniprot_df = collect_all_uniprot_data(CYTOKINES, organism=ORGANISM)
    uniprot_path = os.path.join(OUTPUT_DIR, "uniprot_annotations.csv")
    uniprot_df.to_csv(uniprot_path, index=False)
    print(f"  -> Sauvegardé dans {uniprot_path} ({len(uniprot_df)} entrées)")

    # --- 2. Réseau STRING ---
    print("\n[2/3] Récupération du réseau d'interactions STRING...")
    network_df = get_string_network(CYTOKINES, species=SPECIES_TAXID, required_score=700)
    network_path = os.path.join(OUTPUT_DIR, "string_network.csv")
    network_df.to_csv(network_path, index=False)
    print(f"  -> Sauvegardé dans {network_path} ({len(network_df)} interactions)")

    # --- 3. Enrichissement fonctionnel ---
    print("\n[3/3] Récupération de l'enrichissement fonctionnel STRING...")
    enrichment_df = get_string_functional_enrichment(CYTOKINES, species=SPECIES_TAXID)
    enrichment_path = os.path.join(OUTPUT_DIR, "string_enrichment.csv")
    enrichment_df.to_csv(enrichment_path, index=False)
    print(f"  -> Sauvegardé dans {enrichment_path} ({len(enrichment_df)} termes)")

    print("\n" + "=" * 70)
    print("Collecte terminée. Fichiers disponibles dans data/raw/ :")
    print("  - uniprot_annotations.csv  (fonctions biologiques par cytokine)")
    print("  - string_network.csv       (interactions entre cytokines)")
    print("  - string_enrichment.csv    (fonctions/voies enrichies)")
    print("=" * 70)