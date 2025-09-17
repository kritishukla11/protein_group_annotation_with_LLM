import requests
import pandas as pd
import io
import time

def map_gene_to_uniprot(gene_symbols, from_db="Gene_Name", to_db="UniProtKB", batch_size=100):
    """
    Map gene symbols to UniProtKB accessions using UniProt ID mapping API.
    Handles async jobs + pagination.
    """
    base_url = "https://rest.uniprot.org/idmapping"
    results = {}

    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i+batch_size]
        ids_str = ",".join(batch)

        # Step 1: submit job
        r = requests.post(
            f"{base_url}/run",
            data={"from": from_db, "to": to_db, "ids": ids_str}
        )
        r.raise_for_status()
        job_id = r.json()["jobId"]

        # Step 2: poll until finished
        while True:
            status = requests.get(f"{base_url}/status/{job_id}")
            status.raise_for_status()
            status_json = status.json()
            if status_json.get("jobStatus") == "FINISHED":
                break
            elif status_json.get("jobStatus") == "FAILED":
                print(f"Batch {i} failed.")
                break
            time.sleep(2)

        # Step 3: fetch results (may be paginated)
        url = f"{base_url}/results/{job_id}"
        while url:
            res = requests.get(url)
            res.raise_for_status()
            res_json = res.json()
            for item in res_json.get("results", []):
                results[item["from"]] = item["to"]["primaryAccession"]
            url = res_json.get("next")  # next page, if any

    return results


def fetch_uniprot_annotations(protein_ids, batch_size=50):
    """
    Fetch GO terms, Pfam, and keywords from UniProt for a list of protein IDs.
    Uses the supported fields for the new REST API.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    results = []

    for i in range(0, len(protein_ids), batch_size):
        batch = protein_ids[i:i+batch_size]
        query = " OR ".join([f"accession:{pid}" for pid in batch])

        params = {
            "query": query,
            "fields": "accession,go_id,xref_pfam,keywords",
            "format": "tsv"
        }

        response = requests.get(base_url, params=params, timeout=120)

        if response.status_code == 200:
            df = pd.read_csv(io.StringIO(response.text), sep="\t")
            results.append(df)
        else:
            print(f"Failed batch {i}: {response.status_code} {response.text[:200]}")

        time.sleep(1)  # avoid rate-limit

    if results:
        return pd.concat(results, ignore_index=True)
    else:
        return pd.DataFrame()

# Load input dataframe
df_neighbors = pd.read_csv('example_input_top10.csv')

# Collect all protein IDs in one list
prot_list = df_neighbors['protein'].tolist()

# Get uniprot IDs
gene_to_uniprot = map_gene_to_uniprot(prot_list)

# Fetch annotations for all proteins
annotations = fetch_uniprot_annotations(uniprot_list)

# Save annotations as csv
pd.DataFrame(annotations).to_csv('example_output_labelmapping.csv',index=False)