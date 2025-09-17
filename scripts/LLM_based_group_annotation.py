import pandas as pd
from collections import Counter
import ast
from openai import OpenAI
import os

# Load OpenAI LLM
os.environ["OPENAI_API_KEY"] = "your-API-key-here"
client = OpenAI()

# LLM prompt function
def make_llm_prompt(protein, neighbors, mapping, shared):
    prots = [protein] + neighbors
    m = mapping[mapping["Gene Names"].isin(prots)]

    text_blocks = []
    for _, row in m.iterrows():
        lines = [f"- {row['Gene Names']}"]
        for col in [
            "Gene Ontology (biological process)",
            "Gene Ontology (cellular component)",
            "Gene Ontology (molecular function)",
            "Protein families",
            "Reactome"
        ]:
            if col in m.columns and pd.notna(row.get(col, None)):
                lines.append(f"  {col}: {row[col]}")
        text_blocks.append("\n".join(lines))
    
    return f"""
Group proteins: {', '.join(prots)}

Annotations:
{chr(10).join(text_blocks)}

Task: Based on the annotations, assign a concise, human-readable group name 
(like a Gene Ontology or Pfam family label) that best captures their shared function.
"""

# Load data
top10 = pd.read_csv("../inputs/example_input_top10.csv")
shared = pd.read_csv("../inputs/example_input_commonpathways.csv")
mapping = pd.read_csv("../outputs/example_output_labelmapping.csv")

# Fix neighbors column in top10
top10["closest_proteins"] = top10["closest_proteins"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith("[") else x)

# Fix common_pathways column in shared
shared["common_pathways"] = shared["common_pathways"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith("[") else x)

# Generate prompts for all rows
prompts = []
for _, row in top10.iterrows():
    prompt = make_llm_prompt(
        protein=row["protein"],
        neighbors=row["closest_proteins"],
        mapping=mapping,
        shared=shared
    )
    prompts.append(prompt)

# Run OpenAI on prompt
results = []
for prompt in prompts:
    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role":"user","content":prompt}],
        temperature=0.3
    )
    group_name = resp.choices[0].message.content
    results.append(group_name)

# Save
pd.DataFrame(results).to_csv('example_output_LLMassignedgroup.csv',index=False)
