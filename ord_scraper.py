import sys
import os
import json
import tempfile
from urllib.parse import urlparse
import requests
from tqdm import tqdm

ROLE_MAP = {0: "UNSPECIFIED", 1: "REACTANT", 2: "REAGENT", 3: "SOLVENT", 4: "CATALYST", 5: "LIGAND"}

def extract_dataset_id(arg):
    text = arg.strip()
    if text.startswith("http"):
        parts = [p for p in urlparse(text).path.split("/") if p]
        return parts[-1] if parts else None
    if text.startswith("ord_dataset-"):
        return text
    return None

def load_ord_dataset_from_file(path):
    from ord_schema.message_helpers import load_message
    from ord_schema.proto import dataset_pb2
    return load_message(path, dataset_pb2.Dataset)

def load_ord_dataset_by_id(dataset_id):
    from ord_schema.message_helpers import fetch_dataset
    return fetch_dataset(dataset_id)

def get_identifier_entry(compound):
    from ord_schema import message_helpers
    try:
        smiles = message_helpers.smiles_from_compound(compound, canonical=False) or None
    except Exception:
        smiles = None
    if smiles:
        return smiles, "SMILES"
    first = next(iter(compound.identifiers), None)
    if first is not None:
        return first.value, "UNKNOWN"
    return "", ""

def role_value(compound):
    try:
        return int(getattr(compound, "reaction_role", 0))
    except Exception:
        return 0

def is_metal_smiles(smiles, name):
    s = smiles or ""
    n = (name or "").lower()
    if any(tag in s for tag in ["[Pd]", "[Ni]", "[Cu]", "[Fe]", "[Co]", "[Ru]", "[Rh]", "[Ir]", "[Pt]"]):
        return True
    if any(tag in n for tag in ["palladium", "nickel", "copper", "iron", "cobalt", "ruthenium", "rhodium", "iridium", "platinum", "pd", "ni", "cu", "fe", "co", "ru", "rh", "ir", "pt"]):
        return True
    return False

def extract_reaction_fields(dataset):
    out = []
    for rxn in dataset.reactions:
        rec = {
            "dataset_id": getattr(dataset, "dataset_id", ""),
            "reaction_id": getattr(rxn, "reaction_id", ""),
            "base": [],
            "solvent": [],
            "amine": [],
            "aryl_halide": [],
            "metal": [],
            "ligand": [],
            "carboxylic_acid": [],
            "additive": [],
            "activation_agent": [],
            "m1": [],
            "m2": [],
            "m3": [],
            "m4": [],
            "m5": [],
            "m6": [],
            "m7": [],
            "m1_m4": [],
            "m1_m3": [],
            "m2_m3": [],
            "m1_m2_m3_m6": [],
            "m1_m6_m2": [],
            "m1_m2_m7": [],
            "m2_m4": []
        }
        for name, inp in rxn.inputs.items():
            for comp in inp.components:
                val, id_type = get_identifier_entry(comp)
                rval = role_value(comp)
                entry = {"value": val, "reaction_role_name": ROLE_MAP.get(rval, ""), "input_name": name}
                lname = (name or "").lower()
                if lname == "solvent":
                    rec["solvent"].append(entry)
                elif lname == "base":
                    rec["base"].append(entry)
                elif lname == "amine":
                    rec["amine"].append(entry)
                elif lname == "aryl halide":
                    rec["aryl_halide"].append(entry)
                elif lname == "carboxylic acid":
                    rec["carboxylic_acid"].append(entry)
                elif lname == "additive":
                    rec["additive"].append(entry)
                elif lname == "activation agent":
                    rec["activation_agent"].append(entry)
                elif lname in ("m1", "m2", "m3", "m4", "m5", "m6", "m7"):
                    rec.setdefault(lname, []).append(entry)
                elif lname in ("m1_m4", "m1_m3", "m2_m3", "m1_m2_m3_m6", "m1_m6_m2", "m1_m2_m7", "m2_m4"):
                    rec.setdefault(lname, []).append(entry)
                elif lname == "metal and ligand":
                    smi = val if id_type == "SMILES" else ""
                    name_lower = (val if id_type != "SMILES" else "").lower()
                    if is_metal_smiles(smi, name_lower):
                        rec["metal"].append(entry)
                    else:
                        rec["ligand"].append(entry)
        rec = {k: v for k, v in rec.items() if not (isinstance(v, list) and len(v) == 0)}
        out.append(rec)
    return out

def download_file(session, url):
    with session.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        fd, path = tempfile.mkstemp(suffix=os.path.splitext(urlparse(url).path)[1])
        with os.fdopen(fd, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    return path

def list_ord_datasets_via_github(limit):
    base = "https://api.github.com/repos/Open-Reaction-Database/ord-data/contents/data"
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "Armens-ORD-Scraper/1.0"
    }
    session = requests.Session()
    session.headers.update(headers)

    def walk(path, acc):
        if len(acc) >= limit:
            return
        r = session.get(f"{base}{path}", timeout=30)
        r.raise_for_status()
        items = r.json()
        for it in items:
            if len(acc) >= limit:
                break
            if it.get("type") == "file" and (it["name"].endswith(".pb.gz") or it["name"].endswith(".pb")):
                raw_url = f"https://github.com/Open-Reaction-Database/ord-data/raw/main/{it['path']}"
                acc.append(raw_url)
            elif it.get("type") == "dir":
                walk(f"/{it['name']}", acc)
        return

    urls = []
    walk("", urls)
    return urls

def find_reaction_by_id(reaction_id):
    api = "https://api.github.com/search/code"
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "Armens-ORD-Scraper/1.0"
    }
    r = requests.get(api, params={"q": f"{reaction_id} repo:Open-Reaction-Database/ord-data"}, headers=headers, timeout=30)
    if r.status_code != 200:
        return None
    items = r.json().get("items", [])
    for it in items:
        raw_url = f"https://github.com/Open-Reaction-Database/ord-data/raw/main/{it['path']}"
        try:
            sess = requests.Session()
            tmp_path = download_file(sess, raw_url)
            ds = load_ord_dataset_from_file(tmp_path)
            for rxn in ds.reactions:
                if getattr(rxn, "reaction_id", "") == reaction_id:
                    return rxn
        except Exception:
            continue
    return None

def scrape(limit=5, dataset_id=None, dataset_ids=None):
    session = requests.Session()
    out_path = os.path.join(os.getcwd(), "ord_browse_extract.json")
    if dataset_ids:
        results = []
        for did in dataset_ids:
            try:
                ds = load_ord_dataset_by_id(did)
                if len(ds.reactions) == 0 and len(ds.reaction_ids) > 0:
                    for rid in ds.reaction_ids:
                        rxn = find_reaction_by_id(rid)
                        if rxn is None:
                            continue
                        class DWrap:
                            dataset_id = did
                            reactions = [rxn]
                        results.extend(extract_reaction_fields(DWrap))
                else:
                    results.extend(extract_reaction_fields(ds))
            except Exception as e:
                print("error", str(e))
                continue
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(results, f, ensure_ascii=False, indent=2)
        return out_path, len(results)
    if dataset_id:
        try:
            ds = load_ord_dataset_by_id(dataset_id)
            if len(ds.reactions) == 0 and len(ds.reaction_ids) > 0:
                results = []
                for rid in ds.reaction_ids:
                    rxn = find_reaction_by_id(rid)
                    if rxn is None:
                        continue
                    class DWrap:
                        dataset_id = dataset_id
                        reactions = [rxn]
                    results.extend(extract_reaction_fields(DWrap))
            else:
                results = extract_reaction_fields(ds)
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(results, f, ensure_ascii=False, indent=2)
            return out_path, len(results)
        except Exception as e:
            print("error", str(e))
    urls = list_ord_datasets_via_github(limit)
    results = []
    for file_url in tqdm(urls):
        try:
            local_path = download_file(session, file_url)
            ds = load_ord_dataset_from_file(local_path)
            results.extend(extract_reaction_fields(ds))
        except Exception as e:
            print("error", str(e))
            continue
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    return out_path, len(results)

def main():
    args = sys.argv[1:]
    dataset_ids = []
    limit = 5
    for arg in args:
        dsid = extract_dataset_id(arg)
        if dsid:
            dataset_ids.append(dsid)
        else:
            try:
                limit = int(arg)
            except Exception:
                pass
    path, n = scrape(limit=limit, dataset_ids=dataset_ids if dataset_ids else None)
    print(path)
    print(n)

if __name__ == "__main__":
    main()
