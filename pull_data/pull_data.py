import os, pandas as pd, time, tqdm
from mp_api.client import MPRester 
from secret import MP_API_KEY
OUT_DIR = "data_raw"
os.makedirs(OUT_DIR, exist_ok=True)
query = {
    "elements": "O",
    "nelements": 3,
    "formula": "ABX3",
}
with MPRester(MP_API_KEY) as mpr:
    docs=mpr.materials.summary.search(**query,
                            fields=["material_id","band_gap"])
    
records=[]
for d in tqdm.tqdm(docs,desc="Downloading"):
    struct = mpr.get_structure_by_material_id(d.material_id)
    cif_path = os.path.join(OUT_DIR, f"{d.material_id}.cif")
    struct.to(filename=cif_path)
    records.append({
            "mp_id": d.material_id,
            "cif_path": cif_path,
            "bandgap_pbe": d.band_gap
        })

df = pd.DataFrame(records)
df.to_csv("raw_perovskites.csv", index=False)
print(f"✅ 共 {len(df)} 条钙钛矿，CIF 与 CSV 已保存")