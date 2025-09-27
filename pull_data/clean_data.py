
import os, pandas as pd, tqdm
from pymatgen.core.structure import Structure

CSV_PATH = os.path.join(os.path.dirname(__file__), "raw_perovskites.csv")
df = pd.read_csv(CSV_PATH)

print("一次性加载结构...")
structs = [Structure.from_file(p) for p in tqdm.tqdm(df["cif_path"], desc="Load")]

seen_comp = set()
keep_idx  = []
for i, s in enumerate(tqdm.tqdm(structs, desc="Dedup")):
    comp = tuple(sorted(s.composition.as_dict().items()))   # 仅成分
    if comp not in seen_comp:
        seen_comp.add(comp)
        keep_idx.append(i)

df_uni = df.iloc[keep_idx].reset_index(drop=True)
print(f"成分去重：{len(df)} → {len(df_uni)}")
df_uni.to_csv("clean_perovskites.csv", index=False)
print("✅ 成分去重完成，未过滤，仅去重")