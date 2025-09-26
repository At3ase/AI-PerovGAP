import pandas as pd, os, tqdm
from pymatgen.core.structure import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
CSV_PATH = os.path.join(os.path.dirname(__file__), "raw_perovskites.csv")
df = pd.read_csv(CSV_PATH)
matcher = StructureMatcher(ltol=0.1, stol=0.1, angle_tol=5)

unique = []
for _, row in tqdm.tqdm(df.iterrows(), total=len(df)):
    s = Structure.from_file(row.cif_path)
    if not any(matcher.fit(s, Structure.from_file(u.cif_path)) for u in unique):
        unique.append(row)
df_uni = pd.DataFrame(unique).reset_index(drop=True)

def tolerance_factor(row):
    s = Structure.from_file(row.cif_path)
    elem = [s[i].specie for i in range(len(s))]
    rad = {str(sp): sp.ionic_radii["VI"]["+2"] for sp in elem if sp.ionic_radii}
    if len(rad) < 3: return False
    a, b, x = sorted(rad.values(), reverse=True)[:3]
    t = (a + x) / (2**0.5 * (b + x))
    return 0.8 <= t <= 1.05

df_clean = df_uni[df_uni.apply(tolerance_factor, axis=1)].reset_index(drop=True)
df_clean.to_csv("clean_perovskites.csv", index=False)
print(f" 去重+过滤后剩 {len(df_clean)} 条")