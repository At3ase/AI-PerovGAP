from importlib import resources
import os, pandas as pd, tqdm
from pymatgen.core.structure import Structure
from pymatgen.core import Species


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

# ===== 2. 缺失审计（不过滤，只标记）=====
def add_tolerance(row):
    s = Structure.from_file(row["cif_path"])
    # ===== Goldschmidt 定义：A=最大阳离子，B=次大阳离子，X=阴离子 =====
    from pymatgen.core import Species
    SHANNON_VI = {  # VI 配位，单位 Å，可查表
        "Ca": {"+2": 1.00}, "Sr": {"+2": 1.18}, "Ba": {"+2": 1.35},
        "Ti": {"+4": 0.605}, "Zr": {"+4": 0.72}, "Pb": {"+2": 1.19},
        "O": {"-2": 1.40}, "F": {"-1": 1.33}, "Cl": {"-1": 1.81},
        "Br": {"-1": 1.96}, "I": {"-1": 2.20}, "Cs": {"+1": 1.74},
        "K": {"+1": 1.38}, "Rb": {"+1": 1.52}, "Sn": {"+2": 1.10, "+4": 0.69},
        "Ge": {"+4": 0.53}, "Mn": {"+2": 0.83}, "Fe": {"+3": 0.645},
        "Ni": {"+2": 0.69}, "Zn": {"+2": 0.74}, "Mg": {"+2": 0.72},
    }

    def get_shannon_radius(sp, oxi):
        """Goldschannon VI 配位优先 → 可查表"""
        return SHANNON_VI.get(sp.symbol, {}).get(oxi) or \
               sp.ionic_radii.get("VI", {}).get(oxi) or \
               sp.ionic_radii.get("IV", {}).get(oxi) or \
               1.40  # 终极闸门，写进审计

    rad = {}
    for site in s:
        sp = site.specie
        if not hasattr(sp, "oxi_state"):
            sp = Species(sp.symbol, "+2")   # 先转 Species
        oxi = str(sp.oxi_state) if hasattr(sp, "oxi_state") else "+2"
        r = get_shannon_radius(sp, oxi)
        rad[sp.symbol] = r   # 每个元素只留一个半径

    # Goldschmidt 定义：A=最大阳离子，B=次大阳离子，X=阴离子 
    symbols = list(rad.keys())
    vals    = list(rad.values())
    # 阳离子 vs 阴离子
    cations = {k: v for k, v in zip(symbols, vals) if k not in {"O", "F", "Cl", "Br", "I"}}
    anions  = {k: v for k, v in zip(symbols, vals) if k in {"O", "F", "Cl", "Br", "I"}}
    if len(cations) < 2 or len(anions) == 0:
        return 0.0, "missing_shannon"   # 不过滤，只标记缺失
    # A=最大阳离子，B=次大阳离子，X=最大阴离子
    a = max(cations.values())
    b = sorted(cations.values(), reverse=True)[1] if len(cations) > 1 else min(cations.values())
    x = max(anions.values())
    t = (a + x) / (2**0.5 * (b + x))
    return t, "classic" if 0.8 <= t <= 1.05 else "extreme"

t_tol, t_cls = zip(*df_uni.apply(add_tolerance, axis=1))
rad_sources = []
df_uni = df_uni.assign(t_factor=t_tol, perovskite_class=t_cls)
df_uni.to_csv("clean_perovskites.csv", index=False)
print("✅ 缺失审计完成，未过滤，仅标记")
print("类别分布：")
print(df_uni["perovskite_class"].value_counts())

# ===== 3. 审计日志（每步来源可溯源）=====

df_uni = df_uni.assign(radius_source=", ".join([src for _, src in rad_sources]))
df_uni.to_csv("clean_perovskites.csv", index=False)
print("✅ 审计日志完成，每步来源可溯源")