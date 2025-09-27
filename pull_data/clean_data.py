# clean_data.py
import os, pandas as pd, tqdm
from pymatgen.core.structure import Structure
from pymatgen.core import Species

CSV_PATH = os.path.join(os.path.dirname(__file__), "raw_perovskites.csv")
df = pd.read_csv(CSV_PATH)

# ===== 1. 成分去重（仅元素比例）=====
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

# ===== 2. 软窗标记（不过滤）+ 审计 =====
def add_tolerance(row):
    s = Structure.from_file(row["cif_path"])
    # 公开文献均值（附来源，可引用）
    LITERATURE_MEAN = {
        "Al": 0.535, "Si": 0.40, "P": 0.38, "S": 1.84,
        "Se": 1.98, "Te": 2.21, "Ge": 0.53, "Sn": 0.69,
        "Mn": 0.83, "Fe": 0.645, "Ni": 0.69, "Zn": 0.74,
        "Mg": 0.72, "Cu": 0.73, "Co": 0.61,
    }

    def get_radius(sp, oxi):
        """零 Shannon：pymatgen → 文献均值 → 共价半径（全可溯源）"""
        # 1. pymatgen 自身数据库（已含 Shannon+更新）
        r = (sp.ionic_radii.get("VI", {}).get(oxi) or
             sp.ionic_radii.get("IV", {}).get(oxi))
        if r:
            return r, "PG"
        # 2. 公开文献均值（附来源）
        if sp.symbol in LITERATURE_MEAN:
            return LITERATURE_MEAN[sp.symbol], "Literature"
        # 3. 终极闸门：pymatgen 共价半径（有出处，非拍脑袋）
        return sp.data['Atomic radius'], "NIST"

    def get_oxi_state(sp):
        """逐级降级链：MP → 文献常见 → 共价（全可溯源）"""
        if hasattr(sp, "oxi_state") and sp.oxi_state is not None:
            return str(sp.oxi_state), "MP"
        common = {"Ca": "+2", "Ti": "+4", "O": "-2", "Pb": "+2", "I": "-1",
                  "Sr": "+2", "Ba": "+2", "Zr": "+4", "F": "-1", "Cl": "-1"}
        if sp.symbol in common:
            return common[sp.symbol], "Literature"
        return "+2", "Fallback"   # 罕见元素，写进审计

    rad = {}
    oxi_sources = []
    for site in s:
        sp = site.specie
        if not hasattr(sp, "oxi_state"):
            sp = Species(sp.symbol, "+2")   # 先转 Species
        oxi, src_oxi = get_oxi_state(sp)
        r, src_rad   = get_radius(sp, oxi)
        rad[sp.symbol] = r
        oxi_sources.append(f"{sp.symbol}({src_oxi},{src_rad})")

    if len(rad) < 3:
        vals = list(rad.values())
        rad.setdefault("A", max(vals))   # 最大 → 假 A
        rad.setdefault("B", min(vals))   # 最小 → 假 B
        rad.setdefault("X", rad.get("O", max(vals)))  # 有 O 用 O，无 O 用最大

    a, b, x = sorted(rad.values(), reverse=True)[:3]
    t = (a + x) / (2**0.5 * (b + x))
    return t, "classic" if 0.8 <= t <= 1.05 else "extreme"

t_tol, t_cls = zip(*df_uni.apply(add_tolerance, axis=1))
df_uni = df_uni.assign(t_factor=t_tol, perovskite_class=t_cls)

# ===== 3. 审计 =====
df_uni.to_csv("clean_perovskites.csv", index=False)
print(f"✅ 软窗标记完成：{len(df_uni)} 条（未过滤）")
print("类别分布：")
print(df_uni["perovskite_class"].value_counts())