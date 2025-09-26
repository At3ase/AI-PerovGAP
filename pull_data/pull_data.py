"""
引用工具包
"""
import os, pandas as pd, time, tqdm
from mp_api.client import MPRester 
from secret import MP_API_KEY
"""
创建独立的输出目录 data_raw"""
OUT_DIR = "data_raw"
os.makedirs(OUT_DIR, exist_ok=True)
"""
使用 MPRester 下载钙钛矿结构,注意替换 MP_API_KEY。
解释一下 query 变量：
"elements": "O" 表示材料中必须包含氧元素。
"nelements": 3 表示材料中总共有3种不同的元素。"
"formula": "ABX3" 表示材料的化学式必须符合钙钛矿的通式，其中 A 和 B 是两种不同的阳离子，X 是阴离子（通常是氧）。
同时提取的材料资源，包含material_id和band_gap字段。
"""
query = {
    "elements": "O",
    "nelements": 3,
    "formula": "ABX3",
}
with MPRester(MP_API_KEY) as mpr:
    docs=mpr.materials.summary.search(**query,
                            fields=["material_id","band_gap"])

"""
下载 CIF 文件并保存到本地
"""    
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
"""
保存为 CSV 文件
"""
df = pd.DataFrame(records)
df.to_csv("raw_perovskites.csv", index=False)
print(f"✅ 共 {len(df)} 条钙钛矿，CIF 与 CSV 已保存")