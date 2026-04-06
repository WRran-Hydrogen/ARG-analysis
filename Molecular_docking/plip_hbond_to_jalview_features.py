#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, glob
from pathlib import Path
import pandas as pd
import xml.etree.ElementTree as ET

# PLIP XML 中常见的相互作用条目标签（element tag）
# 你的示例XML里：hydrophobic_interaction / hydrogen_bond / salt_bridge 都存在且字段齐全。[1](https://ueanorwich-my.sharepoint.com/personal/hcb24wcu_uea_ac_uk/Documents/Microsoft%20Copilot%20Chat%20Files/1bls_P99_H_pose_report.txt)
PLIP_TAG_TO_FEATURETYPE = {
    "hydrogen_bond": "HBOND",
    "hydrophobic_interaction": "HYDROPHOBIC",
    "salt_bridge": "SALT_BRIDGE",
    "water_bridge": "WATER_BRIDGE",
    "pi_stack": "PI_STACK",
    "pi_cation_interaction": "PI_CATION",
    "halogen_bond": "HALOGEN_BOND",
    "metal_complex": "METAL_COMPLEX",
}

# 给 Jalview features 文件写颜色定义（FeatureType -> color）
# Jalview 允许每个 FeatureType 在文件开头用 “type<TAB>color” 定义显示风格。[2](https://vicfero.github.io/trimal/)
DEFAULT_COLORS = {
    "HBOND": "128,128,128",           # 灰
    "HYDROPHOBIC": "255,200,0",       # 黄
    "SALT_BRIDGE": "255,0,0",         # 红
    "WATER_BRIDGE": "0,180,255",      # 蓝
    "PI_STACK": "160,0,255",          # 紫
    "PI_CATION": "255,0,160",         # 粉
    "HALOGEN_BOND": "0,200,0",        # 绿
    "METAL_COMPLEX": "120,80,40",     # 棕
}

def iter_files(items):
    files=[]
    for it in items:
        it=os.path.expanduser(it)
        m=glob.glob(it)
        if m: files.extend(m)
        elif os.path.exists(it): files.append(it)
    seen=set(); out=[]
    for f in files:
        if f not in seen:
            out.append(f); seen.add(f)
    return out

def load_seqmap(path: Path):
    """
    seqmap.tsv: seq_pos, chain, resnr, icode, resname3 ...
    返回: (chain,resnr) -> seq_pos（重复则保留第一个）
    """
    df=pd.read_csv(path, sep="\t")
    mp={}
    for _,r in df.iterrows():
        key=(str(r["chain"]).strip(), int(r["resnr"]))
        mp.setdefault(key, int(r["seq_pos"]))
    return mp

def parse_plip_xml_all(xml_path: Path, enabled_featuretypes: set):
    """
    解析 PLIP XML，提取所有相互作用残基信息：
      返回: list of dicts [{feature_type, chain, resnr, restype, extra_desc}, ...]
    PLIP 每条相互作用记录一般含 <resnr>, <reschain>, <restype> 等字段（你的XML氢键/疏水/盐桥均如此）。[1](https://ueanorwich-my.sharepoint.com/personal/hcb24wcu_uea_ac_uk/Documents/Microsoft%20Copilot%20Chat%20Files/1bls_P99_H_pose_report.txt)
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    hits = []
    for elem in root.iter():
        tag = elem.tag.lower().strip()
        if tag not in PLIP_TAG_TO_FEATURETYPE:
            continue
        ftype = PLIP_TAG_TO_FEATURETYPE[tag]
        if enabled_featuretypes and ftype not in enabled_featuretypes:
            continue

        # 提取子节点字段
        sub = {c.tag.lower().strip(): (c.text.strip() if c.text else "") for c in list(elem)}
        chain = sub.get("reschain")
        resnr = sub.get("resnr")
        restype = sub.get("restype")  # 可用于desc

        if not chain or not resnr:
            continue
        try:
            resnr_i = int(resnr)
        except:
            continue

        # 额外信息放入description（可选）
        # 比如 hbond 有 dist_h-a / don_angle 等；疏水/盐桥也有 dist 等。[1](https://ueanorwich-my.sharepoint.com/personal/hcb24wcu_uea_ac_uk/Documents/Microsoft%20Copilot%20Chat%20Files/1bls_P99_H_pose_report.txt)
        extra = []
        for k in ("dist", "dist_h-a", "dist_d-a", "don_angle", "protisdon"):
            if k in sub and sub[k]:
                extra.append(f"{k}={sub[k]}")
        extra_desc = ";".join(extra) if extra else ""

        hits.append({
            "feature_type": ftype,
            "chain": chain.strip(),
            "resnr": resnr_i,
            "restype": (restype.strip() if restype else ""),
            "extra_desc": extra_desc,
        })

    return hits

def choose_enabled_types(types_arg: str):
    """
    --types 支持：
      all
      HBOND,HYDROPHOBIC,SALT_BRIDGE,...
      或 tag名：hbond,hydrophobic,salt_bridge...
    """
    if not types_arg or types_arg.lower() == "all":
        return set(DEFAULT_COLORS.keys())  # 全部
    parts = [x.strip() for x in types_arg.split(",") if x.strip()]
    enabled=set()
    # 允许用户写小写简称
    alias = {
        "hbond":"HBOND",
        "hydrogen":"HBOND",
        "hydrophobic":"HYDROPHOBIC",
        "salt":"SALT_BRIDGE",
        "salt_bridge":"SALT_BRIDGE",
        "water":"WATER_BRIDGE",
        "water_bridge":"WATER_BRIDGE",
        "pi_stack":"PI_STACK",
        "pistack":"PI_STACK",
        "pi_cation":"PI_CATION",
        "pi_cation_interaction":"PI_CATION",
        "halogen":"HALOGEN_BOND",
        "halogen_bond":"HALOGEN_BOND",
        "metal":"METAL_COMPLEX",
        "metal_complex":"METAL_COMPLEX",
    }
    for p in parts:
        key = alias.get(p.lower(), p.upper())
        enabled.add(key)
    return enabled

def main():
    ap = argparse.ArgumentParser(
        description="从 PLIP XML 提取所有(或指定)非共价相互作用残基，并映射到序列坐标(seq_pos)，输出 Jalview features 文件。"
    )
    ap.add_argument("--seqmaps", nargs="+", required=True, help="seqmap.tsv 列表（来自 pdb_to_fasta_seqmap.py）")
    ap.add_argument("--plip", nargs="+", required=True, help="PLIP XML 报告列表（可通配符），建议用 *_report.xml")
    ap.add_argument("-o","--out", default="plip_noncovalent.features", help="输出 Jalview features 文件名")
    ap.add_argument("--match_by_prefix", action="store_true",
                    help="按文件名前缀自动配对：seqmap stem 与 plip stem 互相包含则匹配")
    ap.add_argument("--types", default="all",
                    help="要输出的相互作用类型：all 或逗号分隔，如 HBOND,HYDROPHOBIC,SALT_BRIDGE")
    ap.add_argument("--merge_same_residue", action="store_true",
                    help="同一残基若有多种相互作用，合并成一个 NONCOV feature（description里列出类型）")
    ap.add_argument("--merged_type_name", default="NONCOV", help="合并模式下的FeatureType名称")
    ap.add_argument("--merged_color", default="0,0,0", help="合并模式下的颜色（默认黑）")
    args = ap.parse_args()

    enabled = choose_enabled_types(args.types)

    seqmap_files=[Path(x) for x in iter_files(args.seqmaps)]
    plip_files=[Path(x) for x in iter_files(args.plip)]

    seqmap_by={f.stem: f for f in seqmap_files}
    plip_by={f.stem: f for f in plip_files}

    pairs=[]
    if args.match_by_prefix:
        for sm_stem, sm in seqmap_by.items():
            pp = next((f for pstem,f in plip_by.items() if (sm_stem in pstem) or (pstem in sm_stem)), None)
            if pp:
                pairs.append((sm, pp))
    else:
        if len(seqmap_files) != len(plip_files):
            raise SystemExit("未启用 --match_by_prefix 时，--seqmaps 与 --plip 数量必须相同以便一一对应。")
        pairs=list(zip(seqmap_files, plip_files))

    if not pairs:
        raise SystemExit("没有匹配到任何 seqmap ↔ plip 对，请检查文件命名或使用 --match_by_prefix。")

    out_path=Path(args.out)
    lines=[]

    # 1) 写颜色定义块（Jalview features 文件允许先定义每个FeatureType的显示风格）。[2](https://vicfero.github.io/trimal/)
    if args.merge_same_residue:
        lines.append(f"{args.merged_type_name}\t{args.merged_color}")
    else:
        for ft in sorted(enabled):
            color = DEFAULT_COLORS.get(ft, "128,128,128")
            lines.append(f"{ft}\t{color}")

    # 2) 写 feature 行
    for sm_path, plip_path in pairs:
        # seq_id 必须与对齐FASTA header（第一个token）一致，否则 Jalview 不会显示。[2](https://vicfero.github.io/trimal/)
        seq_id = sm_path.stem.replace(".seqmap","")
        mp = load_seqmap(sm_path)

        if plip_path.suffix.lower() != ".xml":
            raise SystemExit(f"当前脚本仅推荐解析 XML：{plip_path}（TXT全类型解析不稳定）")

        hits = parse_plip_xml_all(plip_path, enabled_featuretypes=(set() if args.merge_same_residue else enabled))

        if args.merge_same_residue:
            # (chain,resnr)-> set(types)
            bucket={}
            extra_info={}
            for h in hits:
                key=(h["chain"], h["resnr"])
                bucket.setdefault(key, set()).add(h["feature_type"])
                if h["extra_desc"]:
                    extra_info.setdefault(key, set()).add(h["extra_desc"])
            for (chain,resnr), typeset in sorted(bucket.items()):
                if (chain,resnr) not in mp:
                    continue
                pos = mp[(chain,resnr)]
                desc = f"PLIP_NONCOV;chain={chain};resnr={resnr};types={','.join(sorted(typeset))}"
                if extra_info.get((chain,resnr)):
                    desc += f";info={'|'.join(sorted(extra_info[(chain,resnr)]))}"
                # Jalview feature行：description sequenceId sequenceIndex start end featureType [score] [2](https://vicfero.github.io/trimal/)
                lines.append(f"{desc}\t{seq_id}\t-1\t{pos}\t{pos}\t{args.merged_type_name}")
        else:
            # 多类型分别输出
            for h in hits:
                chain,resnr,ftype = h["chain"], h["resnr"], h["feature_type"]
                if ftype not in enabled:
                    continue
                if (chain,resnr) not in mp:
                    continue
                pos = mp[(chain,resnr)]
                restype = h["restype"]
                extra = h["extra_desc"]
                desc = f"PLIP_{ftype};chain={chain};resnr={resnr}"
                if restype:
                    desc += f";aa={restype}"
                if extra:
                    desc += f";{extra}"
                lines.append(f"{desc}\t{seq_id}\t-1\t{pos}\t{pos}\t{ftype}")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[OK] wrote Jalview features -> {out_path}")

if __name__=="__main__":
    main()
