#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
P7b_fix_order.py — 修复图表编号顺序与表格重编号
问题：
  1) 图 6 块在图 5 块之前出现（文档序）
  2) 图 8 块在图 9 块之后出现（文档序）
  3) 表 4（多因素）定义在表 3（外部汇总）之前 → 重编号：多因素=表 3，外部汇总=表 4
  4) 外部汇总表注与后续正文粘连（缺空行）
"""
import re, sys

MD = r"D:\projects\OS_matrix_dualKO\manuscript\OS_manuscript_v1_cn.md"
t = open(MD, encoding="utf-8").read()
orig = len(t)

def must(old, new, tag, count=1):
    global t
    n = t.count(old)
    if n != count:
        sys.exit(f"[FAIL] {tag}: 命中 {n} 次（应为 {count}）\n{old[:150]}")
    t = t.replace(old, new)

def move_block(fig_no, img_file, en_caption, new_anchor):
    """把图 N 块（从图片行到英文图注行）摘出并插入 new_anchor 之后（规范化空行）。"""
    global t
    pat = re.compile(
        r"\n!\[图 %d[｜|][^\]]*\]\(\.\./results/figures/%s\)\n\n.*?\*%s\*\n"
        % (fig_no, img_file, re.escape(en_caption)), re.S)
    m = pat.search(t)
    if not m:
        sys.exit(f"[FAIL] move_block 图 {fig_no}: 未找到块")
    block = m.group(0)
    t = t[:m.start()] + t[m.end():]
    must(new_anchor, new_anchor + "\n\n" + block.strip("\n") + "\n",
         f"reinsert 图 {fig_no}")
    print(f"[move] 图 {fig_no} 已移动")

# ---------------------------------------------------------------- 1) 粘连修复
must("而非过拟合。*在共享训练终点",
     "而非过拟合。*\n\n在共享训练终点",
     "粘连修复")

# ---------------------------------------------------------------- 2) 表格重编号（多因素→表3，外部汇总→表4）
# 2a. 多因素表定义改名 + 恢复"（表 3）"引用（两处相邻，一次替换）
must("，OS终点的一致性指数为0.786。\n\n**表 4｜风险评分的多因素独立性分析**",
     "，OS终点的一致性指数为0.786（表 3）。\n\n**表 3｜风险评分的多因素独立性分析**",
     "多因素表改名+引用")
# 2b. 外部汇总表定义改名
must("**表 3｜14 基因签名在训练与四个独立队列中的表现汇总**",
     "**表 4｜14 基因签名在训练与四个独立队列中的表现汇总**",
     "外部汇总表改名")
# 2c. 外部汇总引用改名
must("四个独立队列的外部验证总结于表 3。",
     "四个独立队列的外部验证总结于表 4。",
     "外部汇总引用")
# 2d. GSE39055 多因素段落引用指向表 3
must("（图 8；表 4）。",
     "（图 8；表 3）。",
     "GSE39055多因素引用")
# 2e. 图 8 图注内的表引用
must("多因素独立性结果见表 4。",
     "多因素独立性结果见表 3。",
     "图8图注表引用")
# 2f. 图表清单两行对调
must("| 表 3 | 训练与四个独立队列的表现汇总 |",
     "| 表 99 | 训练与四个独立队列的表现汇总 |", "清单-临时")
must("| 表 4 | 风险评分的多因素独立性分析 |",
     "| 表 3 | 风险评分的多因素独立性分析 |", "清单-多因素")
must("| 表 99 | 训练与四个独立队列的表现汇总 |",
     "| 表 4 | 训练与四个独立队列的表现汇总 |", "清单-外部汇总")

# ---------------------------------------------------------------- 3) 图块移动
# 图 8: 从 GSE21257 敏感性段（Para C）末尾 → GSE39055 多因素段（Para A）末尾
move_block(8, "F8.png", "Figure 8. External validation in the GSE39055 cohort.",
           "同时坏死也显著（HR = 0.23，p = 0.0228）（图 8；表 3）。")
# 图 5: 从正交性段末尾 → 双 KO 首段（设计描述）末尾
move_block(5, "F5.png", "Figure 5. Perturbation landscape of the parallel dual-gene virtual knockout.",
           "每个区室采样1,000个细胞，降维后保留1,501个基因（图 5A-C）。")

# ---------------------------------------------------------------- 4) 空行规范化
t = re.sub(r"\n{3,}", "\n\n", t)

# ---------------------------------------------------------------- 5) 校验
def order_of(pattern):
    return [int(m) for m in re.findall(pattern, t)]

fig_def_order = order_of(r"!\[图 (\d+)｜")
tbl_def_order = order_of(r"\*\*表 (\d)｜")
print("图定义顺序:", fig_def_order)
print("表定义顺序:", tbl_def_order)
assert fig_def_order == sorted(fig_def_order), "图定义顺序不对"
assert tbl_def_order == sorted(tbl_def_order), "表定义顺序不对"

# 首次引用顺序
def first_cite(pat):
    return [int(m.group(1)) for m in re.finditer(pat, t)]
fig_cites = first_cite(r"（图 (\d+)[A-Z\-C]")
tbl_cites = first_cite(r"(?<!定)（?表 (\d)")
print("图首次引用(截取):", sorted(set(fig_cites)))
print("表引用出现序:", tbl_cites[:12])

# 关键数值完好性
KEY = ['0.796', '0.696', '48.4%', '1.02e-11', '5.43', '2.69', '0.0342', '0.0021',
       '0.0017', '3.49', '0.308', '0.0084', '0.0255', '84.2%', '49.2%', '73.6%',
       '11.23', '0.385', '2.15e-4', '9.1e-7', '8/17', '4/20', '5.57', '2.92', '0.786']
miss = [k for k in KEY if k not in t]
print("关键数值缺失:", miss if miss else "无")
assert not miss

open(MD, "w", encoding="utf-8").write(t)
print(f"[write] {MD}  ({orig} -> {len(t)} chars)")
