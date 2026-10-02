# -*- coding: utf-8 -*-
"""
按方案2改写手稿标题与卖点段落。
原则：只改题名、摘要首句、引言最后一段（创新点陈述）与结尾；正文所有数据一字不动。
"""
import io, os, re, shutil

SRC = r"D:/projects/OS_matrix_dualKO/manuscript/OS_manuscript_v1_cn.md"
BAK = r"D:/projects/OS_matrix_dualKO/manuscript/OS_manuscript_v1_cn_pre_title.md"
shutil.copyfile(SRC, BAK)
t = io.open(SRC, encoding="utf-8").read()
orig_len = len(t)

OLD_TITLE = "# 基于单细胞与多组学整合的骨肉瘤预后分层：BUB1与RUNX2双枢纽的机制解析及14基因LASSO-Cox模型构建"
NEW_TITLE = "# 骨肉瘤基质重塑与增殖程序的区室解析：双基因虚拟敲除与14基因事件无进展生存评分的整合研究"

assert OLD_TITLE in t, "标题未找到"
t = t.replace(OLD_TITLE, NEW_TITLE, 1)

# ---- 摘要首句：把"整合多组学数据"换成区室解析 + 双基因虚拟敲除导向 ----
OLD_ABS = "本研究旨在整合多组学数据，鉴定骨肉瘤肿瘤上调且与不良结局相关的候选基因池及两个分属不同功能区室的枢纽基因，解析其细胞来源、恶性属性与下游调控程序，并构建可外部验证的预后风险评分。"
NEW_ABS = ("骨肉瘤的基质重塑与增殖程序在肿瘤生态系统中分工不清，且缺乏在治疗前稳定识别高危患者的分子工具。"
           "本研究旨在以区室解析为框架，鉴定骨肉瘤中分属不同功能区室的两个关键基因，"
           "通过双基因虚拟敲除解析其下游转录程序的分工与独立性，并构建可外部验证的事件无进展生存评分。")
assert OLD_ABS in t, "摘要首句未找到"
t = t.replace(OLD_ABS, NEW_ABS, 1)

# ---- 引言：在首段后补一句"区室分工不清"的缺口，避免红海开头 ----
OLD_INTRO = "高通量转录组研究已在骨肉瘤中报道了若干预后相关基因集和风险模型，但多数研究受限于单队列训练、缺乏外部独立验证或未系统校正平台与终点差异"
NEW_INTRO = ("更重要的是，骨肉瘤中驱动基质重塑与驱动增殖的程序分别由哪些细胞区室承担、二者在下游转录调控层面是否相互独立，"
             "迄今未被系统解析。高通量转录组研究虽已报道了若干预后相关基因集和风险模型，但多数研究受限于单队列训练、"
             "缺乏外部独立验证或未系统校正平台与终点差异")
assert OLD_INTRO in t, "引言段未找到"
t = t.replace(OLD_INTRO, NEW_INTRO, 1)

# ---- 引言末段（创新点陈述）：把"整合...构建递进式分析框架"提到区室对比架构 ----
OLD_INNO = "基于此，本研究整合TARGET-OS多队列bulk转录组、五个公共微阵列/RNA-seq队列、单细胞转录组（GSE162454）、体细胞突变与ESTIMATE微环境评分，构建“候选池筛选—枢纽基因鉴定—单细胞定位—通讯网络—虚拟敲除—预后建模—外部验证”的递进式分析框架。"
NEW_INNO = ("基于此，本研究整合TARGET-OS多队列bulk转录组、五个公共微阵列/RNA-seq队列、单细胞转录组（GSE162454）、"
            "体细胞突变与ESTIMATE微环境评分，构建以“区室解析”为主线的递进式框架："
            "候选池筛选—关键基因鉴定—单细胞区室定位—通讯网络—双基因虚拟敲除—预后建模—外部验证。")
assert OLD_INNO in t, "引言创新点段未找到"
t = t.replace(OLD_INNO, NEW_INNO, 1)

io.open(SRC, "w", encoding="utf-8").write(t)
print("OK  title rewritten")
print("len:", orig_len, "->", len(t), "(delta %+d)" % (len(t) - orig_len))
print("backup:", os.path.basename(BAK))
print("NEW TITLE:", NEW_TITLE)
