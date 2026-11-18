# UESTC-001 · 平面有向反馈顶点集核化改进

- 学校 / 研究者：电子科技大学 · 肖鸣宇（Mingyu Xiao）
- 领域：参数化算法
- 证据类型：author-stated
- 当前状态：literature-check
- 记录日期：2026-11-18

## 原论文
Zimo Sheng, Mingyu Xiao. *A Polynomial Kernel for Planar Directed Feedback Vertex Set*. arXiv:2609.23669 (2026).

- https://arxiv.org/abs/2609.23669
- 原文位置：Section 7, Conclusion，尤其 Improving the kernel size / Connections to directed cut problems

## 已知结论
平面有向反馈顶点集（DFVS）要求删除至多 k 个顶点，使余图不含有向环。论文建立了 O(k^66 log^2 k) 规模的多项式核。中间对偶增广实例规模更小，返回原问题的归约显著放大了指数。

## 论文提出的方向
- 降低多项式核规模，尤其优化从中间压缩实例返回 PDFVS 的参数代价。
- 与平面有向割、终端分离有关的推广。
- 一般有向 DFVS 是否有多项式核是更长期的问题。

## 切入
重写 GSCA→3-CNF→平面化→PDFVS 全部参数传递，独立核对平面保持、正确性和指数来源；尝试替换一段归约并严格证明其参数改善。

- [ ] 检查2026年更新版本及引用
- [ ] 拆分 k^16、k^32、k^66 指数来源
- [ ] 给出可复核的改进定理或反例
