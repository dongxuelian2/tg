# NWPU-001 · 有向图拉普拉斯谱矩极值

- 学校 / 研究者：西北工业大学 · 王力工（Ligong Wang）
- 领域：谱图论
- 证据类型：author-stated
- 当前状态：literature-check
- 记录日期：2026-11-18

## 原论文
Xiuwen Yang, Hajo Broersma, Ligong Wang. *Sharp bounds for Laplacian spectral moments of digraphs with a fixed dichromatic number* (2023).

- https://arxiv.org/abs/2305.06362
- 原文位置：Section 4, Concluding remarks，Problem 4.1、4.2、4.3

## 对象
设 L(G) 为论文规定的有向图拉普拉斯矩阵，特征值 λ_1,...,λ_n，则第 k 谱矩 M_k(G)=Σ λ_i^k。给定有向色数（dichromatic number），寻找谱矩的极值图。

## 作者提出的三个问题（意译）
- Problem 4.1：对固定有向色数的 join digraphs，刻画第三谱矩的全部极值图。
- Problem 4.2：推广为 k≥4 的谱矩。
- Problem 4.3：推广到固定有向色数的一般有向图、k≥3。

## 先攻范围
先处理 Problem 4.1：固定 r=3 与小 n，将谱矩写成分块规模 n_i 的整数多项式优化；重现实验表2并尝试证明极值块分配规律。对矩阵 L 的方向约定须以论文为准。

- [ ] 排查2023–2026后续文献
- [ ] 重现 r=3 的小规模极值
- [ ] 提出并检验可归纳的结构命题
