# HEBUT-001 · 六环图的饱和数精确值

- 学校 / 研究者：河北工业大学 · 兰永新（Yongxin Lan）
- 领域：极值图论
- 证据类型：gap-in-result
- 当前状态：literature-check
- 记录日期：2026-11-18

## 原论文
Yongxin Lan, Yongtang Shi, Yiqiao Wang, Junxue Zhang. *The saturation number of C6*. Discrete Mathematics (2025).

- https://doi.org/10.1016/j.disc.2025.114504
- https://cfc.nankai.edu.cn/_upload/article/files/d0/e0/14bed57849dab44905b374017424/8730f932-22b5-4963-9f2e-fafd6e88338b.pdf
- 论文定位：Theorem 1.5

## 定义与已知
G 是 C6-饱和图，当且仅当 G 不含六环，添加任意一条缺失边都会产生六环。sat(n,C6) 为 n 顶点这种图的最少边数。原文给出 n≥9 时的统一界：

\[4n/3-2 \leq \operatorname{sat}(n,C_6) \leq (4n+1)/3.\]

按 n mod 3 的更精细构造见原文 Theorem 1.5。

## 挑战
确定足够大 n 的精确公式，改进上下界，或刻画达到极值的饱和图。首先筛查2025–2026后续成果。

## 验证与推进
- 以 SAT/SMT/整数规划或带同构消除的搜索计算小 n；每个输出独立检查无 C6、添加任一缺失边产生 C6。
- 从 n mod 3 分组观察结构；检验原文下界取等条件。
- 证明推广时明确图为简单无向图、n 的范围与所有取整细节。

## 未完成
- [ ] 后续文献核验
- [ ] 构造、下界和小规模基线复现
