# NUAA-001 · 突发障碍下的事件触发多机编队避碰

- 学校 / 研究者：南京航空航天大学 · 陈谋（Mou Chen，2018年国家杰青）
- 领域：非线性控制、多智能体系统、无人机避碰
- 证据类型：future-work（论文明确陈述的未来工作；其中编队切换已存在后续研究）
- 当前状态：literature-check
- 记录日期：2026-11-18
- 优先试点日期：2026-11-19

## 官方身份
- 南京航空航天大学个人主页：https://faculty.nuaa.edu.cn/cm/zh_CN/index.htm
- 国家杰青名录：https://faculty.nuaa.edu.cn/jcrcjs.jsp?id=2202&lang=zh_CN&urltype=tsites.HonorTeacherList&wbtreeid=1002

## 原论文与作者提出的后续任务
Hongzhen Guo, Mou Chen, Mihai Lungu, Baomin Li.
*Distributed event-triggered collision avoidance coordinated control for QUAVs based on flexible virtual tubes*.
Chinese Journal of Aeronautics 38(2), 103300 (2025).

- DOI：https://doi.org/10.1016/j.cja.2024.11.010
- 原文位置：Section 5 / Conclusion 的 Future Work 句。
- 原文已实现：柔性虚拟管（virtual tubes）、事件触发控制、扰动观测、闭环信号有界性分析。
- 作者提出：研究降低通信负担的新事件触发机制；通过路径重规划、编队缩放、编队切换处理突发障碍物。

## 已发现的后续研究：必须先排雷
Hongzhen Guo, Mou Chen, Shuang Shi, Zengliang Han, Mihai Lungu.
*Distributed Coordinated Control for QUAVs With Switching Formation Strategy*.
IEEE Transactions on Intelligent Vehicles 10(6):3841–3851 (2025).
DOI：https://doi.org/10.1109/TIV.2024.3466152

该论文已经研究了任务/环境变化下的编队切换，故不能将“研究编队切换”本身当作尚无人做过的新成果。应进一步核查其与柔性虚拟管、突发障碍、事件触发通信、严格避碰不变量的结合情况。

相关2026论文：*Event-Triggered Distributed Formation Control for QUAVs Based on Zero-Sum Differential Game and Adaptive Dynamic Programming*，https://doi.org/10.1002/rnc.70574 。这一论文研究事件触发编队控制，也需比对其适用场景。

## 建议收缩后的候选数学命题（derived，尚未确认新颖性）
对于 N 架四旋翼的简化二维双积分器模型，给定时间变化的障碍区域和有限通信事件：
1. 在虚拟管边界允许在线缩放与切换时，构造分布式控制律；
2. 给出保证两机间及机—障碍之间最小安全距离的充分条件；
3. 在有界外扰和通信触发误差下证明安全集的前向不变性；
4. 检查排除 Zeno 行为所需的触发阈值与假设。

第一步先复现原文稳定性与避碰条件；第二步查明2025年编队切换论文已解决的情形；第三步才考虑新定理。不要先堆叠神经网络或设计大型模拟器。

## 2026-11-19 交付目标
- 精确的原文及后续论文比对，标明已经解决的部分；
- 一个定义完整、尚未被已知论文覆盖的候选命题，含假设和失败情形；
- 若有研究突破：独立核验的证明、反例或可复现实验，以 PR 提交；
- 如尚无突破，如实记录最新文献状态及阻碍，不宣称已解决。
