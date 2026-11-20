# NUAA-001 · 突发障碍下的事件触发多机编队避碰

- 学校 / 研究者：南京航空航天大学 · 陈谋（Mou Chen，2018年国家杰青）
- 领域：非线性控制、多智能体系统、无人机避碰
- 证据类型：future-work（论文明确陈述的未来工作；其中编队切换已存在后续研究）
- 当前状态为 in-progress（已有分模型结论、反例和构造，原始广义问题仍待解决）
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

以下为立项时的候选方向。2026-11-19 的研究选择保留旋翼推力分配及俯仰惯量，所得结果列于本文末节。
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

## 2026-11-19 · 已证明的衍生结果

原论文 §3.3 定理 1 在假设 1–3 及其指定的势场、观测器、量化器和事件触发律下给出信号有界、机间避碰与无 Zeno 行为。式 (21) 把管道势能、跟踪及扰动观测误差、边界调节辅助状态共同纳入 Lyapunov 函数。该定理没有给出方向信息晚到时所有因果策略的能力上界。后续针对未知通道的工作已处理在线感知、缩放与跟踪误差。[原论文](https://doi.org/10.1016/j.cja.2024.11.010)；[Zhu、Lian、Chen、Lin，2026](https://doi.org/10.1016/j.isatra.2026.10.032)

1. 对保留推力与俯仰力矩耦合的平面机群，求出一比特障碍方向消息迟到时、可保证净空的精确一阶损失。传感器对阵风的观测滞后会改变这一损失。给出达到界的分布式规则；两组风史对任何因果策略给出匹配的不可能性界。[定理、证明及新颖性对照](NUAA-001/delay-price.md)
2. 固定水平机动时长、飞行器推力界、移动墙的速度和加速度界，构造任意整数阶的一族实例。它们有正鲁棒净空，却能迫使所有安全的分段常值控制至少切换指定次数；无风时只需要一张单调、凸的墙面位置曲线。另一方面，给定净空精度，可用脉冲矩匹配得到最大保底净空的二阶夹逼。对有可计算连续模的已知可计算墙面，给出对严格正负两侧均终止的穷举判定法。[构造、唯一性证明、有理数证书及计算界](NUAA-001/moving-wall.md)

在仓库根目录运行以下命令可复核。脚本使用 NumPy 和 SciPy，无网络依赖。此次使用 Python 3.13.1、NumPy 2.2.5、SciPy 1.16.1。两旋翼组的 32 段独立直推力优化允许反向倾转，以三个随机初值搜索反例；另核对时域差分、无重置的三个检查区间、24 个单机风史树线性规划、128 条双机联合风史及非线性前后段连续动力学积分。固定归一化 $A=B=H=1$ 时，独立时间差分与积分给出 $K=0.129636212$；双机线性规划的最坏残差为 $0.91$。

~~~sh
python problems/NUAA-001/verify.py
~~~

移动墙的反例搜索取两个检查点、16 格切换时间和 5 档连续控制幅值的离散子集，穷举 13,125 组至多两次切换的指令。最好的候选仍在至少一个检查点不足 $2.2675\times10^{-5}$（归一化位移）。这个数字只描述所搜网格；无网格假设的不可能性由唯一性和紧性证明。脉冲近似另在四种网格上核查格点的角度、角速度及格间位移误差。

该结果限于短时预警、指定通信和观测结构、标明的空间约束。对一般三维障碍、固定电机带宽及任意次重规划的无限时域递归可行性，尚无完整判据。

## 2026-11-20 · 阶段研究成果

本次整理收录 58 份研究笔记与 33 个核验脚本。各笔记列出适用的动力学、观测和执行约束。部分笔记包含待核实文献优先权的命题或仍需补证的工作推导。原问题状态仍为 in-progress。[全部阶段材料](NUAA-001/)

- **三维动力学、重规划与换形。** 推导三维刚体和电机状态的轨迹拼接与可达尺度，以及缩放速度、控制余量和换管净空的判据。[三维跟踪](NUAA-001/three-dimensional-tracking.md) · [编队时间律](NUAA-001/formation-clock.md) · [换管证书](NUAA-001/tube-handover.md) · [短时悬停交接](NUAA-001/hover-handover-width.md)
- **通信延迟与私有信息。** 给出接收历史和剩余控制能力的联合条件，以及两机纠偏的负环证书。三机反例表明，在固定不可行差距下，任意给定数量的情景也可能不足以证伪。另给出一个可靠的一比特广播能够恢复可行性的模型，并核验其有限时限临界值。[擦除与截止时间](NUAA-001/erased-trigger-deadline.md) · [两机纠偏](NUAA-001/private-correction-cycles.md) · [三机反例](NUAA-001/three-agent-witness.md) · [一比特修复](NUAA-001/private-bit-repair.md)
- **扰动、遥测和遮挡。** 推导共同阵风的格间精确包络和相对观测的地速歧义，再将遮挡后验、通信擦除和真实交接状态纳入可检验的条件。[阵风消元](NUAA-001/common-wind-quotient.md) · [相对遥测](NUAA-001/drag-sampling-radius.md) · [空报告传播](NUAA-001/visibility-history.md) · [预警交接](NUAA-001/warning-handover.md)
- **无限时域与复杂性。** 在指定信息与鲁棒余量条件下给出递归可行性构造，并记录有限证书的适用边界、异步私有决策造成的困难以及需要继续证明的条件。[鲁棒递归构造](NUAA-001/nonconvex-recovery.md) · [无限时域底稿](NUAA-001/infinite-horizon.md) · [无公共刷新证书](NUAA-001/no-refresh-certificates.md) · [移动门复杂性](NUAA-001/repeated-gate-hardness.md)

使用 Python 3.13.1、NumPy 2.2.5 和 SciPy 1.16.1 运行了 33 个新增脚本及原有的 `verify.py`，均通过。脚本检验报告枚举、线性规划、精确算术证书及刚体方程的数值积分。每项结论的证明和适用边界见对应笔记。可以在仓库根目录运行以下命令复核。

~~~sh
python problems/NUAA-001/verify.py
python -c "from pathlib import Path; import subprocess, sys; [subprocess.run([sys.executable, str(script)], check=True) for script in sorted(Path('problems/NUAA-001').glob('check_*.py'))]"
~~~
