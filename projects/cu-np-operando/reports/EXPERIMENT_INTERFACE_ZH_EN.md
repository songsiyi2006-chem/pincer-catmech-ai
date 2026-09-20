# 前瞻实验接口 / Prospective experimental interface

**状态：未执行；未提供定量实验预测。Status: unexecuted; no numerical experimental prediction registered.**

当前缺少完整反应映射、校准电位、真实 Cu 位点动力学与连续流数据。因此本文件
定义将来应收集什么，而不伪造已被本项目验证的催化剂性能。

The present isolated-molecule DFT and neutral numerical fixtures do not supply
Cu-site rate constants. The following discriminating measurements remain
conditional on a physically calibrated target-reaction model.

| 问题 / Question | 测量 / Measurement | 关键对照 / Controls | 反证 / Contrary outcome |
|---|---|---|---|
| 产物是否竞争占位？ / Product inhibition? | 固定电位下的初始速率与添加产物浓度 / Initial rate versus added product | 分析回收率、溶液导电性、位点密度、相同电极 / Recovery, conductivity, matched electrode and sites | 效应小于可分辨误差 / No resolvable effect |
| 流量是否改变排序？ / Flow-dependent ranking? | 匹配CuN4/CuN3P1的流量—浓度扫描 / Matched flow and concentration scan | 局部电位、几何面积、负载、产物分离与停留时间 / Potential, area, loading, isolation and residence time | 在稳定工况内无可信交叉 / No credible crossing |
| 活性中心是否保留？ / Active-site retention? | 运行时间序列、溶出与前后结构 / Time series, leaching and structure | 空白载体、浸出液、独立合成批次 / Support blank, filtrate and independent batches | 重构/均相活性解释更好 / Reconstruction or dissolved species explain activity |

实验记录必须含原始色谱/谱图、校准曲线、检测限、内部标准、物料平衡、电极尺寸、
催化剂负载、参比电极与iR处理、温度、完整电流/槽压时间序列、分离收率与批次标识。
不能用出口浓度直接代替分离产率，不能用分离收率直接代替 FE。

Retain raw analysis, calibration curves, detection limits, internal standards,
material balances, electrode dimensions/loading, reference and iR treatment,
temperature, full current/cell-voltage time series, isolation yield and batch IDs.
Randomize run order, use independent batches, and determine repetition from the
effect size and variance. Lock numerical predictions, tolerances and failure
criteria before accessing held-out results.

No chemical experiment, researcher contact or data-sharing action has been
performed by this computational delivery.
