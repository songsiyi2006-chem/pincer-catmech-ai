# 独立代码审查与修复复核

日期：2026-09-20。审查者：独立 literature_audit agent。审查全程未修改被审文件；复现脚本、临时证据和日志写在本审计目录。

初审对象：`work/kinetics_engine/src/neutral_transport/core.py` 与 `work/active-repo/projects/cu-np-operando/scripts/evidence_gate.py`。修复由主代理/实现代理完成。最终复查核心为准备交付的 `work/active-repo/projects/cu-np-operando/src/neutral_transport/core.py`；原 agent 工作目录保持冻结，不能误把其中未更新的副本当最终交付。

## 结论

在本次覆盖范围内，发现的科学门控问题和输入一致性问题已经在 **active-repo 最终副本**闭环。最终已执行并通过 **20 个 neutral-transport 测试与 8 个 evidence-gate 测试**。这支持工程正确性检查，不证明目标 Cu 电还原机理、能垒、排名反转或工业性能。

## 可操作发现与修复状态

| 编号 | 初始严重性 | 可复现问题 | 最终状态 |
|---|---|---|---|
| R1 | P1 | `evaluate` 仅依赖自报标签与哈希。9 条要求均指向同一个明确写着 `SYNTHETIC_UNIT_FIXTURE` 的本地文本，只要贴上 `accepted`、`CALCULATION` 与非空 reviewer/scope，就曾返回 `physical_ranking_ready=true`。这把元数据完整性暴露成了可用于解锁科学预测的布尔值 | **已修复**：完整性改用 `evidence_manifest_complete`，`physical_ranking_ready` 固定 false，注明尚无经验证的目标反应模型。相同重贴标签探针当前只可能使完整性为 true，不能解锁科学排名 |
| R2 | P2 | `gate` CLI 在要求全部不通过时仍正常退出 0，下游按退出码执行的脚本可误认为门槛通过 | **已修复**：被拒绝时退出码 2。以空要求文件在隔离临时目录实际运行，返回 2 |
| R3 | P2 | `requirements: null` 或某条 entry 为布尔值时触发 `AttributeError`，而非完整的拒绝报告 | **已修复**：当前统一产生不通过的结构化结果，且继续封闭物理排名；新增 malformed 测试通过 |
| R4 | P2 | 手工或导入 `RateConstants` 可以绕过生成工厂：`RateConstants(1,1,1,1,1,1,-1,100)` 在初始 `steady_surface` 中被接受，尽管平衡常数为负、循环残差巨大 | **已修复于交付副本**：要求 Keq 为正，按六个速率重算对数循环并与 log(Keq) 比较，校验所报残差。上述探针当前抛 `ValueError`；三类伪造速率元数据测试通过 |
| R5 | P2，低现实发生率的数值防御问题 | 有限但极端电压/电流相乘溢出后，初始 `cell_energy` 可返回无限能耗 | **已由实现代理修复并复核**：乘法/积分使用 NumPy 错误策略，输出也校验有限性。极端有限输入当前明确拒绝 |

R1 的修复没有“认证”外部材料的真实性。当前 `evidence_manifest_complete=true` 只表示声明字段/文件完整，不能转述成证据可信或科学验收通过。将来引入真实目标模型时，须重新设计该模型自己的参数、域、路径连通与前瞻验证门槛。

## 科学、单位与数值核验覆盖

- 人工检查 3 态单个位点循环的 matrix-tree 权重，三个边通量符号及位点归一化公式；未发现所写公式的代数错误。
- 确认所有吸附速率乘以无量纲活度 `c/c_standard`。速率单位为 s^-1，位点密度为 mol sites/m²，面通量为 mol/(m²·s)。PFR 源项由几何反应面积、体积流量和面通量组合，量纲一致。
- 检查薄膜通量可行区间 `[-km_P*cP, km_A*cA]`、逆反应负通量、零位点与零体积极限；该传质近似被准确标为中性物种双膜模型，并未冒充 Nernst-Planck 或真实电解槽。
- 测试覆盖独立线性方程求稳态、对称速率解析解、薄膜解析解、传质渐近极限、PFR 一阶网格收敛、Radau 独立积分与源项积分物料平衡。
- 确认产品进口量不被计入净生成量；能耗使用全槽 `∫UI dt` 并按分离质量归一化；FE 接受的是分离前生成摩尔数与已知电子计量。当前中性网络不能提供该电子计量。
- 自建平衡态探针在约 10^11 s^-1 的快正逆速率下，残余净 TOF 约 6.1×10^-12 s^-1，边通量差约 1.1×10^-11 s^-1，属于该例的消减舍入量级。本次未将它解释成非零化学反应，也未声称穷尽所有刚性参数范围。
- `audit_group_split([])` 返回“无冲突”是空集上的一致性结果；现有 scope 已写明仅检查声明的 family。它不能证明数据量足够、划分科学合理或模型无泄漏。当前没有使用它放行 AI 性能主张。

## 实际执行记录与环境问题

最终指定环境为 `C:/Users/HUIWEI/miniconda3/python.exe`。只读执行结果：

| 文件 | 结果 | unittest 计时 |
|---|---|---|
| `tests/test_neutral_transport.py` | 20/20 通过 | 0.559 s |
| `tests/test_evidence_gate.py` | 8/8 通过 | 0.039 s |

此前直接使用 `phase2ff/python.exe` 与工作台 venv 时，运行至 `numpy.linalg.solve` 出现原生退出 `-1066598273`，未生成 Python traceback。对 phase2ff 单进程补齐 `Library/bin` 的 PATH 与 `os.add_dll_directory` 后，同一 20 项测试通过；随后按照主代理要求改用 base 环境并再次通过。因此将失败归为**该启动方式下的本地 DLL 环境问题**，不据此判定代码失败。未安装、升级或覆盖任何环境。

日志：`neutral_tests_review.log`、`neutral_tests_review_chem_ai4s.log` 保留失败；`neutral_tests_review_dll_path.log` 保留恢复验证；最终采用 `neutral_tests_review_final_base.log` 与 `gate_tests_review_final_base.log`。

## 最终审查哈希

SHA256 见 `reviewed_code_hashes.json`；以下是本轮复查时读取的代码字节，不是签名或真实性认证。

| 文件 | SHA256 |
|---|---|
| 冻结的 agent 核心 `kinetics_engine/src/neutral_transport/core.py` | `8fde79a7e5b4b87cf834dd882e0b15eca91c3026fedc8be6840fa746bdfbde8d` |
| 交付核心 `active-repo/projects/cu-np-operando/src/neutral_transport/core.py` | `e2c9ff20b97f75b891836f9568f5f7022e332c47642f142e5f249557550a2a44` |
| 交付门控 `active-repo/projects/cu-np-operando/scripts/evidence_gate.py` | `385bb61d92c5613d6dd92335310f89c2f0d3b9b277f3e4709343b75740a33168` |

## 局限

这是两个核心模块、相关测试与定向异常输入的独立审查，不是全仓库安全审计、形式证明、完整浮点域验证或科学同行评审。没有检验真实 Cu 位点、反应结构、恒电位电子结构、真实流动实验、寿命或 TEA 输入。`probe_review.py` 是工作区复现辅助脚本，其路径针对当前审计目录，不应作为面向用户的项目入口。
