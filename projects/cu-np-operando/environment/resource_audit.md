# Resource and evidence audit / 计算资源与既有证据审计

Date: 2026-09-20. This report records live checks, not installation-time promises.

## 实际可用资源

本机为 Intel i7-13700H，14 核、20 逻辑线程，Windows 可见内存约 15.73 GiB；首次探查可用约 2.67 GiB，随后降至约 2.48 GiB。显卡为 Intel Iris Xe，没有检测到 CUDA 显卡。必须按**当前可用内存**约束新计算，不能以总内存作为可独占预算。

实际导入通过：Psi4 1.11（phase7）、RDKit 2026.03.5、NumPy 2.4.6、SciPy 1.18.0（phase2ff）。xTB 可执行程序报告 6.7.1。`def2-svpd.gbs`、`def2-svpd-ri.gbs` 和 `def2-universal-jkfit.gbs` 已存在，无需安装。

当前 PATH/已登记 Conda 环境搜索未找到 Quantum ESPRESSO、JDFTx 或调度器 sbatch/squeue/qsub。phase2ff/phase7 中 GPAW 不可导入；ASE 的 `gpaw_out.py` 只是输出解析器。phase7 未安装 PySCF。系统存在 `wsl.exe`，但列举命令返回用法文字，没有核实可运行的 Linux 发行版。

SSH 客户端存在，但用户 `.ssh/config` 与 `.ssh/known_hosts` 均不存在。没有获得可核实的 HPC 主机、调度器、额度或预算；没有尝试连接未知主机，也没有读取/展示任何私钥或凭据。这一结果说明**本任务尚无可使用的已验证 HPC 路径**，不能推断用户没有任何外部账户。

## 既有 Cu–N–P 项目可复用内容

源项目：`new-chat/work/research-repo/projects/cu-np-electroreduction`；源仓库提交 `6197929f41074cb4b21fb2792c7dcd55772a86c0`。本次重新执行文件哈希验证，20 次 xTB 原生记录及对应实际运行源码全部一致；31 项测试全部通过，测试运行 2.553 s。测试通过只说明软件与存档自洽。

分子是 3-methylquinazolin-4(3H)-one，C9H8N2O，SMILES `CN1C=Nc2ccccc2C1=O`。其坐标由 RDKit 独立生成后经 GFN2-xTB 优化，不是正文/SI 的 Cu 位点或 DFT 坐标。初始结构存在 −50.49 cm⁻¹ 振动，旧项目已保留并拒绝；甲基转子修复后的 54 个振动模均为正，最低 78.64 cm⁻¹。认可的几何 SHA256：`29351eee2bde9d27f5d9857c9d24d5ed74faa327a9e0a69ed52117ae34e688b9`。

GFN1/GFN2 同核几何充电能分歧为气相 0.781648 eV、平衡 ALPB-DMF 0.831980 eV，超过旧项目预设 0.20 eV 门槛。这个负结果和原始日志可复用；不能把它们当作还原电位、电子转移势垒、真实溶液垂直电子亲和能或 Cu 催化剂排序。单分子预筛没有验证任何实际 Cu–N–P 活性位点。

## 有价值的本地下一步

按顺序运行单个小分子的固定核 DFT 单点：先 PBE0/def2-SVPD 中性单重态作短基准，随后同核几何负离子双重态；资源允许再加入 B3LYP 同基组的中性/阴离子对。每作业 2 线程，Psi4 内存设置 512 MiB，监控进程树 RSS，达到 1 GiB 或 300 s 即终止。并行作业数为 1；保留未收敛/超时/内存终止证据。

这可回答“xTB 充电能方法敏感性是否在两种常见杂化 DFT 下仍然显著”，并提供可测量的时间/内存基准。它不能解决电极、参考电位、溶剂非平衡极化或 Cu 位点排序。def2-SVPD 包含弥散函数，但单一基组不能证明阴离子束缚性或消除自相互作用误差；即使 SCF 收敛和轨道能为负，也要保留该限制。没有优化或频率计算时不得称为 DFT 极小值。

## Precise execution boundary

The existing tools support bounded molecular diagnostics and software validation. They do not establish a periodic constant-potential interface workflow. Reproduction of CuN4/CuN3P1 chemistry still requires source-supported structures, complete reaction conditions, a verified periodic engine/HPC allocation, calibrated electronic-structure benchmarks and independent experimental/flow data. A portable submission package is an unexecuted workflow, not a completed HPC calculation.

The companion JSON records live versions, prior hashes, negative results and exact resource constraints. New quantum results, when available, belong under `dft_pilot/runs/` and must be interpreted separately from the old xTB Hamiltonians.
