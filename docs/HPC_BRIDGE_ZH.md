# 从本地研究包到超算与专业量化软件：人类协作操作指南

本指南与 `prepare_hpc_campaign.py`、`collect_hpc_results.py` 及 `pincer_catmech.hpc.bridge` 配套。当前用户尚未开通超算，预计可能在大二前获得许可；因此当前工作包括离线输入生成、核验和申请准备，**没有建立远端连接、没有提交付费作业，也没有得到 ORCA 计算结果**。将来取得账号时需重新核实当时的版本、费用与平台规定。

已有 24 个条目是计算设计目录，不是 24 个经实验确证的催化剂。本桥接对默认设计结构标记 `computed design geometry; not experiment-matched`；若采用另行核验的公开晶体结构衍生候选，则显式标记 `public_crystal_derived_candidate` 并保留出处与转换过程。晶体来源不自动证明其为溶液活性物种或确定基态自旋。同一 XYZ 上的不同自旋单点只属于垂直自旋诊断。计算成功不等于结构合成可行、真实活性物种确认、自由能面验证或达到某一期刊标准。

## 1. 现在没有超算账号，先完成什么

先找愿意实际指导的计算化学或有机金属催化老师。带去一页研究问题、目前失败与成功记录，以及本指南；讨论应优先研究哪一个可以联系到文献或实验体系的问题。不要以“已经有 24 个高效催化剂”或“Nature 级 AI 模型”介绍项目。

向老师说明需要的具体支持：确认分子结构、电荷及可能的自旋态；选定一个足够小的基准体系；指导量化方法误差与实验对照；通过课题组或学校合规申请计算资源。学校计算中心、学院共享平台、老师现有超算项目账号都可能是渠道，但本项目未核实任何一方已经承诺资源。没有授权时不要借用他人账号或转交凭证。

账号申请与研究设计可并行推进。现在能够做的包括学习 Linux 文件操作和 Slurm 日志、核对公开文献原始结构与使用许可、复现一个小分子的输入输出流程、整理基准与独立测试划分。没有账号期间不必购买大额算力，也不需要在个人电脑上反复堆积无法收敛的大计算。

## 2. 需要从导师或管理员得到的信息

填写 `configs/hpc/site_information.template.json` 的非敏感信息，分享给研究协作者即可。**不要填写或发送密码、私钥、访问令牌、身份证件或私人 SSH 配置。** 最重要的字段如下。

| 类别 | 必须核实 | 用途 |
|---|---|---|
| 调度 | 是否 Slurm，账户、分区、是否额外需要 QOS、数组上限 | 本包只实现 Slurm 单节点；其他调度器需要适配 |
| 资源 | 每节点核数和内存、最长时限、同时运行作业数 | 防止超出配额和内存不足 |
| 软件 | 获合法授权的 ORCA 确切版本、完整可执行路径、配套 MPI 与模块命令 | 输入语法、MPI 启动与结果可追溯 |
| 运行时 | Python 3.10 或更新版本的完整路径 | 桥接工具只需 Python 标准库 |
| 存储 | 共享项目目录、节点本地 scratch、配额、回收规则 | 计算在新 scratch 目录执行，回传原始文件 |
| 计费 | 按核时、节点时、内存时还是套餐收费，最低计费粒度 | 首次只做一个 pilot，测量后再申请更多预算 |
| 责任 | 指导老师和管理员认可的运行方式、预算上限 | 明确专业判断与资源使用责任 |

如果 ORCA 是平台提供的软件，使用管理员给出的模块，不自行替换 MPI。若需安装，由获得授权的使用者按官方渠道及机构许可办理；脚本不会下载或安装 ORCA，也不包含商业软件二进制。Gaussian 同样必须通过合法许可和管理员提供的安装使用。

## 3. 当前桥接协议到底算什么

默认且目前唯一实现的量化协议是 **ORCA 6.x / r2SCAN-3c / CPCM(toluene) / TightSCF / SP**。`r2SCAN-3c` 是 ORCA 官方提供的组合方法，内部包括专用基组与修正，不能再随意加一套基组或叠加 D3/D4 后仍称原协议。[ORCA 组合方法手册](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/3cmethods.html)

选择它是为先得到可核验的结构与自旋电子能诊断，并非认定它足以定量决定过渡金属自旋排序。单重态当前用 RKS，多重态用 UKS；尚未覆盖破缺对称单重态、多参考态、轨道稳定性或不同初猜的系统搜索。缺少这些检验，不能把一个成功单点称为真实基态。

**CPCM 不是 SMD，也不是先前 xTB 的 ALPB。** 输入写 `CPCM(toluene)` 不会自动获得 SMD 的完整溶剂模型。SMD 有独立启用方式和参数，改变模型必须建立新的协议及运行目录，不能混在同一能量表中。[ORCA 隐式溶剂手册](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html)

383.15 K 是反应参考温度，**当前 SP 不计算振动热化学，不把电子能当作该温度的自由能**。ORCA 内置甲苯参数也没有在本包中进行 383.15 K 的介电常数校正。5 mol% t-BuOK、0.01–0.20 当量扫描以及 t-BuOH 活度是后续物种平衡与动力学问题；单点输入没有凭空加入它们。

桥接程序校验元素组成、元数据分子式、电荷和电子数与多重度奇偶性，但这些检查不能证明键连关系、氧化态或催化身份正确。输入前仍需专业人员查看三维结构和化学模型。

## 4. 本地准备与离线核验

准备一个显式的选择文件。路径相对 `--source-root`，不能包含 `..`、绝对路径或反斜线。下面是已有计算设计几何的示例，不是实验结构认证。

```json
{
  "schema": "pincer_hpc_selection_v1",
  "records": [
    {
      "id": "Fe_bipyridine_pnnoh_iPr_active_M1",
      "xyz": "data/structures/Fe_bipyridine_pnnoh_iPr/active/best_found.xyz",
      "metadata": "data/structures/Fe_bipyridine_pnnoh_iPr/active/metadata.json",
      "charge": 0,
      "multiplicity": 1
    }
  ]
}
```

在仓库目录下，用已有工作环境运行；示例资源用于说明格式，实际值须由管理员和 pilot 决定。

若使用公开晶体衍生候选，在该 record 中增加 `"structure_origin": "public_crystal_derived_candidate"`，并在原 metadata 中提供 `provenance` 对象，其中 `source_url`、`source_identifier`、`structure_derivation` 必须为非空字符串；可以追加 DOI、许可、原 CIF 哈希及氢原子/离子/溶剂处理记录。生成器保存其原始信息，不把来源标签自动升级为经过核验的溶液几何、自旋态或催化身份。`experimental_identity_verified`、`solution_geometry_verified`、`ground_spin_verified` 仍为 false。

```text
python scripts/prepare_hpc_campaign.py --selection selection.json --source-root . --output work/hpc_bundle_001 --nprocs 4 --memory-mb 8192 --maxcore-mb 1000 --walltime 01:00:00 --concurrency 2
python work/hpc_bundle_001/bridge.py verify --bundle work/hpc_bundle_001
```

所有资源参数显式必填。程序要求 `nprocs × MaxCore ≤ 0.75 × memory_mb`，为主进程、MPI 和额外工作区留下余量；这是本项目采用的保守上限，不是内存绝不会超出的保证。ORCA 的 `%maxcore` 是每进程主要工作区设置，并非整个作业硬性峰值上限。[ORCA 内存说明](https://www.faccts.de/docs/orca/6.1/tutorials/first_steps/memory.html)

生成的新目录包含 `manifest.json`、可移植 `bridge.py`、`array.sbatch`、每个任务的 `job.inp`、原始 XYZ 与原始 metadata。每份文件记录 SHA256。已有目录不可覆盖；重试、新版本、不同溶剂或方法必须另建目录。离线验证只证明包内部一致，**不证明 ORCA 已安装或任务可收敛**。

## 5. 取得账号后的人工桥接

用平台允许的 SFTP、SCP、WinSCP 或其他方式把完整包送至共享项目目录。登录与双因素验证由你自行完成，无需向 AI 提供凭证。在超算上先加载管理员指定的软件模块，再重新执行包内的离线 `verify`，比较本地和远端 `manifest.json` 的 SHA256。

需要在超算 shell 中显式设置以下值；尖括号只是说明，必须换成真实配置后才能执行。

```bash
export HPC_ACCOUNT='<authorized-project-account>'
export HPC_PARTITION='<cpu-partition>'
export ORCA_EXE='<absolute-path-to-orca>'
export HPC_PYTHON='<absolute-path-to-python-3.10-or-newer>'
export HPC_SCRATCH='<existing-node-local-scratch-parent>'
```

scratch 路径若只在计算节点存在，不要在登录节点用不存在的路径强行通过检查。先请管理员确认合适的节点分配或共享 scratch 约定；当前提交预检要求该父目录在提交环境可见。这个限制需要按站点情况适配，不能把临时绕过当作已完成部署。

从作业包目录开始，先查看即将提交的参数，再执行 **一个任务**。

```bash
"$HPC_PYTHON" bridge.py submit --bundle . --mode pilot
"$HPC_PYTHON" bridge.py submit --bundle . --mode pilot --execute
```

不带 `--execute` 只打印参数。带它才在当前已登录且获授权的机器上调用 `sbatch`。提交器不读取私人 SSH 配置，不自动连接其他主机。账户和分区不允许默认猜测。若站点必须提供 QOS、约束或特殊 MPI 启动参数，先请管理员审核并适配后再用，当前版本不会悄悄填入。

单节点数组申请 `--ntasks=nprocs`、`--cpus-per-task=1`，输入中 `%pal nprocs` 与之对应，并设 OMP/MKL/OpenBLAS 线程数为 1。程序检查实际 Slurm 分配的核数和内存，发现超配即拒绝运行。**ORCA 主程序直接以完整路径调用，不能再用 `mpirun orca` 包裹；由 ORCA 驱动启动并行模块。** MPI 版本和平台集成必须由管理员验证。[ORCA 并行运行手册](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html)

用平台工具查看队列和账单记录。例如，在相应平台确认支持后，记录 `sacct` 中的 JobID、State、ExitCode、Elapsed、AllocCPUS、TotalCPU、MaxRSS。调度器显示 `COMPLETED` 并不自动等于量化结果有效，必须检查原生输出的任务状态。

## 6. Pilot 验收与费用估算

你与导师首先检查 `results/<task-id>/run.json`、`native/job.out` 和 `native/job.err`：确切 ORCA 版本、实际方法及基组、溶剂设置、电荷/多重度、SCF 是否收敛、最后电子能、S²、峰值内存与墙钟时间，以及结构身份。正常结束但 S² 缺失时，自旋质量仍记为 unknown，不自动通过。

只有 pilot 运行完成、正常收敛且人工判断资源与模型合理时，才写入审查记录。例如 `pilot_review.json` 包含：

```json
{
  "approved": true,
  "reviewer": "真实审核者姓名或项目内署名",
  "reviewed_utc": "真实审核时间",
  "manifest_sha256": "当前manifest.json的真实SHA256",
  "pilot_run_sha256": "pilot run.json的真实SHA256",
  "notes": "结构/模型/资源/预算核验说明；缺失项及限制"
}
```

这是一份人工责任记录，不是数字签名或防伪认证。不要复制本示例直接声称已经审查。程序核对两个哈希与 pilot 的数值收敛状态，随后才允许提交剩余数组；它不会代替专家判断。

```bash
"$HPC_PYTHON" bridge.py submit --bundle . --mode production --pilot-review pilot_review.json
"$HPC_PYTHON" bridge.py submit --bundle . --mode production --pilot-review pilot_review.json --execute
```

数组排除已完成的第 0 个 pilot，按生成包时的并发上限运行。Slurm `%N` 语法表示同时运行任务数上限；它不减少总任务量或总费用。[Slurm 数组文档](https://slurm.schedmd.com/job_array.html)

若平台按分配核时收费，粗估费用为 `任务数 × 每任务核数 × 每任务计费小时 × 每核时价格`。4 核、实际 0.5 小时的单次 pilot 是 2 核时；如果按申请时限、整节点或最低计费块结算，就必须按平台规则重算。这只是算式示例，并非本项目测得耗时或实时报价。用 pilot 实测记录估计后续范围、失败重试预算和存储成本；先取得额度，再扩批。

## 7. 需要返还什么，如何核验

保留整个包和 `results/`，包括失败的原生文件、各任务 `run.json`、Slurm 标准输出/错误、`manifest.json`、人类 pilot 审查记录，以及 `sacct` 资源记录。不要只发截图或最终能量。脚本在独立 scratch 中执行后将原生文件全部复制回结果目录，保留输入、输出、波函数与其他产生的文件；它不会自动删除 scratch，也不假装存储开销为零。

如果作业被节点故障、硬性时限或 OOM 强制终止，`run.json` 可能停留在 `running`，部分文件留在节点 scratch。这样的记录不能收集为完成任务。请及时联系管理员按平台保留规则恢复原始文件，记录终止原因，并在新的包或任务目录中重试，不能手改 `executed=true`。

在超算上或下载完整目录后执行：

```bash
"$HPC_PYTHON" bridge.py collect --bundle . --output ../hpc_return_001
```

本地等价入口：

```text
python scripts/collect_hpc_results.py --bundle path/to/returned_bundle --output work/hpc_return_001
```

收集程序逐文件检查输入、输出、几何及运行回执的 SHA256，保留正常执行但 SCF 失败的诊断结果，不把它们筛掉；未返回的任务标记 `not_returned`。缺少执行状态、输出哈希不符、版本缺失、非有限数值等会阻止建立有效记录。没有任何任务真正执行时，不能创建“结果包”。输出为 `collection.json` 与 `hpc_return.zip`，ZIP 每个成员另有 SHA256 并核对 CRC；收集器只读取目录，**不会自动解压来路不明的归档**。

哈希证明的是传输和版本一致性，不能证明科学结论正确，也不能证明远端操作者没有伪造文件。需要研究者保留平台作业记录、原始数据和可复算条件。单元测试的原生样式输出明确标记 SOFTWARE FIXTURE，默认收集器拒绝把它们当研究结果。

## 8. 何时升级到专业机理计算

先用已有实验或公开可靠文献确定研究对象与反应条件，再按专家审查逐级升级：活性物种与离子对/构象平衡；关键几何优化和 Hessian；竞争反应路径及过渡态频率与 IRC；合理方法、基组和溶剂处理下的稳健性；热化学与标准态修正；必要时的多参考性质、破缺对称、MECP 和自旋轨道耦合。每一步有单独验收条件。

这些后续任务尚未由当前桥接器自动生成或认证。NEB 最高点、一个负频、一个自旋能量交点或单次低能单点，都不足以直接生成实验速率、TOF 或动力学参数证书。需要方法适用性、路径连通性、物种与标准态的一致证据。

如导师选择 Gaussian，先依据其已授权版本核验方法支持、基组/赝势、积分网格、SCF、溶剂模型、温度和频率设置；它不是把 ORCA 输入换后缀。相同方法名字也可能存在实现与默认值差异。当前仓库没有把 Gaussian 结果称为与本 ORCA 协议等价，也未声称两套程序交叉核验已完成。

本桥接的现实目标是把“需要算什么、输入从哪里来、如何在获授权平台运行、如何返还和追溯失败”具体化。Nature 主刊级贡献最终仍取决于发现的重要性和独立证据，软件桥接与超算额度本身不构成发表保证。
