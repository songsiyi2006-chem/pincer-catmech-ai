"""Write bilingual diagnostic notes from the validated archive summary."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
summary = json.loads((ROOT / "summary.json").read_text(encoding="utf-8"))
table = ["| Run | Method/charge | Status | Energy / Eh | Wall / s | Peak RSS / MiB |", "|---|---|---|---:|---:|---:|"]
for r in summary["runs"]:
    energy = f'{r["energy_hartree"]:.12f}' if "energy_hartree" in r else "unavailable"
    elapsed = f'{r["elapsed_seconds"]:.2f}' if r.get("elapsed_seconds") is not None else "not launched"
    table.append(f'| {r["run_id"]} | {r["method"]}, q={r["charge"]} | {r["status"]} | {energy} | {elapsed} | {r["peak_rss_mib"]:.2f} |')
table_text = "\n".join(table)
pair_rows = []
for p in summary["energy_pairs"]:
    pair_rows.append(f'- {p["method"]}/def2-SVPD: ΔE = E(anion)−E(neutral) = {p["delta_E_anion_minus_neutral_eV"]:.6f} eV; A = −ΔE = {p["electronic_attachment_energy_neutral_minus_anion_eV"]:.6f} eV. Anion alpha-HOMO = {p["anion_alpha_homo_hartree"]:.6f} Eh.')
pairs_text = "\n".join(pair_rows) if pair_rows else "No complete charge pair: ΔE, attachment energy and inter-functional spread are **unavailable**."
status_text = f'{summary["completed_quantum_runs"]} completed; {summary["failed_or_terminated_runs"]} failed/terminated; {summary["preflight_rejections_not_counted_as_quantum_runs"]} prelaunch resource rejections.'
en = f'''# Local DFT diagnostic report

**Evidence class: CALCULATION.** {status_text}

This is a bounded, real calculation on independently generated C9H8N2O. It is not a completed CuN4/CuN3P1 catalytic study. The accepted prior GFN2-xTB neutral geometry was used as fixed input; no original SI coordinates or Cu active sites are implied.

{table_text}

## Charge-pair result

{pairs_text}

Only a completed neutral/anion pair with identical nuclear hashes, basis content, Hamiltonian, DF settings, grid and convergence tolerances supplies a difference. Restricted/unrestricted spin treatments differ as required by charge state. A documented SOSCF rescue, if present, changes the iterative solver only. No last SCF iteration from a timeout is accepted as an energy label. No xTB total energy appears in the DFT differences.

## Interpretation and limits

At fixed nuclei R, ΔE = E(anion; R) − E(neutral; R), and the electronic attachment energy is A = −ΔE. Positive A means a lower anion energy within the finite-basis model; negative A means the charged state is higher. Neither establishes a calibrated physical electron affinity or an electrode potential.

SCF convergence is not SCF stability. No orbital-Hessian stability analysis was performed, and a single SAD initial guess was used for the standard runs. A second-order SCF solver does not substitute for a stability analysis. Spin expectations, if a run completes, describe the corresponding Kohn–Sham determinant. Positive anion occupied orbital energies are warnings; finite diffuse basis functions can confine a state that is not physically bound. A negative orbital energy would likewise not prove binding.

The single def2-SVPD basis and 75×302 integration grid were not convergence-tested. No DFT optimization, vibrational Hessian, non-equilibrium solvent polarization, electrode electron reservoir, potential-reference calibration, Cu coordination site, kinetic network or experiment was included. Agreement between two approximate functionals would not eliminate their shared systematic error.

## Resource and next-decision evidence

The native inputs, output, software version, spin populations, final energies, orbital coefficients (where completed), wall time, sampled memory history and hashes are preserved. Each quantum launch used two threads, a 512 MiB Psi4 memory setting, and a 300 s timeout. The process-tree RSS cap was tightened as host memory changed; no quantum jobs ran in parallel. Prelaunch RAM failures are recorded separately from native runs. A 0.1 s RSS monitor is a sampled guard, not an OS-enforced absolute memory quota.

The next defensible work is a separately budgeted SCF/basis/anion-binding diagnostic or a sourced electrode model after the missing SI/HPC requirements are resolved. These results do not support catalyst ranking, reduction potentials, selectivity or industrial productivity. See [Psi4 SCF/stability documentation](https://psicode.org/psi4manual/master/scf.html#stability-analysis) for the distinction between a converged SCF solution and a stable one.

## Reproduce

Run `python summarize_dft.py`, followed by `python -m unittest discover -s . -p "test_dft_archive.py" -v`. Rebuild the bilingual notes with `python write_summary_reports.py`. These commands require no Psi4 calculation. The computational rerun command and exclusion policy are in README.md. Rebuildable scratch is retained locally and listed by hash, but excluded from Git.
'''
zh = f'''# 本地 DFT 诊断技术说明

**证据类型：CALCULATION。** 本次归档：成功 {summary["completed_quantum_runs"]} 个原生计算，失败/终止 {summary["failed_or_terminated_runs"]} 个，启动前资源拒绝 {summary["preflight_rejections_not_counted_as_quantum_runs"]} 次。

本项工作确实运行了局部量子计算，但研究对象只是独立生成的 C9H8N2O 分子，不能代替 CuN4/CuN3P1 催化研究。坐标继承此前经过 GFN2-xTB 振动验证的中性几何，本轮全部固定核位置；不是原论文/SI 坐标，也没有 Cu 活性位点。

{table_text}

## 中性/阴离子配对结果

{pairs_text}

只有原生输出明确收敛，且中性/阴离子具有相同几何哈希、基组文件、哈密顿量、密度拟合设置、数值积分网格和收敛阈值时，才计算能量差。自旋处理按电荷态分别使用 RKS/UKS；如果有明确记录的 SOSCF 救援，其只改变迭代求解器。超时输出的最后一个 SCF 迭代值不构成能量标签。DFT 与 xTB 总能量没有交叉相减。

## 正负号与物理边界

固定几何 R 下定义 ΔE = E(阴离子; R) − E(中性; R)，电子附着能 A = −ΔE。A 为正只表示这个有限基组电子模型中的阴离子能量更低；A 为负表示带电态更高。不能据此宣称已校准的真实电子亲和能或还原电位。

SCF 收敛不等于电子解稳定。本次没有轨道 Hessian 稳定性分析，标准计算每个状态只用了一个 SAD 初猜；二阶 SCF 收敛算法也不能替代稳定性检查。成功计算中的 S² 是相应 Kohn–Sham 行列式诊断。阴离子正的占据轨道能是警告，但有限弥散基组可能人为局域不束缚的电子；即使轨道能为负也不能独自证明物理束缚态。

没有验证 def2-SVPD 基组、75×302 积分网格和固定几何的收敛性，也没有进行 DFT 优化/频率、非平衡溶剂化、电极电子库、参比标定、Cu 位点计算、真实动力学或实验。两种近似泛函一致仍可能共享系统误差，不能称为校准。

## 资源决策与下一步

每次实际启动使用 2 线程、Psi4 内部内存设置 512 MiB 和 300 s 超时；进程树 RSS 上限随主机可用内存收紧，不并行运行量子作业。启动前内存不足与已启动但超时分别存档。0.1 s 采样的 RSS 监控不等于操作系统硬内存配额，可能漏过极短峰值。

原始输入、输出、软件版本、轨道文件（成功时）、时间、采样内存和哈希均保留。实际证据支持小分子方法/资源诊断；尚不支持催化剂排名、还原电位、选择性或工业生产率。后续应单独预算 SCF/基组/阴离子束缚性验证，或在获得 SI/HPC 后转向有来源的电极模型。SCF 收敛和稳定性的区别可参见 [Psi4 官方文档](https://psicode.org/psi4manual/master/scf.html#stability-analysis)。

## 重现

运行 `python summarize_dft.py` 后执行 `python -m unittest discover -s . -p "test_dft_archive.py" -v`，不需要重跑 Psi4。运行 `python write_summary_reports.py` 重建本中英文说明。量子计算重跑命令见 README.md。大型可重建 scratch 留在本地，另有文件大小和哈希清单，不推送到 Git。
'''
(ROOT / "DFT_REPORT_EN.md").write_text(en, encoding="utf-8")
(ROOT / "DFT_REPORT_ZH.md").write_text(zh, encoding="utf-8")
print("Wrote DFT_REPORT_EN.md and DFT_REPORT_ZH.md from summary.json")
