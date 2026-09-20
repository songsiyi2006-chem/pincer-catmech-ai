# Neutral kinetics/transport verification kernel / 中性动力学—传质验证内核

**Status: software verification only. This is not the Cu-catalysed electroreductive coupling/ring-opening mechanism, and no fixture number is a prediction for that reaction.**

**状态：仅完成软件验证。此内核不是铜催化电还原偶联开环机理，所有测试参数均为人为构造的中性数学夹具，不得写入真实化学能量表、AI 训练标签或催化剂排名。**

## Scope / 范围

The package implements the neutral, single-site cycle

```text
A(solution) + * <=> A* <=> P* <=> P(solution) + *
```

All six Eyring rates are derived from the same state energies and three shared transition-state energies. Forward and reverse constants cannot be fitted independently through the energy API. Ideal solution activities are `c/c_standard`, with default `c_standard = 1000 mol/m³` (1 mol/L). Constants are in s⁻¹ because adsorption multiplies dimensionless activity. Matrix-tree weights give nonnegative stationary coverages; an independent linear solve verifies them.

正逆速率由同一组态能和过渡态导出；用无量纲活度保持标准态一致。三态稳态概率由矩阵树公式求解，并使用独立线性方程求解作对照。所有态必须具有一致的原子库存和能量参考；本程序无法仅凭数字核查这一化学前提。

The film calculation jointly solves

```text
J = km_A (c_A,bulk - c_A,surface)
  = km_P (c_P,surface - c_P,bulk)
  = Gamma * TOF(surface activities)
```

`Gamma` is mol of active sites per geometric m²; `J` is mol/(geometric m² s). The scalar root is bracketed by nonnegative surface concentrations. This is **effective stagnant-boundary-layer transport**, not a selective ionic membrane or Nernst–Planck simulation.

塞流模型将局部表面反应与边界层传质联立，采用保正、守恒的一阶隐式空间离散，并用独立 Radau 自适应积分验证。入口已经含有的产物不计入新生成产物；逆反应与负净生成保留。反应器假设等温、单液相、恒定体积流量、均匀可访问位点与有效面积、无轴向扩散。

For equal rate constants and equal film coefficients, an independent closed form exists. Let `Ct = c_A + c_P`, `lambda = Gamma*k/[3(c_standard+Ct)]`, `k_eff = lambda/(1+2*lambda/km)` and `area = volume*area_density`:

```text
c_A,out = Ct/2 + (c_A,in - Ct/2) * exp(-2*area*k_eff/flow)
```

This analytic limit tests kinetics and transport together; it is not an empirical result.

## Electrical accounting / 电耗核算

`cell_energy` integrates full-cell `Ucell(t)*I(t)` and divides by **isolated** product mass. It includes all sampled current, including capacitive current, in electrical consumption. Pumping, separation and auxiliary loads are explicitly excluded. Temporal sampling convergence is the caller's responsibility.

`faradaic_efficiency` separately uses **formed product moles before isolation**, measured/supplied total charge and externally validated electron stoichiometry. The neutral A→P model cannot provide electron stoichiometry or current. The independent electrical unit fixture must not be linked to its synthetic product flux.

电耗使用全槽电压和总电流，不用阴极过电位替代；分离收率影响每千克分离产品能耗，但不改变由分离前形成量定义的 FE。电子计量缺失时不能调用真实反应 FE 核算。超过 100% 的 FE 会明确报错。

## Reproduce / 复现

Requires Python ≥3.10, NumPy ≥1.24, SciPy ≥1.10. Existing environments should be reused. Actual run versions are recorded in the audit JSON; NumPy 1.x compatibility is implemented using the `trapz` fallback, while the recorded execution used NumPy 2.5.2.

From the project root (`projects/cu-np-operando`):

```powershell
python scripts/run_audit.py
```

Outputs:

- `results/neutral_transport_audit.json`: executed versions, source hashes, check counts, thermodynamic/site/film/source-outlet residuals, analytic and grid convergence values, electrical unit verification, and explicit prediction rejection.
- `results/neutral_transport_test_log.txt`: complete unittest output.

The final audit reruns 20 checks, including analytical limits, equilibrium, independent stationary and axial solvers, forward/reverse integrated source conservation, zero reaction, first-order grid convergence, units, overflow rejection, imported-rate consistency and evidence gates. Passing tests establish implementation behavior under the tested fixtures; they do not establish chemical validity.

## Recorded verification / 已执行验证

The first implementation run obtained 19/19 passing tests; the reviewed final version has 20/20 passing tests. With 20, 40, 80 and 160 axial cells, outlet absolute errors versus the analytic solution were approximately 0.2200, 0.1108, 0.05560 and 0.02785 mol/m³, consistent with first-order convergence. The 160-cell error normalized by 100 mol/m³ inlet concentration was 2.785×10⁻⁴; the independent Radau outlet error was approximately 5.40×10⁻¹³ mol/m³. See the regenerated JSON for authoritative values.

## Scientific stop / 科学结论门控

`target_ranking_gate` always rejects a Cu-electrosynthesis ranking from this neutral kernel, even if someone relabels the metadata as a real calculation. Missing work includes the actual multi-reactant, electron/proton-balanced network; validated constant-potential free energies and barriers; potential/current self-consistency; competing reactions; ionic transport/conduction; deactivation; uncertainty; and experiments. No physical catalyst ordering, industrial yield, FE advantage or economic claim is supported by this module.

此模块的成熟度仅为 A：可运行工程原型。把测试通过称为 G2 机理验证、G3 完整跨尺度验收或 G4 排名反转证据，均不符合本项目要求。
