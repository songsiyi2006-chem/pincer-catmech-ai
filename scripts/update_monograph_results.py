"""Insert frozen source-backed results into both scientific monographs."""
from pathlib import Path
import csv
import json
import re

REPO = Path(__file__).resolve().parents[1]


def replace_section(text, number, title, body):
    pattern = rf"## {number}\. .*?(?=\n## \d+\. |\Z)"
    result, count = re.subn(pattern, lambda _: f"## {number}. {title}\n\n{body.strip()}\n", text, flags=re.S)
    if count != 1:
        raise ValueError(f"Expected exactly one section {number}")
    return result


def add_figure(text, before_section, caption, stem):
    line = f"![{caption}](../examples/plots/{stem}.svg)"
    if line in text:
        return text
    return text.replace(f"\n## {before_section}.", f"\n{line}\n\n## {before_section}.", 1)


def main():
    status = json.loads((REPO / "data/CAMPAIGN_STATUS.json").read_text(encoding="utf-8"))
    with (REPO / "data/datasets/catalyst_descriptors.csv").open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    with (REPO / "data/datasets/thermochemistry/reaction_free_energies.csv").open(encoding="utf-8-sig", newline="") as stream:
        reactions = [row for row in csv.DictReader(stream) if float(row["temperature_K"]) == 383.15]
    reaction_map = {(row["catalyst_id"], row["reaction"]): float(row["delta_G_kcal_mol"]) for row in reactions}
    thermal = status["certified_thermochemistry"]
    n, active, states = status["alpb_conformer_jobs"], status["best_found_states"]["active"], status["best_found_state_total"]
    accepted_ts = status["accepted_dehydrogenation_TS_records"]
    elapsed = status["actual_main_window_seconds"] / 3600
    start, finish = status["campaign_started_utc"], status["main_finished_utc"]
    backbone = {"macho_pnp": "PNP", "pyridine_pnn": "pyridine PNN", "bipyridine_pnnoh": "bipyridine PNN(O)"}
    tables = {}
    for language in ("EN", "ZH"):
        table = ("| Metal | Backbone | R | Angle / deg | Vbur / % | ΔG‡ / kcal mol−1 |\n" if language == "EN" else
                 "| 金属 | 骨架 | R | 夹角／度 | Vbur／% | ΔG‡／kcal mol−1 |\n")
        table += "| --- | --- | --- | ---: | ---: | --- |\n"
        for row in rows:
            angle = f"{float(row['terminal_donor_angle_deg']):.3f}" if row.get("terminal_donor_angle_deg") else "—"
            volume = f"{float(row['buried_volume_percent']):.3f}" if row.get("buried_volume_percent") else "—"
            table += f"| {row['metal']} | {backbone[row['backbone']]} | {row['substituent']} | {angle} | {volume} | — |\n"
        tables[language] = table
    abstracts = {
        "EN": f"""**Frozen bounded-campaign record.** Real local quantum computation began at {start} and the main search ended at {finish}, spanning {elapsed:.4f} hours including initial gas-phase preparation and the subsequent ALPB-toluene campaign. This is a report of the actual calculations, including unsuccessful searches. The requested complete transition-state, kinetic, and learned-prediction outcomes were not achieved within this run.

Borrowing-hydrogen N-alkylation of aniline with benzyl alcohol was investigated through a factorial library of 24 scaffold-informed designs: Ru(II), Mn(I), Fe(II), and Co(I), three proton-responsive ligand backbones, and Ph/iPr phosphine substitution. The solution-model conformer controller completed {n} native optimization jobs. Its full chemical-identity and donor-coordination screens retained best-found geometries for {states} of 72 intended catalyst/state combinations, including {active} active designs. Both Co bipyridine active designs rearranged through ligand oxygen–carbonyl coupling and are excluded from the intended active-state set. These designed analogues are not a library of experimentally established catalysts, and sampled electronic minima do not establish global minima or electronic ground spins.

Independent full-Hessian auditing selected {thermal['selected_catalyst_states']} catalyst-state minimum certificates and {thermal['selected_reference_species']} molecular reference species. Their consistent 13-temperature records contain {thermal['catalyst_temperature_rows']} catalyst rows and {thermal['reference_temperature_rows']} reference rows; balanced dehydrogenation, imine hydrogenation and separated-ion activation reactions yield {thermal['reaction_temperature_rows']} reaction-temperature rows. A separately certified K–OtBu contact pair and an N-deprotonated hemiaminal cluster provide additional explicitly scoped molecular evidence. Reaction free energies are thermodynamic differences, not activation barriers.

The user-selected conditions are ALPB toluene, 383.15 K baseline, potassium tert-butoxide at 0.05 equivalent, and a requested total-base scan of 0.01–0.20 equivalent. Actual paths were attempted for two Ru, two Mn, and two Co representatives and for explicit condensation clusters. At freeze, {accepted_ts} dehydrogenation transition-state records passed the complete configured acceptance gates. No complete six-step source-verified reaction network was available, so no real TOF grid, steady-state Campbell degree of rate control, trained surrogate, or Pareto catalyst ranking is claimed. The deliverable includes executable evidence gates, real failed and successful numerical records, endpoint thermodynamics, and an explicit list of the scientific targets still unmet.""",
        "ZH": f"""**已冻结的限时计算记录。** 本地真实量子计算从 {start} 开始，主构象搜索于 {finish} 结束，覆盖 {elapsed:.4f} 小时，包含初始气相预优化及随后按用户条件开展的 ALPB 甲苯计算。本报告陈述实际得到的数据及失败结果；本次运行未在时限内完成全部过渡态、微观动力学与训练模型目标。

研究对象为苯胺与苄醇的借氢 N-烷基化。四种金属分支 Ru(II)、Mn(I)、Fe(II)、Co(I)，三种质子响应骨架与两种膦取代基组成二十四项设计。溶液模型控制器实际完成 {n} 次原生优化任务，经完整化学身份与供体配位检查，在七十二个催化剂—状态组合中获得 {states} 个合格的最佳已发现几何，其中包括 {active} 个活性设计。两种 Co 联吡啶活性模型均发生配体氧与羰基连接重排，因此不计入目标活性态集合。这里的结构属于按文献骨架启发的设计类似物，不能表述为二十四种已在实验中确证的催化剂；有限构象搜索也不证明全局最低点或电子基态自旋。

独立的完整 Hessian 审核最终选出 {thermal['selected_catalyst_states']} 个催化剂状态证书和 {thermal['selected_reference_species']} 种分子参照。每个证书使用自己的几何与完整振动谱，覆盖十三个温度，共导出 {thermal['catalyst_temperature_rows']} 条催化剂热化学温度行、{thermal['reference_temperature_rows']} 条参照温度行，以及 {thermal['reaction_temperature_rows']} 条满足原子与电荷守恒的反应温度行。另有单独认证的 K–OtBu 接触离子对和去 N-质子化半缩醛胺簇模型。反应自由能只说明指定分子参照之间的热力学差异，不能替代活化自由能。

条件为 ALPB 甲苯、383.15 K 基准温度、叔丁醇钾基准零点零五当量及零点零一至零点二零当量的请求扫描域。两 Ru、两 Mn、两 Co 代表以及显式缩合簇均接受真实路径探索；冻结时有 {accepted_ts} 条脱氢过渡态记录通过全部规定门槛。由于缺乏完整六步、来源一致的反应网络，实际 TOF 网格、稳态 Campbell 速率控制度、训练代理模型和 Pareto 催化剂排名均未得到。本交付保留数值上成功与失败的记录、端点热力学、可执行证据门槛，以及尚未达到的科学目标，不能把软件测试通过解读为催化机理或实验活性已被验证。"""}
    results = {
        "EN": f"""The following frozen table uses only ALPB-toluene best-found active geometries from `data/datasets/catalyst_descriptors.csv`. There are {active} identity-preserving entries among 24 designs. Angles denote the actual terminal donors: P–M–P for PNP, P–M–N for both PNN backbones. Buried-volume values use 50,000 Monte Carlo samples; per-geometry Wilson intervals and full-precision coordinates are retained in the CSV and `data/structures/`. A dash means unavailable, never zero. Every activation-free-energy cell remains unavailable because no corresponding source-verified free-energy barrier dataset has been produced.

{tables['EN']}
The two missing Co bipyridine active entries have a specific computational explanation: an unintended ligand O–carbonyl C connection persists across alternative starting geometries. The Ph case includes an O–C distance of 1.4963 Å and native xTB Wiberg bond order 0.8337 in the saved diagnostic structure. This is evidence of rearrangement on this approximate potential, not proof of a solution reaction yield or experimental decomposition rate. Those structures may possess positive Hessians as mathematical minima while still failing the intended catalyst identity.

Electronic-energy selection and thermochemical certification serve different purposes. The descriptor table reports the lowest sampled eligible electronic geometry, while the thermochemistry CSV selects the lowest certified G at 383.15 K and retains that same certified geometry at all thirteen temperatures. These choices are explicitly linked by geometry hashes and are not merged as if they were necessarily the same conformer. Unexpected Co PNN angles and several iPr/Ph volume reversals therefore remain geometrical observations whose chemical explanation requires additional calculations.

The append-only search records distinguish assembly failures, electronic nonconvergence, donor loss, covalent rearrangement, Hessian rejection, and finite-budget path failure. Current-round controller counters differ from cumulative log counts after a resumed run; the frozen status uses the full historical record. Native evidence archives retain the actual optimization and gradient output, including failed calculations. `data/CAMPAIGN_STATUS.json` and the source-verified thermal export record the final counts without erasing the early gas-phase preparation history.

![Active-state buried volume; hatched cells failed intended identity](../examples/plots/steric_buried_volume_map.svg)""",
        "ZH": f"""下表完全使用 `data/datasets/catalyst_descriptors.csv` 冻结的 ALPB 甲苯活性态最佳已发现几何。二十四项设计中有 {active} 项保留目标身份。PNP 角度为 P–M–P，两类 PNN 角度为 P–M–N；埋藏体积采用五万个 Monte Carlo 点，逐项 Wilson 区间及完整精度坐标见 CSV 与 `data/structures/`。破折号表示没有可用结果，绝不等于零。由于尚未形成来源验证完整的活化自由能数据集，本表全部势垒单元格仍为空缺。

{tables['ZH']}
两种缺失的 Co 联吡啶活性设计有明确的计算原因：更换起始几何后，非预期配体 O—羰基 C 连接仍反复出现。苯基案例所保存诊断结构的 O—C 距离为 1.4963 Å，原生 xTB Wiberg 键级为 0.8337。它支持该近似势能面上发生重排的解释，却不能确定溶液中的副反应产率或实验分解速率。重排结构即便具有全正 Hessian，也只是另一种化学身份的数学极小点，不能计入原定活性态。

电子能选择与热化学认证具有不同任务。描述符表使用抽样中电子能最低的合格几何；热化学表则在已获得完整证书的构象中，按 383.15 K 的 G 选择，并在十三个温度统一沿用该证书几何。二者通过文件路径、原子顺序和几何哈希明确绑定，不假定必然来自同一个构象。Co PNN 的特殊夹角以及若干异丙基与苯基体积次序反转，只属于结构观测，仍需进一步计算解释。

追加式记录分别保留装配错误、电子不收敛、供体脱落、共价重排、Hessian 拒绝与路径预算耗尽。恢复运行后控制器当前轮次计数可能不同于累计日志，最终统计按完整记录汇总。原生证据归档包含真实优化、梯度及失败输出；`data/CAMPAIGN_STATUS.json` 和热化学导出清单给出冻结数目，同时保留早期气相准备的历史。

![活性态埋藏体积；斜线格表示目标身份未保留](../examples/plots/steric_buried_volume_map.svg)"""}
    for language in ("EN", "ZH"):
        path = REPO / "docs" / f"MONOGRAPH_BORROWING_HYDROGEN_{language}.md"
        text = path.read_text(encoding="utf-8")
        text = replace_section(text, 1, "Abstract and frozen campaign status" if language == "EN" else "摘要与冻结计算状态", abstracts[language])
        text = replace_section(text, 9, "Frozen structural results and acceptance boundaries" if language == "EN" else "冻结结构结果与接受边界", results[language])
        text = text.replace("at this interim snapshot", "at the frozen evidence level").replace("present interim evidence level", "present frozen evidence level")
        text = text.replace("The early calculations demonstrate", "The completed bounded calculations demonstrate")
        text = text.replace("本阶段性快照", "本次冻结记录").replace("当前阶段性证据", "当前冻结证据")
        text = text.replace("早期真实计算表明", "已完成的限时真实计算表明")
        text = text.replace("后续必须使用新的 ALPB 能量与力，统一参考态", "本次已经重新计算 ALPB 能量与力并统一参考态")
        text = text.replace("本阶段记录", "本次冻结记录")
        text = text.replace("At this reporting snapshot, dimer continuation initialized from the diagnosed mode remained a separate search in progress. Its eventual outcome must be taken from its own force, Hessian, identity, and downhill-connectivity records; the present diagnostic establishes neither success nor impossibility of the Mn pathway.",
                            "By the final freeze, dimer continuation initialized from this diagnosed mode had ended without TS certification. Its final full-Hessian diagnostic appears in the table below, alongside the original candidate; each retains its own force, mode and identity records. These unsuccessful bounded searches do not establish that the Mn pathway is impossible.")
        text = text.replace("在本段记录快照中，以该模式初始化的二聚体延续仍是进行中的独立搜索，其最终状态必须由自身的力、Hessian、身份和双向下坡连接证据决定。这份诊断既不宣告后续成功，也不能证明 Mn 路径不存在。",
                            "至最终冻结时，以该模式初始化的二聚体续算已结束，仍未获得过渡态认证。续算末帧的完整 Hessian 诊断与原候选一起列于下表，各自保留独立的真实力、振动模式和身份记录。有限时限内的搜索失败不能证明 Mn 路径不存在。")
        thermal_title = "Certified catalyst reaction thermodynamics" if language == "EN" else "认证的催化剂反应热力学"
        reaction_table = ("| Design | Alcohol dehydrogenation ΔG | Imine hydrogenation ΔG | Separated-ion deprotonation ΔG |\n" if language == "EN" else
                          "| 设计 | 醇脱氢反应 ΔG | 亚胺加氢反应 ΔG | 分离离子去质子化 ΔG |\n")
        reaction_table += "| --- | ---: | ---: | ---: |\n"
        for row in rows:
            identifier = row["catalyst_id"]
            if (identifier, "alcohol_dehydrogenation") in reaction_map:
                label = f"{row['metal']} / {backbone[row['backbone']]} / {row['substituent']}"
                reaction_table += f"| {label} | {reaction_map[(identifier, 'alcohol_dehydrogenation')]:+.3f} | {reaction_map[(identifier, 'imine_hydrogenation')]:+.3f} | {reaction_map[(identifier, 'separated_ion_deprotonation')]:+.3f} |\n"
        discussion = ("All table entries are kcal/mol at 383.15 K, for separate 1 M molecular species and the fixed ALPB potential. Alcohol dehydrogenation denotes Cat + benzyl alcohol → Cat–H₂ + benzaldehyde; activation denotes protonated Cat(+) + tBuO(−) → Cat + tBuOH. Each sum is atom- and charge-balanced and uses the complete certificate attached to its own geometry.\n\n"
                      "Both Ru PNP alcohol-dehydrogenation reactions are exergonic in this model, whereas both Mn PNP alcohol-dehydrogenation reactions are endergonic. The imine-hydrogenation endpoint reactions show the opposite sign pattern for these same four designs. This supports a specific thermodynamic distinction between these chosen molecular states. It does not establish a lower Ru barrier or faster Ru turnover: the association basins, true saddles, imine reduction, condensation and catalyst populations have not been fully determined. Several pyridine-PNN reactions are more favorable than their PNP counterparts; replacing the responsive atom and ligand topology also changes the molecular state, so this cannot be assigned solely to metal or substituent electronics.\n\n"
                      "The strongly negative separated-ion proton-transfer values describe neutralization of deliberately isolated charged references in a low-dielectric model. They must not be interpreted as measured pKa values, a dissolved tBuOK activation equilibrium, or evidence that 0.05 equivalent total salt generates a specified free-anion concentration. Contact pairing and the counterion chemical potential are absent from that reaction definition. The complete 13-temperature series is available for thermal-model sensitivity, while the omitted solvent temperature derivatives remain a distinct limitation.\n\n" if language == "EN" else
                      "本表统一为 383.15 K、分离的一摩尔每升分子物种以及固定 ALPB 势能模型下的 kcal/mol。醇脱氢反应指 Cat＋苄醇→Cat–H₂＋苯甲醛；活化参照指质子化 Cat(+)＋叔丁醇盐(−)→Cat＋叔丁醇。所有加和均满足原子与电荷守恒，并使用各自几何绑定的完整证书。\n\n"
                      "两个 Ru PNP 设计的醇脱氢反应在本模型中均为放能，两个 Mn PNP 设计的醇脱氢反应均为吸能；这四个设计的亚胺加氢端点反应则呈相反的正负号分布。这只支持所选分子状态之间的热力学差别，不能推出 Ru 势垒更低或周转更快；缔合盆地、真实鞍点、亚胺还原、缩合及催化剂分布尚未全部确定。若干吡啶 PNN 反应比 PNP 更有利，但响应原子和骨架拓扑也同时改变，因此不能将其全部归因于金属或取代基电子效应。\n\n"
                      "较大的负去质子化自由能对应低介电模型中刻意分开的带电参照发生中和。它们不是测量的 pKa、真实溶解叔丁醇钾的活化平衡，也不能证明零点零五当量总盐对应某个自由阴离子浓度；该反应定义没有包含接触离子配对与反离子化学势。十三温度完整序列可用于既定热模型内的敏感性研究，尚未计算的溶剂温度导数仍属于另一项模型局限。\n\n")
        closure_note = ("Imine hydrogenation uses Cat–H₂ + imine → Cat + product amine. Adding this reaction to alcohol dehydrogenation and the two certified condensation reference reactions cancels the catalyst states and intermediates exactly. The independently exported cycle residual checks this thermodynamic bookkeeping at every temperature. Favorable hydrogen uptake can be accompanied by less favorable hydrogen delivery; neither isolated endpoint difference establishes a rate or selects a catalyst.\n\n" if language == "EN" else
                        "亚胺加氢定义为 Cat–H₂＋亚胺→Cat＋产物胺。它与醇脱氢及两项已认证缩合参照反应相加时，催化剂状态和中间体必须准确相消。独立导出的循环残差在每个温度检查该热力学闭合关系。较有利的储氢可以伴随较不利的放氢；任一孤立端点差均不能确定速率或选择最优催化剂。\n\n")
        addition = f"### {thermal_title}\n\n{discussion}{closure_note}{reaction_table}\n"
        if f"### {thermal_title}" in text:
            text = re.sub(rf"### {re.escape(thermal_title)}.*?(?=\n## 13\.)", lambda _: addition, text, flags=re.S)
        else:
            text = text.replace("\n## 13.", f"\n{addition}\n## 13.", 1)
        text = add_figure(text, 13, "Certified endpoint free energies; no transition-state barrier" if language == "EN" else "认证端点自由能；不是过渡态势垒", "pes_comparison_ru_vs_mn")
        diagnostic_title = "Actual nonstationary curvature diagnostics" if language == "EN" else "真实非驻点曲率诊断"
        diagnostic_rows = []
        for diagnostic_path in sorted((REPO / "data/campaign/neb").glob("*/diagnostic.json")):
            diagnostic = json.loads(diagnostic_path.read_text(encoding="utf-8"))
            if diagnostic.get("status") != "diagnostic_complete":
                continue
            modes = diagnostic["negative_mode_transfer_coordinates"]
            mode = min(modes, key=lambda row: row["frequency_cm1"]) if modes else None
            label = diagnostic["catalyst_id"].replace("_macho_pnp_", " PNP ")
            label += " / " + diagnostic["finished_utc"][11:19] + " UTC"
            cells = [label, f"{diagnostic['true_fmax_eV_A']:.4f}", str(diagnostic['negative_mode_count'])]
            cells += [f"{mode[key]:.3f}" if mode else "—" for key in ("frequency_cm1", "proton_overlap", "hydride_overlap")]
            diagnostic_rows.append("| " + " | ".join(cells) + " |")
        diagnostic_text = ("These fixed-candidate full-Hessian calculations are diagnostics at nonstationary geometries. Columns report unmodified maximum force (eV/Å), negative-mode count, and the most negative mode's frequency (cm⁻¹) and separate proton/hydride overlaps. They are not stationary-point certificates or activation free energies. Later mode-guided searches retain separate results in the final NEB inventory.\n\n| Candidate | Force | Negative modes | Lowest frequency | Proton overlap | Hydride overlap |\n" if language == "EN" else
                           "以下为固定候选几何的完整 Hessian 非驻点诊断。依次列出未经修改的最大力（eV/Å）、负模式数、最负模式频率（cm⁻¹）及其质子／氢负离子重叠。它们不是驻点证书或活化自由能；随后以模式为依据的搜索另有结果文件，见最终 NEB 清单。\n\n| 候选 | 最大力 | 负模式数 | 最负频率 | 质子重叠 | 氢负离子重叠 |\n")
        diagnostic_text += "| --- | ---: | ---: | ---: | ---: | ---: |\n" + "\n".join(diagnostic_rows)
        diagnostic_addition = f"### {diagnostic_title}\n\n{diagnostic_text}\n"
        if f"### {diagnostic_title}" in text:
            text = re.sub(rf"### {re.escape(diagnostic_title)}.*?(?=\n## 14\.)", lambda _: diagnostic_addition, text, flags=re.S)
        else:
            text = text.replace("\n## 14.", f"\n{diagnostic_addition}\n## 14.", 1)
        text = add_figure(text, 18, "Kinetic grid availability; no numerical TOF result" if language == "EN" else "动力学网格可用性；未产生数值 TOF", "microkinetic_tof_heatmap")
        path.write_text(text, encoding="utf-8", newline="\n")
        assert len(re.findall(r"^## \d+\.", text, re.M)) == 20
    print(json.dumps({"updated": ["EN", "ZH"], "source_frozen_utc": status["frozen_at_utc"]}))


if __name__ == "__main__":
    main()
