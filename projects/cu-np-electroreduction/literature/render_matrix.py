"""Regenerate the human-readable audit from the curated JSON, using only stdlib."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def main():
    data = json.loads((HERE / 'evidence_matrix.json').read_text(encoding='utf-8'))
    rows = data['sources']
    if len(rows) != 12 or len({r['doi'].lower() for r in rows}) != 12:
        raise ValueError('Expected twelve unique original-study DOI records.')
    required = ['system', 'central_claim', 'interface_treatment', 'evidence',
                'limitations', 'project_overlap', 'potential_novel_contribution']
    for row in rows:
        if not all(row.get(k, {}).get('text') for k in required):
            raise ValueError(f"Incomplete matrix entry: {row['id']}")
        if not row['author_roles'].get('source'):
            raise ValueError(f"Missing author-role provenance: {row['id']}")
    out = [
        '# Literature evidence and novelty audit / 文献证据与创新边界审计', '',
        f"Audit date / 检索截止：{data['audit_date']}。", '',
        '`LITERATURE FACT` records what a source reports; it does not mean independent reproduction. '
        '`MODEL INFERENCE` denotes this audit’s interpretation. `WORKING HYPOTHESIS` denotes an unproven contribution. '
        '/ 文献事实不等于独立复现；推断与待检验假说均单独标示。', '',
        '## Count and relevance / 数量与相关性', '',
        '| Class | Count | Entries |', '|---|---:|---|',
        '| Primary Cu–N–P/quinazolinone system | 1 | A0 |',
        '| Other organic-electroreduction comparators | 3 | N1, N5, N7 |',
        '| Interface/speciation/method neighbors | 4 | N2, N3, N4, N6 |',
        '| Mandatory conceptual anchors | 4 | A1–A4 |', '',
        'The matrix contains 12 original studies. Eight are relevant under an explicit '
        'primary-plus-organic-plus-methodological definition; four conceptual anchors are excluded from that count. '
        'Only four directly study organic electroreduction, a shortfall of four if the eight-study threshold '
        'is interpreted to require that narrower scope. One study concerns the exact primary reaction; '
        'no additional same-system original study was verified. The contract does not require eight studies on an identical reaction.', '',
        '共12项原创研究；仅在明确包含界面、位点和方法近邻的定义下，核心及近邻合计8项。'
        '直接研究有机电还原的只有4项；若将“8项直接相关”限定于有机电还原，则尚缺4项。'
        '同一Cu–N–P/喹唑啉酮体系仅1项，没有将热催化或CO₂研究称作同反应研究。', '',
        '## Evidence matrix / 证据矩阵', '',
        '| ID / source | System/electrolyte: LITERATURE FACT | Claim/evidence: LITERATURE FACT | '
        'Interface: LITERATURE FACT | Limits/overlap: MODEL INFERENCE | Possible contribution: WORKING HYPOTHESIS |',
        '|---|---|---|---|---|---|',
    ]
    esc = lambda text: text.replace('|', r'\|').replace('\n', ' ')
    for r in rows:
        cells = [f"{r['id']} [{r['doi']}]({r['publisher_url']})", r['system']['text'],
                 r['central_claim']['text'] + ' ' + r['evidence']['text'],
                 r['interface_treatment']['text'],
                 r['limitations']['text'] + ' ' + r['project_overlap']['text'],
                 r['potential_novel_contribution']['text']]
        out.append('| ' + ' | '.join(map(esc, cells)) + ' |')
    out += ['', '## Source access, locators and author roles / 访问范围、证据定位与作者角色', '']
    for r in rows:
        ar = r['author_roles']
        out += [f"### {r['id']}: {r['title']}", '',
                f"- Classification: `{r['counting_class']}`.",
                '- Article: ' + r['access_status']['article'],
                '- Supporting information: ' + r['access_status']['si'],
                '- Evidence locations: ' + '; '.join(r['evidence_locators']),
                '- Sources: ' + '; '.join(f'[{i+1}]({s})' for i, s in enumerate(r['sources'])),
                '- Authors in order: ' + '; '.join(ar['authors']) + '.',
                '- Corresponding authors: ' + ('; '.join(ar['corresponding'])
                  if ar.get('corresponding') else 'NOT VERIFIED in accessible source') + '.',
                '- Role verification: ' + ar.get('contribution_scope', ar.get('gxnu_roles', '')),
                '- Role source: [original source](' + ar['source'] + '). ' + ar['verification_limit']]
        if ar.get('equal_contribution'):
            out.append('- Equal contribution explicitly reported: ' + '; '.join(ar['equal_contribution']) + '.')
        out.append('')
    out += [
        '## Novelty boundary / 创新边界', '',
        '1. **MODEL INFERENCE:** N1 establishes a Cu/P organic-electroreduction precedent. '
        'Adding phosphorus or redrawing a static adsorption diagram is insufficient novelty.',
        '2. **MODEL INFERENCE:** N2–N4 and N6 already address dynamic Cu sites, proton-associated '
        'instability or clusters. Their existence motivates controls but does not prove that these states occur in DMF/TEOA.',
        '3. **MODEL INFERENCE:** N2 and N6 report different reconstruction reversibility. '
        'Support, electrolyte and potential protocol must be examined before transferring either result.',
        '4. **MODEL INFERENCE:** N5’s dechlorination ratio and FE are different observables. '
        'Its fixed added-charge description does not independently establish constant-potential feedback.',
        '5. **WORKING HYPOTHESIS:** A useful target-specific contribution would distinguish '
        'coordination-only adsorption from changing site populations or interfacial transfer and predict a '
        'condition where the Cu–N–P advantage weakens. No validated catalyst ranking is supplied by this audit.', '',
        'Reproduction, application to a new system, mechanistic discovery, predictive methodology and '
        'prospective experimental validation are separate achievements. DFT labels do not supply independent experimental validation.', '',
        '## Retrieval and search / 检索与原文追溯', '',
        'See [search_log.json](search_log.json) for query scope and the bounded follow-up search. '
        'The 2026 original N7 predates the audit cutoff; it is not a same-system follow-up. '
        'Reviews and unrelated citing works are not counted. Absence from a bounded search is not proof of absence.', '',
        'See [README.md](README.md) for access levels, source manifests and local-only third-party '
        'files. No inaccessible SI parameters have been reconstructed from abstracts.', '',
    ]
    (HERE / 'evidence_matrix.md').write_text('\n'.join(out), encoding='utf-8')
    print('Rendered 12 unique source records with explicit relevance classes and role provenance.')


if __name__ == '__main__':
    main()
