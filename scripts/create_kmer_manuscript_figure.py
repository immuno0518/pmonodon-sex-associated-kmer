from release_config import *
"""Render a local manuscript figure from completed, fixed k31 results."""
import csv
import gzip
import hashlib
import json
import os
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

ROOT = RUN_ROOT
os.environ['MPLCONFIGDIR'] = str(ROOT / 'work/matplotlib_cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

BASE = ROOT / 'outputs/kmer'
OUT = BASE / 'manuscript_main_figure'
COLORS = {'F': '#B53C7D', 'M': '#2375AD'}


def read_tsv(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(p, rows):
    with p.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    chrom_path = BASE / 'reference_mappability/unique_opportunity_normalized_chromosomes.tsv'
    mapping_path = BASE / 'heldout_mapping/all_replicated_candidates_reference_mapping.tsv.gz'
    individual_path = BASE / 'regional_sensitivity/per_individual_region_threshold_sensitivity.tsv'
    chroms = read_tsv(chrom_path)
    assert sum(int(r['unique_k31_starts']) for r in chroms) == 1139681763
    assert [sum(int(r[g + '_unique']) for r in chroms) for g in 'FM'] == [180082, 184371]

    # The last display group pools all unplaced contigs except the separately shown Contig6.
    panel_a = [dict(r) for r in chroms if r['chromosome'].startswith('Chr') or r['chromosome'] == 'Contig6']
    others = [r for r in chroms if r['chromosome'].startswith('Contig') and r['chromosome'] != 'Contig6']
    pooled = {'chromosome': 'Other_contigs', 'unique_k31_starts': sum(int(r['unique_k31_starts']) for r in others)}
    for g in 'FM':
        pooled[g + '_unique'] = sum(int(r[g + '_unique']) for r in others)
    panel_a.append(pooled)
    panel_a = [{k: r[k] for k in ('chromosome', 'unique_k31_starts', 'F_unique', 'M_unique')} for r in panel_a]
    for r in panel_a:
        for g in 'FM':
            r[g + '_per_1000_unique_starts'] = int(r[g + '_unique']) / int(r['unique_k31_starts']) * 1000
    assert len(panel_a) == 46
    assert sum(int(r['unique_k31_starts']) for r in panel_a) == 1139681763

    # Recover the verified concatenated-reference offset from the original manifests.
    names = identity_names()
    seqs = json.loads((BASE / 'reference_mappability/server_download/reference_sequence_manifest.json').read_text())['sequences']
    offsets = {}
    offset = 0
    for r in seqs:
        local = names[r['sequence']]
        assert r['length_bp'] == int(local['length_bp']) and r['sequence_sha256'] == local['sequence_sha256']
        offsets[local['original_id']] = offset
        offset += r['length_bp']
    freq = np.memmap(str(BASE / 'reference_mappability/reference_k31_e0.freq8'), mode='r', dtype=np.uint8)
    assert len(freq) == offset
    bins = []
    for a in range(2500001, 3100001, 10000):
        b = a + 9999
        n = int(np.count_nonzero(freq[offsets['Chr19'] + a - 1:offsets['Chr19'] + b] == 1))
        bins.append({'start_1based': a, 'end_1based': b, 'unique_k31_starts': n, 'F_unique': 0, 'M_unique': 0})
    target = Counter()
    scanned = 0
    with gzip.open(mapping_path, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            scanned += 1
            if r['exact_reference_hits'] != '1' or r['unique_chr'] != 'Chr19':
                continue
            a = int(r['unique_start_1based'])
            if 2500001 <= a <= 3100000:
                bins[(a - 2500001) // 10000][r['candidate_group'] + '_unique'] += 1
            if 2710000 <= a <= 2890000:
                target[r['candidate_group']] += 1
    assert scanned == 2566199 and dict(target) == {'F': 32753, 'M': 4507}
    target_opp = int(np.count_nonzero(freq[offsets['Chr19'] + 2710000 - 1:offsets['Chr19'] + 2890000] == 1))
    assert target_opp == 86367
    for r in bins:
        for g in 'FM':
            r[g + '_per_1000_unique_starts'] = r[g + '_unique'] / r['unique_k31_starts'] * 1000 if r['unique_k31_starts'] else ''

    raw = [r for r in read_tsv(individual_path) if r['region'] == 'Chr19_target' and r['read_count_threshold'] == '2']
    panel_cd = []
    summary = json.loads((BASE / 'regional_sensitivity/regional_threshold_summary.json').read_text())['summary']
    for g in 'FM':
        for sex in 'FM':
            rs = sorted([r for r in raw if r['candidate_group'] == g and r['sex'] == sex], key=lambda r: r['sample_id'])
            assert len(rs) == 92
            s = next(r for r in summary if r['region'] == 'Chr19_target' and r['candidate_group'] == g and r['read_count_threshold'] == 2)
            assert abs(np.median([float(r['detected_fraction']) for r in rs]) - s[sex]['median']) < 1e-12
            for i, r in enumerate(rs):
                panel_cd.append({'candidate_group': g, 'sex': sex, 'individual_within_sex': i + 1,
                                 'panel_kmers': int(r['panel_kmers']), 'detected_kmers': int(r['detected_kmers']),
                                 'detected_percent': float(r['detected_fraction']) * 100})
    write_tsv(OUT / 'panel_A_chromosome_density.tsv', panel_a)
    write_tsv(OUT / 'panel_B_Chr19_10kb_density.tsv', bins)
    write_tsv(OUT / 'panels_CD_heldout_detection.tsv', panel_cd)

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.labelsize': 8,
                         'axes.titlesize': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .7, 'xtick.major.width': .6, 'ytick.major.width': .6,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(7.48, 6.9))
    gs = fig.add_gridspec(3, 2, left=.115, right=.975, bottom=.085, top=.945,
                          height_ratios=[1, 1.12, 1], hspace=.72, wspace=.48)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])

    # Both directions share the same physical opportunity denominator, including zeros.
    xpos = np.array(list(range(44)) + [45, 48])
    for g, dx, marker in [('F', -.14, 'o'), ('M', .14, '^')]:
        yy = [r[g + '_per_1000_unique_starts'] for r in panel_a]
        ax_a.scatter(xpos + dx, yy, color=COLORS[g], s=14, marker=marker, linewidths=0,
                     label=('Female' if g == 'F' else 'Male') + '-associated')
    for i in [18, 33, 44]:
        ax_a.axvspan(xpos[i] - .43, xpos[i] + .43, color='#F0E2B4', alpha=.65, zorder=0)
    ax_a.set_yscale('symlog', linthresh=.03, linscale=.6)
    ax_a.set_ylim(-.006, 65)
    ax_a.set_yticks([0, .03, .1, 1, 10])
    ax_a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: format(v, 'g')))
    ax_a.set_xlim(-.8, 49)
    ax_a.set_xticks(xpos)
    ax_a.set_xticklabels([str(i) for i in range(1, 45)] + ['Ct6', 'Other'], fontsize=6)
    ax_a.tick_params(axis='x', length=2, pad=3)
    ax_a.set_xlabel('1-44: chromosomes; Ct6: Contig6; Other: remaining unplaced contigs pooled', labelpad=5)
    ax_a.set_ylabel('Associated starts /\n1,000 unique starts')
    ax_a.set_title('Genome-wide distribution', loc='left', pad=8)
    ax_a.legend(frameon=False, loc='upper left', bbox_to_anchor=(.0, 1.0), fontsize=7, ncol=2,
                handletextpad=.3, columnspacing=1)
    ax_a.grid(axis='y', linewidth=.35, alpha=.25)

    ax_b.axvspan(2.71, 2.89, color='#F0E2B4', alpha=.7, zorder=0)
    xx = [(r['start_1based'] + r['end_1based']) / 2 / 1e6 for r in bins]
    for g in 'FM':
        yy = [r[g + '_per_1000_unique_starts'] if r[g + '_per_1000_unique_starts'] != '' else np.nan for r in bins]
        ax_b.plot(xx, yy, color=COLORS[g], lw=1.1, marker='o', markersize=2)
        low = [r['unique_k31_starts'] < 1000 for r in bins]
        # A 10-kb display follows the same 10% opportunity criterion as the original 100-kb plot.
        for x, y, flag in zip(xx, yy, low):
            if flag:
                ax_b.scatter([x], [y], facecolors='white', edgecolors=COLORS[g], s=17, linewidths=.6, zorder=5)
    ax_b.set_xlim(2.5, 3.1)
    ax_b.set_ylim(bottom=0)
    ax_b.set_xticks(np.arange(2.5, 3.1001, .1))
    ax_b.set_xlabel('Chr19 position (Mb; 10-kb bins)')
    ax_b.set_ylabel('Associated starts /\n1,000 unique starts')
    ax_b.set_title('Overlap with the Chr19 candidate interval', loc='left', pad=8)
    ax_b.text(.98, .92, 'Chr19: 2.71-2.89 Mb\n32,753 F-associated; 4,507 M-associated',
              transform=ax_b.transAxes, va='top', ha='right', fontsize=7,
              bbox={'facecolor': 'white', 'alpha': .85, 'edgecolor': 'none', 'pad': 2})
    ax_b.grid(axis='y', linewidth=.35, alpha=.25)

    rng = np.random.RandomState(20261021)
    for ax, g in [(ax_c, 'F'), (ax_d, 'M')]:
        for i, sex in enumerate('FM'):
            vals = np.array([r['detected_percent'] for r in panel_cd if r['candidate_group'] == g and r['sex'] == sex])
            ax.scatter(i + rng.uniform(-.16, .16, len(vals)), vals, s=10, alpha=.5,
                       color=COLORS[sex], edgecolors='none', rasterized=False)
            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            ax.plot([i, i], [q1, q3], color='#20282E', lw=2.4, zorder=5)
            ax.plot([i - .13, i + .13], [med, med], color='#20282E', lw=1.3, zorder=6)
            ax.text(i, 108, '{:.2f}%'.format(med), ha='center', va='center', fontsize=7, fontweight='bold')
        n = 32753 if g == 'F' else 4507
        ax.set_title(('Female' if g == 'F' else 'Male') + '-associated panel\n' + format(n, ',') + ' Chr19 interval k-mers', loc='left', pad=8)
        ax.set_xlim(-.48, 1.48)
        ax.set_ylim(-4, 116)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Female\n(n = 92)', 'Male\n(n = 92)'])
        ax.set_ylabel('Panel detected per individual (%)')
        ax.grid(axis='y', linewidth=.35, alpha=.25)
    for letter, ax in zip('ABCD', [ax_a, ax_b, ax_c, ax_d]):
        ax.text(-.105 if ax in [ax_a, ax_b] else -.27, 1.08, letter, transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='bottom')
    fig.text(.115, .017, 'One full-sib family: 16 discovery + 184 heldout individuals. Spatial densities and selected-panel detection are descriptive.',
             fontsize=6.5, color='#4B555C')
    for fmt in ('pdf', 'svg', 'png'):
        fig.savefig(OUT / ('kmer_main_figure.' + fmt), dpi=600 if fmt == 'png' else 300, facecolor='white')
    plt.close(fig)

    caption = '''Figure X. Read-derived k31 associations converge on the Chr19 candidate interval.
(A) Chromosome-level densities of direction-consistent, Bonferroni-significant female-associated (magenta circles) and male-associated (blue triangles) k-mers with exactly one whole-reference match. Counts of unique candidate start positions are divided by all reference-unique k31 start positions in the corresponding sequence. All 44 chromosomes are shown; Contig6 is shown separately, and the remaining 48 unplaced contigs are pooled by summing counts and opportunities. The y axis uses a symmetric-log transformation with a linear region below 0.03 to retain zeros. Shading highlights Chr19, Chr34 and Contig6 for reference, without defining new significance tests.
(B) The same unique-match density in fixed nonoverlapping 10-kb windows on Chr19, shown from 2.50 to 3.10 Mb. Shading denotes the previously specified interval Chr19:2,710,000-2,890,000 (1-based inclusive). The interval contains 32,753 female-associated and 4,507 male-associated unique k-mer start positions, among 86,367 reference-unique k31 opportunities. Counts are assigned by k-mer start coordinate. Hollow points indicate bins containing fewer than 1,000 unique starts (a display flag following the 10% opportunity criterion; no new statistical cutoff).
(C,D) Per-individual detection of the fixed female-associated (C) and male-associated (D) panels mapping uniquely within the Chr19 interval, using retained read count >=2. Each dot represents one heldout individual; vertical bars span the interquartile range and horizontal bars mark the median. Numbers above groups are medians. Panels were selected after the heldout association tests and unique mapping, so these plots summarize the selected candidates rather than providing a further independent validation.
Discovery used 8 females and 8 males; candidate k-mers were frozen before testing in 92 female and 92 male heldout individuals from the same full-sib family. All 7,178,221 fixed candidates were tested by two-sided Fisher tests and Bonferroni correction at family-wise alpha=0.05; 2,566,199 direction-consistent significant candidates were mapped by exact k31 matching on both reference strands, of which 180,082 female-associated and 184,371 male-associated candidates matched uniquely. Spatial densities are descriptive, and overlapping k-mers are not independent loci. Sex-related local depth differences and within-sex depth-dependent detection remain; sequencing batches were unrecorded. These results do not establish generalizable sex specificity, W linkage or causality. Full regional profiles and count-threshold sensitivity are appropriate supplementary figures. TE overlap is pending and is not presented here.

中文使用建议：本图可作为正文主图的新图或多面板的一部分，主旨为“重测序读段来源的k-mer关联信号与Chr19候选区重合”。建议将全基因组空间分布(A)、Chr19局部放大(B)、同家系留出个体的面板检出比例(C,D)结合展示；阈值2/3/5/10敏感性、Chr34和Contig6精细窗口图放补充材料。不能将此图描述为新的外部验证或普适性别特异标记证明。当前仅生成独立图稿，未修改现有论文或主图文件。
'''
    (OUT / 'figure_caption_and_use_notes.txt').write_text(caption, encoding='utf-8')
    inputs = [chrom_path, mapping_path, individual_path,
              BASE / 'regional_sensitivity/regional_threshold_summary.json',
              BASE / 'reference_mappability/server_download/reference_sequence_manifest.json']
    provenance = {'created_at': datetime.now().astimezone().isoformat(),
                  'source_data_status': 'Completed analyses reused; no counting, association tests or mapping rerun',
                  'analysis_scope': 'Display-only descriptive summary; original scientific parameters unchanged',
                  'mapping_rows_read': scanned, 'genome_unique_candidates_F_M': [180082, 184371],
                  'target_unique_candidates_F_M': [32753, 4507], 'target_unique_opportunities': target_opp,
                  'heldout_individuals_per_sex': 92, 'candidate_panel_selected_after_heldout_tests': True,
                  'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
                  'figure_and_source_data_sha256': {p.name: sha(p) for p in OUT.iterdir() if p.is_file() and
                                                   (p.name.startswith(('kmer_main_figure.', 'panel_A_', 'panel_B_', 'panels_CD_')) or p.name == 'figure_caption_and_use_notes.txt')},
                  'limitations': ['Single full-sib family', 'Local depth effects remain', 'Unrecorded batch',
                                  'Correlated overlapping k-mers', 'Post-selection panel summary', 'No TE results included'],
                  'versions': {'numpy': np.__version__, 'matplotlib': matplotlib.__version__}}
    (OUT / 'figure_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'output': str(OUT), 'rows_checked': scanned, 'target_F_M': dict(target)}))


if __name__ == '__main__':
    main()
