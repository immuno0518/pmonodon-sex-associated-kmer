from release_config import *
"""Display-only revision: whole-sequence profiles plus the Chr19 interval zoom."""
from create_kmer_manuscript_figure import *

OLD = BASE / 'manuscript_main_figure'
OUT = BASE / 'manuscript_main_figure_v2'


def main():
    OUT.mkdir(exist_ok=True)
    provenance = json.loads((OLD / 'figure_provenance.json').read_text())
    for name in ('panel_A_chromosome_density.tsv', 'panel_B_Chr19_10kb_density.tsv',
                 'panels_CD_heldout_detection.tsv'):
        assert sha(OLD / name) == provenance['figure_and_source_data_sha256'][name]
    arows = read_tsv(OLD / 'panel_A_chromosome_density.tsv')
    zoom = read_tsv(OLD / 'panel_B_Chr19_10kb_density.tsv')
    individuals = read_tsv(OLD / 'panels_CD_heldout_detection.tsv')
    chrom_path = BASE / 'reference_mappability/unique_opportunity_normalized_chromosomes.tsv'
    window_path = BASE / 'reference_mappability/unique_opportunity_normalized_100kb_windows.tsv'
    chroms = {r['chromosome']: r for r in read_tsv(chrom_path)}
    profiles = []
    for r in read_tsv(window_path):
        if r['chromosome'] in ('Chr19', 'Chr34'):
            profiles.append({k: r[k] for k in ('chromosome', 'start_1based', 'end_1based',
                                              'length_bp', 'unique_k31_starts', 'F_unique', 'M_unique')})
    names = identity_names()
    manifest_path = BASE / 'reference_mappability/server_download/reference_sequence_manifest.json'
    offset = 0
    contig_offset = None
    for r in json.loads(manifest_path.read_text())['sequences']:
        local = names[r['sequence']]
        assert r['length_bp'] == int(local['length_bp']) and r['sequence_sha256'] == local['sequence_sha256']
        if local['original_id'] == 'Contig6':
            contig_offset = offset
        offset += r['length_bp']
    freq = np.memmap(str(BASE / 'reference_mappability/reference_k31_e0.freq8'), mode='r', dtype=np.uint8)
    assert len(freq) == offset and contig_offset is not None
    contig_len = int(chroms['Contig6']['length_bp'])
    contig = []
    for start in range(1, contig_len + 1, 10000):
        end = min(start + 9999, contig_len)
        contig.append({'chromosome': 'Contig6', 'start_1based': start, 'end_1based': end,
                       'length_bp': end - start + 1,
                       'unique_k31_starts': int(np.count_nonzero(freq[contig_offset + start - 1:contig_offset + end] == 1)),
                       'F_unique': 0, 'M_unique': 0})
    mapping_path = BASE / 'heldout_mapping/all_replicated_candidates_reference_mapping.tsv.gz'
    scanned = 0
    with gzip.open(mapping_path, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            scanned += 1
            if r['exact_reference_hits'] == '1' and r['unique_chr'] == 'Contig6':
                contig[(int(r['unique_start_1based']) - 1) // 10000][r['candidate_group'] + '_unique'] += 1
    assert scanned == 2566199
    profiles += contig
    checks = {}
    for c in ('Chr19', 'Chr34', 'Contig6'):
        rs = [r for r in profiles if r['chromosome'] == c]
        checks[c] = {}
        for field in ('unique_k31_starts', 'F_unique', 'M_unique', 'length_bp'):
            total = sum(int(r[field]) for r in rs)
            assert total == int(chroms[c][field]), (c, field, total)
            checks[c][field] = total
        assert int(rs[0]['start_1based']) == 1
        assert all(int(x['end_1based']) + 1 == int(y['start_1based']) for x, y in zip(rs, rs[1:]))
        for r in rs:
            size = 10000 if c == 'Contig6' else 100000
            r['nominal_bin_bp'] = size
            r['low_opportunity_or_partial_bin'] = int(int(r['unique_k31_starts']) < size // 10 or int(r['length_bp']) != size)
            for g in 'FM':
                r[g + '_per_1000_unique_starts'] = int(r[g + '_unique']) / int(r['unique_k31_starts']) * 1000 if int(r['unique_k31_starts']) else ''
    write_tsv(OUT / 'panel_A_chromosome_density.tsv', arows)
    write_tsv(OUT / 'panel_B_whole_sequence_profiles.tsv', profiles)
    write_tsv(OUT / 'panel_C_Chr19_10kb_density.tsv', zoom)
    write_tsv(OUT / 'panels_DE_heldout_detection.tsv', individuals)

    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.labelsize': 8,
                         'axes.titlesize': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .7, 'xtick.major.width': .6, 'ytick.major.width': .6,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(7.48, 9.25))
    gs = fig.add_gridspec(4, 1, left=.115, right=.975, bottom=.075, top=.955,
                          height_ratios=[1, 1.04, 1.05, 1.05], hspace=.85)
    ax_a = fig.add_subplot(gs[0])
    bg = gs[1].subgridspec(1, 3, wspace=.23)
    bax = [fig.add_subplot(bg[0, i]) for i in range(3)]
    ax_c = fig.add_subplot(gs[2])
    dg = gs[3].subgridspec(1, 2, wspace=.48)
    dax = [fig.add_subplot(dg[0, i]) for i in range(2)]
    xpos = np.array(list(range(44)) + [45, 48])
    for g, dx, marker in [('F', -.14, 'o'), ('M', .14, '^')]:
        ax_a.scatter(xpos + dx, [float(r[g + '_per_1000_unique_starts']) for r in arows],
                     color=COLORS[g], s=14, marker=marker, linewidths=0,
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
    ax_a.set_title('Genome-wide chromosome-average density', loc='left', pad=8)
    ax_a.legend(frameon=False, loc='upper left', fontsize=7, ncol=2, handletextpad=.3, columnspacing=1)

    # Shared linear y-axis for the three whole-sequence profiles; every bin retained.
    maximum = max(float(r[g + '_per_1000_unique_starts']) for r in profiles for g in 'FM'
                  if r[g + '_per_1000_unique_starts'] != '')
    ytop = np.ceil(maximum / 100) * 100
    for i, (ax, c) in enumerate(zip(bax, ('Chr19', 'Chr34', 'Contig6'))):
        rs = [r for r in profiles if r['chromosome'] == c]
        xx = [(int(r['start_1based']) + int(r['end_1based'])) / 2e6 for r in rs]
        low = np.array([bool(r['low_opportunity_or_partial_bin']) for r in rs])
        for g in 'FM':
            yy = np.array([float(r[g + '_per_1000_unique_starts']) if r[g + '_per_1000_unique_starts'] != '' else np.nan for r in rs])
            ax.plot(xx, yy, color=COLORS[g], lw=.6, alpha=.22)
            ax.plot(xx, np.where(low, np.nan, yy), color=COLORS[g], lw=.8)
            ax.scatter(np.array(xx)[low], yy[low], facecolors='white', edgecolors=COLORS[g],
                       s=6, linewidths=.35, alpha=.45, zorder=3)
        ax.set_xlim(0, int(chroms[c]['length_bp']) / 1e6)
        ax.set_ylim(0, ytop * 1.06)
        ax.set_yticks(np.arange(0, ytop + 1, 100))
        if i:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel('Associated starts /\n1,000 unique starts')
        ax.set_title(c + (' | 10-kb bins' if c == 'Contig6' else ' | 100-kb bins'), loc='left', fontsize=8, pad=8)
        ax.set_xlabel('Position (Mb)')
        ax.tick_params(labelsize=7)
        if c == 'Chr19':
            ax.set_xticks([0, 10, 20, 30, 40, 50])
            ax.axvspan(2.71, 2.89, color='#D4B454', alpha=.7, zorder=0)
            ax.annotate('2.71-2.89 Mb\nzoom in C', xy=(2.8, ytop * .98), xytext=(.37, .92),
                        textcoords='axes fraction', fontsize=6.5, va='top',
                        arrowprops={'arrowstyle': '-', 'lw': .6, 'color': '#7B6424'})
        elif c == 'Chr34':
            ax.set_xticks([0, 10, 20, 30])
        else:
            ax.set_xticks([0, .1, .2, .3])
    bax[0].text(0, 1.34, 'Whole-sequence profiles (shared y scale; faint / hollow = low opportunity or partial bin)',
                transform=bax[0].transAxes, fontsize=8)

    ax_c.axvspan(2.71, 2.89, color='#F0E2B4', alpha=.7, zorder=0)
    xx = [(int(r['start_1based']) + int(r['end_1based'])) / 2e6 for r in zoom]
    for g in 'FM':
        yy = np.array([float(r[g + '_per_1000_unique_starts']) if r[g + '_per_1000_unique_starts'] != '' else np.nan for r in zoom])
        low = np.array([int(r['unique_k31_starts']) < 1000 for r in zoom])
        ax_c.plot(xx, yy, color=COLORS[g], lw=1.1, marker='o', markersize=2)
        ax_c.scatter(np.array(xx)[low], yy[low], facecolors='white', edgecolors=COLORS[g], s=17, linewidths=.6, zorder=5)
    ax_c.set_xlim(2.5, 3.1)
    ax_c.set_ylim(bottom=0)
    ax_c.set_xticks(np.arange(2.5, 3.1001, .1))
    ax_c.set_xlabel('Chr19 position (Mb; 10-kb bins)')
    ax_c.set_ylabel('Associated starts /\n1,000 unique starts')
    ax_c.set_title('Chr19 candidate interval: local detail', loc='left', pad=8)
    ax_c.text(.98, .92, 'Chr19: 2.71-2.89 Mb\n32,753 F-associated; 4,507 M-associated',
              transform=ax_c.transAxes, va='top', ha='right', fontsize=7,
              bbox={'facecolor': 'white', 'alpha': .85, 'edgecolor': 'none', 'pad': 2})
    rng = np.random.RandomState(20261021)
    for ax, g in zip(dax, 'FM'):
        for i, sex in enumerate('FM'):
            vals = np.array([float(r['detected_percent']) for r in individuals if r['candidate_group'] == g and r['sex'] == sex])
            assert len(vals) == 92
            ax.scatter(i + rng.uniform(-.16, .16, len(vals)), vals, s=10, alpha=.5,
                       color=COLORS[sex], edgecolors='none')
            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            ax.plot([i, i], [q1, q3], color='#20282E', lw=2.4, zorder=5)
            ax.plot([i - .13, i + .13], [med, med], color='#20282E', lw=1.3, zorder=6)
            ax.text(i, 108, '{:.2f}%'.format(med), ha='center', va='center', fontsize=7, fontweight='bold')
        ax.set_title(('Female' if g == 'F' else 'Male') + '-associated panel\n' + ('32,753' if g == 'F' else '4,507') + ' Chr19 interval k-mers', loc='left', pad=8)
        ax.set_xlim(-.48, 1.48)
        ax.set_ylim(-4, 116)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Female\n(n = 92)', 'Male\n(n = 92)'])
        ax.set_ylabel('Panel detected per individual (%)')
    for ax in [ax_a] + bax + [ax_c] + dax:
        ax.grid(axis='y', linewidth=.35, alpha=.25)
    for letter, ax in zip('ABCDE', [ax_a, bax[0], ax_c] + dax):
        pos = ax.get_position()
        fig.text(.025 if letter != 'E' else pos.x0 - .075,
                 pos.y1 + (.045 if letter == 'B' else .015), letter, fontsize=12, fontweight='bold')
    fig.text(.115, .018, 'One full-sib family: 16 discovery + 184 heldout individuals. Densities and selected-panel detection are descriptive.',
             fontsize=6.5, color='#4B555C')
    for fmt in ('pdf', 'svg', 'png'):
        fig.savefig(OUT / ('kmer_main_figure_v2.' + fmt), dpi=600 if fmt == 'png' else 300, facecolor='white')
    plt.close(fig)

    caption = (OLD / 'figure_caption_and_use_notes.txt').read_text(encoding='utf-8')
    start = caption.index('(B)')
    end = caption.index('Discovery used')
    caption = caption[:start] + '''(B) Whole-sequence profiles on Chr19, Chr34 and Contig6, using the same density denominator as A and a shared linear y scale. Fixed nonoverlapping 100-kb bins are used for Chr19 and Chr34; 10-kb bins are used for the shorter Contig6 (375,397 bp). Bins are assigned by 1-based k-mer start coordinate. Faint lines and hollow points flag bins with fewer than 10,000 (100-kb bins) or 1,000 (10-kb bins) reference-unique starts, or a terminal partial bin. All bins are retained; these display flags are not new statistical tests. Bin widths and x-axis ranges differ, so peak widths and heights should be interpreted at the labeled resolution. Gold shading on Chr19 marks 2.71-2.89 Mb, expanded in C.
(C) Chr19 detail from 2.50 to 3.10 Mb in fixed 10-kb bins. Gold shading denotes Chr19:2,710,000-2,890,000 (1-based inclusive). This interval contains 32,753 female-associated and 4,507 male-associated uniquely mapped candidate starts among 86,367 reference-unique k31 opportunities. Hollow points flag bins with fewer than 1,000 unique starts. C uses its own linear y scale to resolve the local peaks.
(D,E) Per-individual detection of the fixed female-associated (D) and male-associated (E) panels mapping uniquely within the Chr19 interval, using retained read count >=2. Each dot represents one heldout individual; vertical bars span the interquartile range and horizontal bars mark the median. Numbers above groups are medians. Panels were selected after the heldout association tests and unique mapping; these are post-selection summaries, not a further independent validation.
''' + caption[end:]
    caption = caption[:caption.index('\n\n中文使用建议')]
    caption = caption.replace('Full regional profiles and count-threshold sensitivity are appropriate supplementary figures.',
                              'Count-threshold sensitivity and detailed low-mappability audits are appropriate supplementary figures.')
    caption += '\n\n中文使用建议：A保留全基因组序列平均密度；B并列展示Chr19、Chr34、Contig6全长分布，使用同一纵轴，Chr19目标区金色突出；C保留目标区10 kb放大；D/E为同家系留出个体的后选择面板检出概览。B中两条染色体采用100 kb窗口，短Contig6采用10 kb窗口，不能将不同窗口下的峰宽直接比较。低可比对性窗口保留并淡化。原始分析参数和有效结果未改，v1图稿保留。\n'
    (OUT / 'figure_caption_and_use_notes.txt').write_text(caption, encoding='utf-8')
    sources = [OLD / n for n in ('panel_A_chromosome_density.tsv', 'panel_B_Chr19_10kb_density.tsv', 'panels_CD_heldout_detection.tsv')]
    sources += [chrom_path, window_path, mapping_path, manifest_path]
    result = {'created_at': datetime.now().astimezone().isoformat(), 'revision': 'v2: three whole-sequence profiles plus Chr19 zoom',
              'source_data_status': 'Completed results reused; no counting, association tests or reference mapping rerun',
              'mapping_rows_read': scanned, 'whole_profile_sum_checks': checks,
              'display_bins_bp': {'Chr19': 100000, 'Chr34': 100000, 'Contig6': 10000, 'Chr19_zoom': 10000},
              'B_shared_linear_y_limit': float(ytop * 1.06),
              'limitations': provenance['limitations'],
              'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sources},
              'figure_and_source_data_sha256': {p.name: sha(p) for p in OUT.iterdir() if p.suffix in ('.pdf', '.svg', '.png', '.tsv', '.txt') and not p.name.startswith('pdf_render')}}
    (OUT / 'figure_provenance.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'output': str(OUT), 'checks': checks, 'B_max_rate': maximum}))


if __name__ == '__main__':
    main()
