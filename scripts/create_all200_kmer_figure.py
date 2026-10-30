from release_config import *
"""Display-only revision: whole-sequence profiles plus the Chr19 interval zoom."""
from create_kmer_manuscript_figure import *

OLD = BASE / 'manuscript_main_figure_v2'
OUT = BASE / 'manuscript_main_figure_v3'


def main():
    OUT.mkdir(exist_ok=True)
    provenance = json.loads((OLD / 'figure_provenance.json').read_text())
    for name, digest in provenance['figure_and_source_data_sha256'].items():
        assert sha(OLD / name) == digest
    validation = json.loads((OUT / 'discovery_panel_validation.json').read_text())
    assert validation['status'] == 'passed'
    assert sha(OUT / 'panels_DE_all200_detection.tsv') == validation['all200_table_sha256']
    arows = read_tsv(OLD / 'panel_A_chromosome_density.tsv')
    zoom = read_tsv(OLD / 'panel_C_Chr19_10kb_density.tsv')
    profiles = read_tsv(OLD / 'panel_B_whole_sequence_profiles.tsv')
    for r in profiles:
        r['low_opportunity_or_partial_bin'] = int(r['low_opportunity_or_partial_bin'])
    individuals = read_tsv(OUT / 'panels_DE_all200_detection.tsv')
    chroms = {r['chromosome']: r for r in read_tsv(BASE / 'reference_mappability/unique_opportunity_normalized_chromosomes.tsv')}
    for name in ('panel_A_chromosome_density.tsv', 'panel_B_whole_sequence_profiles.tsv', 'panel_C_Chr19_10kb_density.tsv'):
        (OUT / name).write_bytes((OLD / name).read_bytes())

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
            assert len(vals) == 100
            selected = [r for r in individuals if r['candidate_group'] == g and r['sex'] == sex]
            assert Counter(r['cohort'] for r in selected) == {'Discovery': 8, 'Heldout': 92}
            for cohort, marker, pointsize, alpha in [('Heldout', 'o', 10, .5), ('Discovery', 'D', 22, 1)]:
                vv = [float(r['detected_percent']) for r in selected if r['cohort'] == cohort]
                ax.scatter(i + rng.uniform(-.16, .16, len(vv)), vv, s=pointsize, alpha=alpha,
                           marker=marker, facecolors=COLORS[sex] if cohort == 'Heldout' else 'none',
                           edgecolors='none' if cohort == 'Heldout' else '#20282E', linewidths=.7,
                           zorder=3 if cohort == 'Heldout' else 4)
            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            ax.plot([i, i], [q1, q3], color='#20282E', lw=2.4, zorder=5)
            ax.plot([i - .13, i + .13], [med, med], color='#20282E', lw=1.3, zorder=6)
            ax.text(i, 108, '{:.2f}%'.format(med), ha='center', va='center', fontsize=7, fontweight='bold')
        ax.set_title(('Female' if g == 'F' else 'Male') + '-associated | ' + ('32,753' if g == 'F' else '4,507') + ' k-mers\nChr19:2.71-2.89 Mb', loc='left', pad=8)
        ax.set_xlim(-.48, 1.48)
        ax.set_ylim(-4, 116)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Female\n(n = 100)', 'Male\n(n = 100)'])
        ax.set_ylabel('Panel detected per individual (%)')
    for ax in [ax_a] + bax + [ax_c] + dax:
        ax.grid(axis='y', linewidth=.35, alpha=.25)
    for letter, ax in zip('ABCDE', [ax_a, bax[0], ax_c] + dax):
        pos = ax.get_position()
        fig.text(.025 if letter != 'E' else pos.x0 - .075,
                 pos.y1 + (.045 if letter == 'B' else .015), letter, fontsize=12, fontweight='bold')
    fig.legend(handles=[Line2D([], [], marker='o', color='none', markerfacecolor='#808890', markeredgecolor='none', markersize=4, label='Heldout (92 / sex)'),
                        Line2D([], [], marker='D', color='none', markerfacecolor='none', markeredgecolor='#20282E', markersize=4, label='Discovery (8 / sex)')],
               loc='lower center', bbox_to_anchor=(.55, .024), ncol=2, frameon=False, fontsize=7)
    fig.text(.115, .010, 'One full-sib family; all 200 individuals shown in D/E. Selected-panel detection is descriptive, not an additional validation.',
             fontsize=6.3, color='#4B555C')
    for fmt in ('pdf', 'svg', 'png'):
        fig.savefig(OUT / ('kmer_main_figure_v3.' + fmt), dpi=600 if fmt == 'png' else 300, facecolor='white')
    plt.close(fig)

    caption = (OLD / 'figure_caption_and_use_notes.txt').read_text(encoding='utf-8')
    a = caption.index('(D,E)'); b = caption.index('Discovery used')
    caption = caption[:a] + """(D,E) Descriptive per-individual detection of the fixed female-associated (D; 32,753 k-mers) and male-associated (E; 4,507 k-mers) panels mapping uniquely within Chr19:2,710,000-2,890,000 (1-based inclusive), at retained read count >=2. All 200 individuals are displayed: 100 females and 100 males, each comprising 8 discovery and 92 heldout individuals. Filled circles identify heldout individuals and open diamonds identify discovery individuals. Vertical bars span the interquartile range across all 100 individuals of the corresponding sex; horizontal bars and numbers indicate the pooled median. Discovery observations were extracted from the existing per-individual k31 databases, and their focal-sex per-k-mer prevalences were independently matched to the original discovery screen. No reads were recounted and no association test was rerun. The panels were selected using discovery screening followed by heldout association and unique mapping; the all-200 display is a post-selection descriptive summary, not a 200-individual validation set or a further independent validation.
""" + caption[b:]
    caption = caption[:caption.index('\n\n')] + '\n'
    (OUT / 'figure_caption_and_use_notes.txt').write_text(caption, encoding='utf-8')
    sources = [OLD / n for n in ('panel_A_chromosome_density.tsv', 'panel_B_whole_sequence_profiles.tsv', 'panel_C_Chr19_10kb_density.tsv')]
    sources += [OUT / 'panels_DE_all200_detection.tsv', OUT / 'discovery_panel_validation.json', OUT / 'all200_display_plan.json']
    result = {'created_at': datetime.now().astimezone().isoformat(), 'revision': 'v3: D/E include all 100 individuals per sex, discovery and heldout distinguished',
              'source_data_status': 'Existing discovery DBs queried; no counting, association tests or reference mapping rerun',
              'whole_profile_sum_checks': provenance['whole_profile_sum_checks'],
              'cohorts_per_sex': {'Discovery': 8, 'Heldout': 92}, 'n_per_sex': 100,
              'target_interval': 'Chr19:2710000-2890000 (1-based inclusive)',
              'fixed_panel_sizes': {'F': 32753, 'M': 4507},
              'limitations': provenance['limitations'] + ['Pooled discovery and heldout samples are not a 200-individual validation cohort'],
              'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sources},
              'figure_and_source_data_sha256': {p.name: sha(p) for p in OUT.iterdir() if p.suffix in ('.pdf', '.svg', '.png', '.tsv', '.txt') and not p.name.startswith('pdf_render')}}
    (OUT / 'figure_provenance.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'output': str(OUT), 'points_per_panel': 200}))


if __name__ == '__main__':
    main()
