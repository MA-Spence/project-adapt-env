"""Static graph descriptors for ANA-072; no learning, fitting or search task.

Neighbor moments are recomputed from fixed maps, in O(L*q**L) using the
complete-clique structure of each single-site fiber. Published RUN-186 graph
summaries supply reachability and neutrality; no graph policy is rerun.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from average_tools import MODES, MODE_LABELS, iter_maps, save_figure

LAYER_COLORS = {'intrinsic': '#247aaf', 'observed': '#cf7026'}
GROUPS = [('none', 'intrinsic'), ('independent_difference', 'intrinsic'),
          ('independent_difference', 'observed'), ('shared_stop_sensitivity', 'intrinsic'),
          ('shared_stop_sensitivity', 'observed')]
GROUP_LABELS = ['No error, r = 4\nIntrinsic', 'Independent error, r = 2\nIntrinsic',
                'Independent error, r = 2\nObserved', 'Shared-stop, r = 0.5\nIntrinsic',
                'Shared-stop, r = 0.5\nObserved']


def edge_moments(values, valid, q=20, dim=4):
    """Moments on the undirected, measured-mask Hamming graph.

    Each single-site fiber is a clique. With valid count n, centered sum s and
    square-sum t, the sum of squared edge differences is n*t-s*s. Each node
    appears n-1 times among the endpoints in that fiber. We symmetrize edge
    orientation, so the Pearson correlation is 1 - mean_delta2/(2*var_endpoint).
    """
    values, valid = np.asarray(values, float), np.asarray(valid, bool)
    assert values.size == valid.size == q ** dim
    assert np.isfinite(values[valid]).all()
    center = float(values[valid].mean())
    v = np.where(valid, values - center, 0.).reshape((q,) * dim)
    mask = valid.reshape((q,) * dim)
    edges = 0
    delta2 = endpoint_sum = endpoint_square_sum = 0.
    for axis in range(dim):
        n = mask.sum(axis=axis)
        s = v.sum(axis=axis)
        t = np.square(v).sum(axis=axis)
        edges += int(np.sum(n * (n - 1) // 2))
        delta2 += float(np.sum(n * t - s * s))
        endpoint_sum += float(np.sum((n - 1) * s))
        endpoint_square_sum += float(np.sum((n - 1) * t))
    assert edges > 0
    edge_delta2 = max(0., delta2 / edges)
    endpoint_variance = endpoint_square_sum / (2 * edges) - (endpoint_sum / (2 * edges)) ** 2
    node_variance = float(values[valid].var())
    return dict(nodes=int(valid.sum()), edges=edges,
                neighbor_correlation=float(1 - edge_delta2 / (2 * endpoint_variance))
                    if endpoint_variance > 0 else float('nan'),
                edge_rms_native=float(np.sqrt(edge_delta2)),
                normalized_edge_rms=float(np.sqrt(edge_delta2 / node_variance))
                    if node_variance > 0 else float('nan'),
                endpoint_variance=endpoint_variance, node_variance=node_variance)


def validate_edge_moments():
    """Compare against explicit edge enumeration on small independent fixtures."""
    checks = []
    for q, dim in [(2, 3), (3, 2), (4, 2)]:
        codes = np.array(np.unravel_index(np.arange(q ** dim), (q,) * dim)).T
        for missing in [False, True]:
            valid = np.ones(q ** dim, dtype=bool)
            if missing:
                valid[[1, -2]] = False
            f = 1000. + np.sin(np.arange(q ** dim) * .7) + np.arange(q ** dim) ** 2 / 7
            edges = [(i, j) for i in range(len(f)) for j in range(i + 1, len(f))
                     if valid[i] and valid[j] and np.sum(codes[i] != codes[j]) == 1]
            i, j = np.array(edges).T
            a, b = np.r_[f[i], f[j]], np.r_[f[j], f[i]]
            result = edge_moments(f, valid, q, dim)
            assert result['edges'] == len(edges)
            np.testing.assert_allclose(result['neighbor_correlation'], np.corrcoef(a, b)[0, 1], atol=1e-12)
            np.testing.assert_allclose(result['edge_rms_native'], np.sqrt(np.mean((f[i] - f[j]) ** 2)), atol=1e-12)
            np.testing.assert_allclose(result['normalized_edge_rms'],
                                       np.sqrt(np.mean((f[i] - f[j]) ** 2)) / np.std(f[valid]), atol=1e-12)
            checks.append(dict(q=q, dim=dim, missing=missing, edges=len(edges)))
    return checks


def graph_table(ctx):
    table = pd.read_csv(ctx.output / 'graph_summary.csv')
    manifest = json.loads((ctx.output / 'manifest.json').read_text())
    assert hashlib.sha256((ctx.output / 'graph_summary.csv').read_bytes()).hexdigest() == manifest['files']['graph_summary.csv']
    assert len(table) == 363
    assert set(table.relative_resolution) == {0., .001, .01}
    assert set(table.nodes) == {159129} and set(table.edges) == {6016433}
    table['inaccessible_top_fraction'] = 1 - table.accessible_top
    # Components restricted to above-floor nodes use an above-floor denominator.
    # Also supply a common all-node denominator for unambiguous table reuse.
    table['above_floor_component_fraction_all_nodes'] = (table.above_floor_largest_fraction
                                                        * table.above_floor_fraction)
    return table


def summarize(frame, metrics, keys):
    records = []
    for values, g in frame.groupby(keys, sort=False):
        if not isinstance(values, tuple):
            values = (values,)
        for metric in metrics:
            x = g[metric].dropna().to_numpy()
            if not len(x):
                continue
            records.append(dict(zip(keys, values)) | dict(metric=metric, n=len(x),
                minimum=float(x.min()), median=float(np.median(x)), maximum=float(x.max()),
                q10=float(np.quantile(x, .1)), q90=float(np.quantile(x, .9))))
    return pd.DataFrame(records)


def all_edge_moments(ctx):
    expected = np.array(np.unravel_index(np.arange(20 ** 4), (20,) * 4)).T
    np.testing.assert_array_equal(ctx.genotypes, expected)
    rows = [dict(seed=-1, cohort='empirical', observation='measured', ratio=-1.,
                 layer='empirical', **edge_moments(ctx.empirical, ctx.valid))]
    for meta, values in iter_maps(ctx):
        rows.append(dict(**meta, **edge_moments(values, ctx.valid)))
    result = pd.DataFrame(rows)
    assert len(result) == 121
    assert set(result.nodes) == {159129} and set(result.edges) == {6016433}
    return result


def zero_threshold_invariance(ctx, graph):
    """Check archived monotone-readout rank descriptors, not a new walk."""
    latent = pd.read_csv(ctx.output / 'latent_graph_summary.csv')
    rows = []
    for g in graph[(graph.layer == 'intrinsic') & (graph.relative_resolution == 0)].itertuples():
        a = latent[(latent.seed == g.seed) & (latent.ratio == g.ratio)].iloc[0]
        for metric in ['no_beneficial_fraction', 'accessible_top', 'accessible_best']:
            difference = abs(float(getattr(g, metric)) - float(a[metric]))
            np.testing.assert_allclose(getattr(g, metric), a[metric], rtol=0., atol=1e-14)
            rows.append(dict(seed=g.seed, observation=g.observation, ratio=g.ratio,
                             metric=metric, absolute_difference=difference))
    return rows


def figure_static(ctx, graph, edges):
    """All fresh seeds are shown; reused points remain separate and subdued."""
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 6.0), sharey=True)
    config = [('neighbor_correlation', 'Neighbor fitness correlation\n(edge-endpoint Pearson r)', edges, 1.),
              ('no_beneficial_fraction', 'No beneficial neighbor (%)\nNonlinear scale',
               graph[graph.relative_resolution == .001], 100.),
              ('inaccessible_top_fraction', 'No uphill path to top 1% (%)\nNonlinear scale',
               graph[graph.relative_resolution == .001], 100.)]
    for ax, (metric, title, data, factor) in zip(axes, config):
        empirical = float(data.loc[data.cohort == 'empirical', metric].iloc[0]) * factor
        ax.axvline(empirical, color='#222222', linestyle='--', lw=1.5, zorder=0)
        ax.text(empirical, 1.025, f'TrpB {empirical:.3g}', transform=ax.get_xaxis_transform(),
                ha='center', va='bottom', fontsize=9)
        for y, (mode, layer) in enumerate(GROUPS):
            select = data[(data.observation == mode) & (data.layer == layer)]
            fresh = select.loc[select.cohort == 'fresh', metric].to_numpy() * factor
            reused = select.loc[select.cohort == 'reused', metric].to_numpy() * factor
            assert len(fresh) == 8 and len(reused) == 16
            ax.scatter(reused, y + np.linspace(-.19, -.07, len(reused)), s=14,
                       facecolors='none', edgecolors='#aaaaaa', linewidths=.65, zorder=1)
            ax.plot([fresh.min(), fresh.max()], [y, y], color=LAYER_COLORS[layer], lw=2, alpha=.5)
            ax.scatter(fresh, y + np.linspace(.02, .16, len(fresh)), s=23,
                       color=LAYER_COLORS[layer], alpha=.8, zorder=3)
            ax.scatter([np.median(fresh)], [y], s=65, color=LAYER_COLORS[layer], marker='|',
                       linewidths=2.4, zorder=4)
        ax.set_xlabel(title, fontsize=10)
        ax.grid(axis='x', alpha=.2)
        ax.set_axisbelow(True)
        ax.set_ylim(4.5, -.55)
    axes[0].set_yticks(range(len(GROUPS)), GROUP_LABELS, fontsize=10)
    axes[0].set_xlim(-.05, 1.0)
    axes[1].set_xscale('symlog', linthresh=.01)
    axes[1].set_xlim(left=0)
    axes[1].set_xticks([0, .01, .1, 1, 10, 40], ['0', '0.01', '0.1', '1', '10', '40'])
    axes[2].set_xscale('symlog', linthresh=.001)
    axes[2].set_xlim(left=0)
    axes[2].set_xticks([0, .001, .01, .1, 1, 10, 40], ['0', '0.001', '0.01', '0.1', '1', '10', '40'])
    axes[1].tick_params(axis='x', labelsize=9)
    axes[2].tick_params(axis='x', labelsize=9)
    fig.suptitle('Static structure: resemblance depends on the observation assumption', fontsize=15, y=.99)
    handles = [Line2D([], [], color=LAYER_COLORS['intrinsic'], marker='o', ls='', label='Intrinsic'),
               Line2D([], [], color=LAYER_COLORS['observed'], marker='o', ls='', label='Observed'),
               Line2D([], [], color='#aaaaaa', marker='o', markerfacecolor='none', ls='', label='16 reused seeds'),
               Line2D([], [], color='#222222', ls='--', label='Empirical TrpB')]
    fig.legend(handles=handles, loc='lower center', ncol=4, bbox_to_anchor=(.56, .055), frameon=False)
    fig.text(.5, .012, 'Filled points: 8 fresh seeds; bars: fresh-seed min–max and median. '
             'Navigation panels use Δfitness > 0.001 × calibration range.\n'
             'Navigation axes use a logarithmic scale above their linear near-zero interval. Same measured graph for every map; these are not task-performance scores.',
             ha='center', va='bottom', fontsize=9)
    fig.subplots_adjust(left=.18, right=.98, top=.86, bottom=.23, wspace=.22)
    save_figure(fig, ctx.analysis / 'figures', '05_static_learning_and_navigation_proxies')
    return fig


def figure_resolution(ctx, graph):
    metrics = [('neutral_edge_fraction', 'Neutral edges (%)'),
               ('all_largest_fraction', 'Largest neutral component (%)\nAll measured nodes'),
               ('above_floor_largest_fraction', 'Largest neutral component (%)\nAbove-floor nodes only')]
    fig, axes = plt.subplots(3, 3, figsize=(12.6, 9.0), sharex=True, sharey='row')
    xs = np.arange(3)
    for col, mode in enumerate(MODES):
        for row, (metric, ylabel) in enumerate(metrics):
            ax = axes[row, col]
            empirical = graph[graph.cohort == 'empirical'].sort_values('relative_resolution')[metric].to_numpy() * 100
            ax.plot(xs, empirical, color='#222222', lw=1.6, linestyle='--', marker='s', markersize=3, zorder=5)
            for layer in ['intrinsic', 'observed']:
                selection = graph[(graph.observation == mode) & (graph.layer == layer)]
                if selection.empty:
                    continue
                fresh = selection[selection.cohort == 'fresh'].groupby('relative_resolution')[metric]
                lo = fresh.min().to_numpy() * 100
                mid = fresh.median().to_numpy() * 100
                hi = fresh.max().to_numpy() * 100
                reused = selection[selection.cohort == 'reused'].groupby('relative_resolution')[metric].median().to_numpy() * 100
                ax.fill_between(xs, lo, hi, color=LAYER_COLORS[layer], alpha=.14)
                ax.plot(xs, mid, color=LAYER_COLORS[layer], marker='o', markersize=4, lw=1.8)
                ax.plot(xs, reused, color=LAYER_COLORS[layer], lw=.8, ls=':', alpha=.75)
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=10)
            if row == 0:
                ratio = next(x['ratio'] for x in ctx.settings if x['observation'] == mode)
                ax.set_title(f'{MODE_LABELS[mode]}\nr = {ratio:g}', fontsize=11)
            ax.grid(alpha=.2)
            ax.set_ylim(-3, max(3., 105 * float(graph[metric].max())))
            if row == 2:
                ax.set_xticks(xs, ['0', '0.001', '0.01'])
                ax.set_xlabel('Resolution / calibration range', fontsize=10)
    fig.suptitle('Apparent neutrality depends strongly on readout and measurement resolution', fontsize=14, y=.995)
    handles = [Line2D([], [], color=LAYER_COLORS['intrinsic'], marker='o', label='Intrinsic: fresh median'),
               Line2D([], [], color=LAYER_COLORS['observed'], marker='o', label='Observed: fresh median'),
               Line2D([], [], color='#222222', ls='--', marker='s', label='Empirical TrpB'),
               Line2D([], [], color='#777777', ls=':', label='Reused median')]
    fig.legend(handles=handles, loc='lower center', ncol=4, bbox_to_anchor=(.51, .055), frameon=False, fontsize=9)
    fig.text(.5, .009, 'Shading: min–max across 8 fresh seeds. All three resolutions were fixed before fitting. '
             'Neutral means |Δfitness| ≤ ε.\n'
             'A connected chain of neutral edges may span > ε in total fitness; a large component is not proof of a flat landscape.\n'
             'Above-floor panels exclude nodes at or below the fixed calibration floor; their denominator differs from the middle row.',
             ha='center', va='bottom', fontsize=9)
    fig.subplots_adjust(left=.10, right=.985, top=.91, bottom=.18, hspace=.20, wspace=.10)
    save_figure(fig, ctx.analysis / 'figures', '06_resolution_and_neutral_connectivity')
    return fig
