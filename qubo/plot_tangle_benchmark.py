import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from qubo_solvers.plotting_helpers import (
    radar_factory, 
    union_axes_bbox,
    get_compare_solver_data, 
    get_compare_annotator_data, 
    get_violin_plot_data,
    relocate_radar_labels,
    place_row_legend,
)


solver_names = {
    "dwave": "D-Wave",
    "mqlib": "MQLib",
    "gurobi": "Gurobi",
    "pathfinder": "Pathfinder"
}
annotator_names = {
    "mg": "Minigraph",
    "km": "Kmer2node",
    "ga": "GraphAligner"
}
solver_colours = {
     "gurobi": '#5366E0',
    "mqlib": '#13553A',
    "pathfinder": '#911449',
    "dwave": "#FD8153"
}

dtype="cons"
data_folder = '/nfs/users/nfs_j/jc59/quantumwork/pangenome/out/pathfinder_copy_numbers_harder_cons_20_5.60.120.300_3_5'
solvers = ["gurobi", "mqlib", "pathfinder"]
# data_folder = '/nfs/users/nfs_j/jc59/quantumwork/pangenome/out/remake_dwave'
# solvers = ["dwave", "pathfinder"]
num_annotators = len(annotator_names.keys())




# ---------------------------
# Main builder: now supports multiple violin panels in bottom row
# ---------------------------
def make_combined_figure(
    plot_radar_top,         # function(axs_list) draws into list of axes
    plot_radar_mid,         # function(axs_list)
    plot_violin_fn,         # function(axs_list)   <-- now accepts list of violin axes
    ncols_top=3,
    ncols_mid=3,
    ncols_violin=3,
    width_in=6.27,
    dpi=300,
    top_height_in=2.8,
    mid_height_in=2.8,
    violin_height_in=2.4,
    caption_space_in=0.22,
    caption_fontsize=9,
    small_gap_in=0.06,
    min_side_margin_in=0.25,  # minimal left/right margin in inches
    fig_title=None,
    savepath=None
):
    """
    Layout where each panel (regardless of row) has the same physical width.
    The panel width is computed from width_in and the maximum number of columns
    among the rows, leaving min_side_margin_in margins on left and right.
    """

    # total height in inches (stacks rows + caption space)
    total_h_in = (top_height_in + caption_space_in + small_gap_in +
                  mid_height_in + caption_space_in + small_gap_in +
                  violin_height_in + caption_space_in)

    # Determine panel width in inches using the max number of columns
    max_ncols = max(ncols_top, ncols_mid, ncols_violin)

    # available width for panels after guaranteeing minimal side margins
    usable_w_in = width_in - 2 * min_side_margin_in
    if usable_w_in <= 0:
        raise ValueError("width_in too small for given min_side_margin_in")

    # panel width such that max_ncols panels fit in usable width
    panel_width_in = usable_w_in / max_ncols

    # If panel_width_in is tiny, reduce min_side_margin_in to allow more width (optional)
    if panel_width_in < 0.6:   # arbitrary lower bound for reasonable panels
        # recompute with smaller side margin
        min_side_margin_in = max(0.05, (width_in - max_ncols * 0.6) / 2.0)
        usable_w_in = width_in - 2 * min_side_margin_in
        panel_width_in = usable_w_in / max_ncols

    # Helper to compute left/right figure coords that center ncols panels with width panel_width_in
    def row_lr(ncols):
        total_panels_w_in = panel_width_in * ncols
        side_margin_in = max(min_side_margin_in, (width_in - total_panels_w_in) / 2.0)
        left_in = side_margin_in
        right_in = width_in - side_margin_in
        # convert to fraction of figure width
        left_frac = left_in / width_in
        right_frac = right_in / width_in
        return left_frac, right_frac

    # Compute top/mid/violin vertical fractions (figure coordinates)
    # We'll use top anchor at 0.95 and work downwards
    fig = plt.figure(figsize=(width_in, total_h_in), dpi=dpi)

    top_top = 1.0
    top_bottom = top_top - (top_height_in) / total_h_in

    mid_top = top_bottom # - small_gap_in / total_h_in
    mid_bottom = mid_top - (mid_height_in) / total_h_in

    violin_top = mid_bottom - 3*small_gap_in / total_h_in
    violin_bottom = violin_top - (violin_height_in) / total_h_in
    print(top_top, top_bottom, mid_top, mid_bottom, violin_top, violin_bottom)

    theta = radar_factory(7, frame='polygon')
    
    # build top row with computed left/right so each panel has same physical width
    left_top, right_top = row_lr(ncols_top)
    gs_top = GridSpec(1, ncols_top, left=left_top, right=right_top, top=top_top, bottom=top_bottom, wspace=1.25*min_side_margin_in)
    axs_top = []
    for i in range(ncols_top):
        ax = fig.add_subplot(gs_top[0, i], projection='radar', sharey=axs_top[0] if len(axs_top) else None)
        axs_top.append(ax)
    axs_top = np.array(axs_top)
    plot_radar_top(axs_top)
    relocate_radar_labels(axs_top, theta, frame='polygon', r_label=1.3, fontsize=8, wrap_len=8)
    handles, labels = axs_top[-1].get_legend_handles_labels()
    place_row_legend(fig, axs_top, handles, labels, loc='upper center', fontsize=8)
    # place_row_legend_inset(fig, axs_top, handles, labels, ncol=1, fontsize=8)

    
    # middle row
    left_mid, right_mid = row_lr(ncols_mid)
    gs_mid = GridSpec(1, ncols_mid, left=left_mid, right=right_mid, top=mid_top, bottom=mid_bottom, wspace=1.25*min_side_margin_in)
    axs_mid = []
    for i in range(ncols_mid):
        ax = fig.add_subplot(gs_mid[0, i], projection='radar', sharey=axs_mid[0] if len(axs_mid) else None)
        axs_mid.append(ax)
    axs_mid = np.array(axs_mid)
    plot_radar_mid(axs_mid)
    relocate_radar_labels(axs_mid, theta, frame='polygon', r_label=1.25, fontsize=8, wrap_len=8)
    handles, labels = axs_mid[-1].get_legend_handles_labels()
    place_row_legend(fig, axs_mid, handles, labels, loc='upper center', fontsize=8)
    
    # violin row
    left_v, right_v = row_lr(ncols_violin)
    gs_violin = GridSpec(1, ncols_violin, left=left_v, right=right_v, top=violin_top, bottom=violin_bottom, wspace=min_side_margin_in)
    axs_violin = [fig.add_subplot(gs_violin[0, 0])]
    axs_violin.extend([fig.add_subplot(gs_violin[0, i], sharey=axs_violin[0]) for i in range(1, ncols_violin)])
    plot_violin_fn(axs_violin)
    for ax in axs_violin:
        ax.tick_params(axis='x', pad=6)

    # Subcaptions: center under each row's union bbox
    bbox_top = union_axes_bbox(axs_top)
    caption_y_top = bbox_top.y0 - (caption_space_in / total_h_in) *1.05
    fig.text(bbox_top.x0 + bbox_top.width * 0.5, caption_y_top,
             "(a) " + (fig_title[0] if fig_title and isinstance(fig_title, (list,tuple)) else "Top radar set"),
             ha='center', va='top', fontsize=caption_fontsize)

    bbox_mid = union_axes_bbox(axs_mid)
    caption_y_mid = bbox_mid.y0 - (caption_space_in / total_h_in) *1.05
    fig.text(bbox_mid.x0 + bbox_mid.width * 0.5, caption_y_mid,
             "(b) " + (fig_title[1] if fig_title and isinstance(fig_title, (list,tuple)) else "Middle radar set"),
             ha='center', va='top', fontsize=caption_fontsize)

    bbox_violin = union_axes_bbox(axs_violin)
    caption_y_violin = bbox_violin.y0 - (caption_space_in / total_h_in) * 1.1
    fig.text(bbox_violin.x0 + bbox_violin.width * 0.5, caption_y_violin,
             "(c) " + (fig_title[2] if fig_title and isinstance(fig_title, (list,tuple)) else "Violin panels"),
             ha='center', va='top', fontsize=caption_fontsize)

    if fig_title and isinstance(fig_title, str):
        fig.suptitle(fig_title, y=0.99, fontsize=10)

    if savepath:
        fig.savefig(savepath, dpi=dpi)
    return fig, (axs_top, axs_mid, axs_violin)


if __name__ == "__main__":
    def plot_compare_annotators(axs):
        data, spoke_labels, labels = get_compare_annotator_data(data_folder, dtype, solvers)
        N = data[0][1].shape[1] 
        num_time_limits = int((data[0][1].shape[0] - 1)/(len(solvers) - 1))

        
        theta = radar_factory(N, frame='polygon')
        
        colors = [solver_colours[solver] for solver in solvers]
        styles = [':', '-.', '--', '-']

        colors_for_lines = list(
            np.array(
                [[colors[i]] * num_time_limits for i in range(len(solvers))]
            ).reshape((num_time_limits * len(solvers),))
        ) + [colors[-1]]
        styles_for_lines = styles[-num_time_limits:] * (len(solvers)-1) + [styles[-1]]
        
        for ax, (title, case_data) in zip(axs.flat, data):
            ax.set_rgrids([0.2, 0.4, 0.6, 0.8, 1])
            ax.set_title(annotator_names[title], weight='bold', size='medium', y=1.23,
                            horizontalalignment='center')
            for i in range(case_data.shape[0]):
                ax.plot(theta, case_data[i,:], color=colors_for_lines[i], linestyle=styles_for_lines[i], label=labels[i])  
            ax.set_varlabels(spoke_labels)
            ax.set_yticklabels([])

        
    
    def plot_compare_solvers(axs):
        data, spoke_labels, labels = get_compare_solver_data(data_folder, dtype, solvers)
        N = data[0][1].shape[1] 
        theta = radar_factory(N, frame='polygon')
        
        colors_for_lines = list(solver_colours.values())[:num_annotators]
        styles_for_lines = ['-'] * num_annotators

        for ax, (title, case_data) in zip(axs.flat, data):
            ax.set_rgrids([0.2, 0.4, 0.6, 0.8, 1])
            for i in range(case_data.shape[0]):
                ax.plot(theta, case_data[i,:], color=colors_for_lines[i], linestyle=styles_for_lines[i], label=labels[i])  
            ax.set_title(title, weight='bold', size='medium', y=1.23,
                horizontalalignment='center')
            ax.set_varlabels(spoke_labels)
            ax.set_yticklabels([])

        # axs[-1].legend(labels, loc='upper left', bbox_to_anchor=(0.85, 1.40 if len(solvers) == 2 else 1.45), # (0.85, 1.55)
        #                     borderaxespad=0.2, labelspacing=0.1, fontsize='small')
    

    def plot_violin(axs):
        violin_width = 0.6
        data, line_labels = get_violin_plot_data(data_folder, dtype, solvers)
        for ax_idx, ax in enumerate(axs):
            col = solver_colours[solvers[ax_idx]]
            slice_data = data[(ax_idx)*num_annotators:(ax_idx+1)*num_annotators]
            xticks = line_labels[(ax_idx)*num_annotators:(ax_idx+1)*num_annotators]
            title = solver_names[solvers[ax_idx]]

            parts = ax.violinplot(slice_data, showmeans=True, widths=violin_width, points=100)
            for pc in parts['bodies']:
                pc.set_facecolor(col); pc.set_edgecolor('black'); pc.set_alpha(0.85)
            for name in ('cmeans', 'cbars', 'cmins', 'cmaxes'):
                if name in parts:
                    parts[name].set_color('black')

            ax.set_xticks(np.arange(1, len(xticks)+1))
            ax.set_xticklabels(xticks, ha='center', fontsize='8')  
            for i, label in enumerate(ax.get_xticklabels()):
                label.set_y(-0.0 if i % 2 == 0 else -0.08)
            ax.set_title(title, weight='bold', size='medium', position=(0.5, 1.5),
                        horizontalalignment='center')

        axs[0].set_ylabel(r'$(\text{Covered} + \text{Used})\ / \ 2$')
        

    # Build and save combined figure
    fig, (axs_top, axs_mid, axs_violin) = make_combined_figure(
        plot_compare_annotators,
        plot_compare_solvers,
        plot_violin,
        ncols_top=num_annotators,
        ncols_mid=len(solvers),
        ncols_violin=len(solvers),
        width_in=6.27,
        dpi=300,
        top_height_in=2.6,
        mid_height_in=2.6,
        violin_height_in=1.6,
        caption_space_in=0.45,
        caption_fontsize=8,
        small_gap_in=0.1,
        min_side_margin_in=0.4,
        fig_title= [
            "Radar charts comparing the performance of each of the annotation strategies.",
            "Radar charts comparing the performance of each of the classical solvers.",
            "Violin plots comparing the performance of each combination of annotation strategy and classical solver."
        ],
        # fig_title= [
        #     "Radar charts comparing the performance of each of the annotation strategies.",
        #     "Radar charts comparing the performance of D-Wave and pathfinder",
        #     "Violin plots comparing the performance of D-Wave and pathfinder for each annotation strategy."
        # ],
        savepath=f'{data_folder}/new_combined.{".".join(solvers)}.png',
    )

    plt.show()

