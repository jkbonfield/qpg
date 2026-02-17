import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
from matplotlib.transforms import Bbox

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


def union_axes_bbox(axes_list):
    bboxes = [ax.get_position() for ax in axes_list]
    bbox_union = bboxes[0]
    for b in bboxes[1:]:
        bbox_union = Bbox.union([bbox_union, b])
    return bbox_union


def radar_factory(num_vars, frame='circle'):
    """
    Create a radar chart with `num_vars` Axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle', 'polygon'}
        Shape of frame surrounding Axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):

        def transform_path_non_affine(self, path):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):

        name = 'radar'
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta


def get_compare_annotator_data(data_folder, dtype, solvers=["gurobi", "mqlib", "pathfinder"]):
    line_labels = []
    spoke_labels = ["%Covered ", "%Used ", "Num. \ncontig ", "Num. break ", "Num. indel ", "Num. \ndiff ", "%Identity "]
    
    data = []
    
    for annotate in ["ga", "km", "mg"]:
        new_data = []
        for solver in solvers:
            file = f"{data_folder}/{solver}.{annotate}.{dtype}.avg.parsed.txt"
            with open(file, 'r') as f:
                lines = f.readlines()
            if solver == "pathfinder":
                for line in lines:
                    new_data.append([float(x) for x in line.strip().split(' ') if len(x)])
                    line_labels.append(solver_names[solver])
            else:
                for line in [lines[0], lines[-1]]:
                    new_data.append([float(x) for x in line.strip().split(' ')[1:]])
                    line_labels.append(f"{solver_names[solver]} {line.split(' ')[0]}")            
        new_data = np.array(new_data)
        data.append([annotate, new_data])
    all_data = np.array([d[1] for d in data])
    
    for i in [0, 1, 6]:
        all_data[:, :, i] = all_data[:, :, i] / 100
        spoke_labels[i] = f"{spoke_labels[i]}\n(100)"
        
    for i in [2, 3, 4, 5]:
        min_val = all_data[:, :, i].min()
        max_val = all_data[:, :, i].max()
        if max_val == 0:
            all_data[:, :, i] = 1
        else:
            all_data[:, :, i] = 1 - all_data[:, :, i] / max_val + min_val/max_val
        spoke_labels[i] = f"{spoke_labels[i]}\n({np.round(min_val,1)})"
        
    for i in range(len(data)):
        data[i][1] = all_data[i, :, :]
    return data, spoke_labels, line_labels


def get_compare_solver_data(data_folder, dtype, solvers=["gurobi", "mqlib", "pathfinder"]):
    line_labels = []
    spoke_labels = ["%Covered ", "%Used ", "Num. \ncontig ", "Num. break ", "Num. indel ", "Num. \ndiff ", "%Identity "]
    
    data = []
    
    for solver in solvers:
        new_data = []
        for annotate in ["ga", "km", "mg"]:
            
            file = f"{data_folder}/{solver}.{annotate}.{dtype}.avg.parsed.txt"
            with open(file, 'r') as f:
                lines = f.readlines()
                        
            new_data.append([float(x) for x in lines[-1].split(' ')[1:]])
            line_labels.append(f"{annotator_names[annotate]}")
                    
        new_data = np.array(new_data)    
        data.append([solver_names[solver], new_data])
        
    all_data = np.array([d[1] for d in data])
    
    for i in [0, 1, 6]:
        all_data[:, :, i] = all_data[:, :, i] / 100
        spoke_labels[i] = f"{spoke_labels[i]}\n(100)"
            
    for i in [2, 3, 4, 5]:
        min_val = all_data[:, :, i].min()
        max_val =  all_data[:, :, i].max()
        all_data[:, :, i] = 1 - all_data[:, :, i] / max_val + min_val/max_val
        spoke_labels[i] = f"{spoke_labels[i]}\n({np.round(min_val,1)})"
    for i in range(len(data)):
        data[i][1] = all_data[i, :, :]
        
    return data, spoke_labels, line_labels


def get_violin_plot_data(data_folder, dtype, solvers: list[str]=["gurobi", "mqlib", "pathfinder"]):
    line_labels = []
    new_data = []
    for solver in solvers:
        for annotate in ["ga", "km", "mg"]:
            to_add = None
            file = f"{data_folder}/{solver}.{annotate}.{dtype}.violin.txt"
            with open(file, 'r') as f:
                lines = f.readlines()

            in_data = False
            data = []        
            for line in lines:
                if not in_data:
                    # time_limit = int(line)
                    in_data = True
                elif line[0] == '-':
                    in_data = False
                    to_add = data
                    data = []
                    
                else:
                    data.append(float(line))
            
            new_data.append(to_add)
            line_labels.append(annotator_names[annotate])
            
                    
    return new_data, line_labels
            
        
def relocate_radar_labels(axs_row, theta, frame='polygon', r_label=None, fontsize=8, wrap_len=18):
    """
    Replace automatic polar/tick labels with manually placed labels at radius r_label.
    - axs_row: list of radar axes (all assumed to share the same theta)
    - theta: array of angles used by radar_factory
    - frame: 'polygon' or 'circle' (affects r_label default)
    - r_label: radial position for labels (in axis units, >1 pushes outside)
    - fontsize: label font size (pt)
    - wrap_len: length at which to insert a newline into long labels (None to disable)
    """
    if r_label is None:
        r_label = 1.06 if frame == 'circle' else 1.12

    for ax in axs_row:
        # collect existing labels (text) if any, else nothing
        old_labels = [t.get_text() for t in ax.get_xticklabels()]
        # fallback: if no old labels, try to get them from set_thetagrids metadata
        if not any(old_labels):
            # If you used ax.set_varlabels(...) earlier, those labels may be present in ax.xaxis properties;
            # best to pass spoke_labels to your plotting function; if not available we'll skip.
            old_labels = None

        # clear automatic labels
        ax.set_thetagrids([])

        labels_to_place = old_labels if old_labels else []
        # if labels empty, skip this axis (caller should ensure plotted axes used set_varlabels before)
        if not labels_to_place:
            continue

        for angle, label in zip(theta, labels_to_place):
            # optionally wrap long labels into two lines
            if wrap_len and len(label) > wrap_len:
                # split on whitespace nearest wrap_len
                parts = label.split()
                # naive wrap: insert newline after first few words
                acc = []
                s = ""
                for p in parts:
                    if len(s) + len(p) + 1 <= wrap_len or not s:
                        s = (s + " " + p).strip()
                    else:
                        acc.append(s)
                        s = p
                acc.append(s)
                label_text = "\n".join(acc) 
            else:
                label_text = label

            # compute readable rotation and alignment
            # For readability, don't rotate, keep centered under tick (use center)
            ha = 'center'
            va = 'center'

            ax.text(angle, r_label, label_text,
                    fontsize=fontsize,
                    horizontalalignment=ha,
                    verticalalignment=va,
                    rotation=0,
                    rotation_mode='anchor',
                    transform=ax.get_xaxis_transform())  # polar axes accept angle, r coords
            
            
            
def place_row_legend(fig, axs_row, handles, labels,
                     loc='upper right',
                     pad_x_frac=0.005,   # small x offset inside figure fraction
                     pad_y_frac=0.005,   # small y offset down from top of row
                     fontsize=8,
                     ):
    bbox = union_axes_bbox(axs_row)   # returns Bbox in fig coords (x0,y0,width,height)
    # compute anchor near top-right of the row bbox, nudged a little left/down so legend doesn't overlap axes border
    anchor_x = bbox.x0 + bbox.width / 2
    anchor_y = bbox.ymin - 0.030
    ncol = len(labels)

    # create a figure-level legend anchored to (anchor_x, anchor_y) in figure coords
    # bbox_transform=fig.transFigure makes bbox_to_anchor interpreted in figure coordinates
    fig.legend(handles, labels,
               loc=loc,
               bbox_to_anchor=(anchor_x, anchor_y),
               bbox_transform=fig.transFigure,
               fontsize=fontsize,
               ncol=ncol,
               frameon=False,
               borderaxespad=0.2,
               labelspacing=0.2)
    
    
def place_row_legend_inset(fig, axs_row, handles, labels,
                           legend_w_frac=0.17,  # width of inset axes as fraction of figure width
                           legend_h_frac=0.035,  # height (fraction) - tune to fit entries
                           pad_frac=0.01,       # gap between row bbox and inset
                           fontsize=8,
                           ncol=1):
    bbox = union_axes_bbox(axs_row)
    # compute inset axes coordinates to the right of the row bbox
    inset_left = min(bbox.x0 + bbox.width + pad_frac, 0.995 - legend_w_frac)  # clamp so it doesn't pass figure right edge
    inset_bottom = bbox.ymax + legend_h_frac  # align top of inset with row top
    # clamp bottom to stay inside figure
    inset_bottom = max(inset_bottom, 0.01)
    inset = fig.add_axes([inset_left, inset_bottom, legend_w_frac, legend_h_frac])
    inset.axis('off')
    # put legend inside this inset axes
    inset.legend(handles, labels, loc='upper left', fontsize=fontsize, frameon=False, ncol=ncol, labelspacing=0.2)
    return inset