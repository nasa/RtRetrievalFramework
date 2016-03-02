from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
from .decorators import call_data_pairs

import numpy
from matplotlib.pyplot import *

from .routines_base import PlotRoutinesBase
from .status_routines import StatusRoutines

class CustomPlotRoutines(PlotRoutinesBase):
    def _plot_counts_histogram(self, count_name, obj_names, data_values, **kwargs):
        ax = self.new_axis(**kwargs)


        names_with_mean = []
        for d_name, d_iters in zip(obj_names, data_values):
            names_with_mean.append("%s (mean %.2f)" % (d_name, numpy.mean(d_iters)))

        plot_data = numpy.transpose(numpy.array(data_values))
        ax.hist(plot_data, align='left', histtype='bar', label=names_with_mean)

        ind = numpy.arange(numpy.max(plot_data))
        ax.set_xticks(ind)
        ax.set_xticklabels([str(i) for i in ind])

        ax.legend(loc='best')

        ax.set_title("%s histogram" % count_name)
        ax.set_xlabel("Number of %s" % count_name)
        ax.set_ylabel("Count")

        return ax

    def plot_iter_histogram(self, obj_names, iterations, **kwargs):
        """Plot histogram of number of iterations for multiple sets of data"""
        return self._plot_counts_histogram("iterations", obj_names, iterations, **kwargs)

    def plot_div_histogram(self, obj_names, diverging_steps, **kwargs):
        """Plot histogram of number of divergence steps for multiple sets of data"""
        return self._plot_counts_histogram("divergences", obj_names, diverging_steps, **kwargs)

    def plot_statevector_diff_trend(self, obj_names, sounding_id, RetrievedStateVector__state_vector_names, RetrievedStateVector__state_vector_result, RetrievedStateVector__state_vector_aposteriori_uncert, **kwargs):
        """Plot the trends for how state vector elements are differing between two sets of runs
        """
        # Call helper where we can use call_data_pairs and not loose inspection of analysis routine into our behavior
        return self._plot_statevector_diff_trend_helper(obj_names, sounding_id, RetrievedStateVector__state_vector_names, RetrievedStateVector__state_vector_result, RetrievedStateVector__state_vector_aposteriori_uncert, **kwargs)

    @call_data_pairs(list(range(1,6)))
    def _plot_statevector_diff_trend_helper(self, obj_names, sounding_id, state_vector_names, state_vector_result, state_vector_aposteriori_uncert, **kwargs):
        rtol = kwargs.get("rtol", 1e-5)
        atol = kwargs.get("atol", 1e-8)
        num_plot_elems = kwargs.get("num_plot_elems", 30)
        absolute = kwargs.get("absolute", False)

        exceed_pos_dict = {}
        exceed_neg_dict = {}
        mean_pos_dict = {}
        mean_neg_dict = {}
        all_sv_names = []
        seen_sv = {}
        for snd_idx in range(state_vector_names[0].shape[0]):
            if numpy.any(state_vector_names[0][snd_idx,:] != state_vector_names[1][snd_idx,:]):
                raise ValueError("State vector names posistions not the same for sounding %s" % sounding_id[0][snd_idx])

            for sv_idx, sv_name in enumerate(state_vector_names[0][snd_idx,:]):
                val_1 = state_vector_result[0][snd_idx,sv_idx]
                val_2 = state_vector_result[1][snd_idx,sv_idx]

                if absolute:
                    vals_diff = abs(val_1 - val_2)
                else:
                    vals_diff = val_1 - val_2

                #cmp_mag = numpy.mean([val_1,val_2])
                cmp_mag =  numpy.mean([state_vector_aposteriori_uncert[0][snd_idx,sv_idx],
                                       state_vector_aposteriori_uncert[1][snd_idx,sv_idx]])
                vals_tol = atol + rtol * numpy.abs(cmp_mag)

                if sv_name not in seen_sv:
                    seen_sv[sv_name] = 1
                    all_sv_names.append(sv_name)

                if(vals_diff > vals_tol):
                    exceed_pos_dict[sv_name] = exceed_pos_dict.get(sv_name,0) + 1
                    mean_pos_dict[sv_name] = mean_pos_dict.get(sv_name,0) + vals_diff

                if(vals_diff < -vals_tol):
                    exceed_neg_dict[sv_name] = exceed_neg_dict.get(sv_name,0) - 1
                    mean_neg_dict[sv_name] = mean_neg_dict.get(sv_name,0) + vals_diff

                    #w_fine = set(list(numpy.where(-vals_tol < vals_diff)[0])).intersection(set(list(numpy.where(vals_diff < vals_tol)[0])))
            #print w_fine, '\n', [ vals_diff[idx] for idx in w_fine ]


        names = []
        exceed_pos = numpy.zeros(len(all_sv_names), dtype=int)
        exceed_neg = numpy.zeros(len(all_sv_names), dtype=int)
        mean_pos = numpy.zeros(len(all_sv_names), dtype=float)
        mean_neg = numpy.zeros(len(all_sv_names), dtype=float)
        for sv_idx, sv_name in enumerate(all_sv_names):
            exceed_pos[sv_idx] = exceed_pos_dict.get(sv_name, 0)
            exceed_neg[sv_idx] = exceed_neg_dict.get(sv_name, 0)

            if exceed_pos[sv_idx] > 0:
                mean_pos[sv_idx] = old_div(mean_pos_dict.get(sv_name,0), exceed_pos[sv_idx])

            if exceed_neg[sv_idx] > 0:
                mean_neg[sv_idx] = old_div(mean_neg_dict.get(sv_name,0), exceed_neg[sv_idx])

            # Make names shorter
            curr_nm = sv_name
            curr_nm = curr_nm.replace(" for Press Lvl","")
            curr_nm = curr_nm.replace(" VMR","")
            curr_nm = curr_nm.replace(" Log(extinction)","")
            curr_nm = curr_nm.replace("Temperature","Temp")
            curr_nm = curr_nm.replace("Surface Pressure","PSURF")
            curr_nm = curr_nm.replace(" factor","")
            curr_nm = curr_nm.replace(" (Kelvin)","")
            curr_nm = curr_nm.replace(" (Pascals)","")
            curr_nm = curr_nm.replace(" Parm","")
            curr_nm = curr_nm.replace("Ground Lambertian ","")
            curr_nm = curr_nm.replace(" Offset","")
            curr_nm = curr_nm.replace("GOSAT Inst ","")

            names.append(curr_nm)

        # Divide the counts into smaller sets for easier plotting
        all_indexes = list(range(len(names)))
        part_indexes_list = []
        while len(all_indexes) > 0:
            part_indexes_list.append(all_indexes[:num_plot_elems])
            all_indexes = all_indexes[num_plot_elems:]

        yrange = min(exceed_neg), max(exceed_pos)

        res = []
        for part_num, part_idx in enumerate(part_indexes_list):
            part_names = [ names[idx] for idx in part_idx ]
            pos_part = numpy.array([ exceed_pos[idx] for idx in part_idx ])
            neg_part = numpy.array([ exceed_neg[idx] for idx in part_idx ])

            ax = self.new_axis(**kwargs)
            ax.get_figure().subplots_adjust(left=0.10, right=0.98, bottom=0.30) # make room for names
            
            ax.xaxis.set_major_locator(matplotlib.ticker.IndexLocator(1,1))
            ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(part_names))

            for label in ax.get_xticklabels():
                label.set_ha("right")
                label.set_rotation(90)

            ax.bar(numpy.arange(len(part_names)), pos_part)
            ax.bar(numpy.arange(len(part_names)), neg_part)

            #for mean, x, y in zip(mean_pos, numpy.arange(len(part_names)), pos_part):
            #    print '%e'%mean, x, y
            #    ax.annotate(' %.1e' % mean, xy=(x,y), xycoords='data', rotation='vertical')

            # So all figures have the same y range
            ax.set_ylim(yrange)

            # So last label tick shows up
            ax.set_xlim([0, len(part_idx)+0.1])

            ax.set_title("Plot %d of Num (%s - %s) Diffs Exceeding\ntol = (%.2e + %2.e * apost_uncert)" % (part_num+1, obj_names[0], obj_names[1], atol, rtol))

            if absolute:
                ax.set_ylabel("# > 0 := diff > tol")
            else:
                ax.set_ylabel("# < 0 := diff < -tol ; # > 0 := diff > tol")

            res.append(ax)

        return res
