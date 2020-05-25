# -*- coding: UTF-8 -*-
from stringr import accumulate_data
from stringr import cohens_d
from stringr import connection_matrix
from stringr import ena
from stringr import ena_accumulate_data
from stringr import ena_accumulate_data_file
from stringr import ena_conversations
from stringr import ena_correlations
from stringr import ena_generate
from stringr import ena_get_stats
from stringr import ena_group
from stringr import ena_make_set
from stringr import ena_optimize_set
from stringr import ena_plot
from stringr import ena_plot_group
from stringr import ena_plot_network
from stringr import ena_plot_points
from stringr import ena_plot_subtraction
from stringr import ena_plot_trajectory
from stringr import ena_plotter
from stringr import ena_rotate_by_mean
from stringr import ena_set
from stringr import ena_set_creator
from stringr import ena_svd
from stringr import ena_unit_group
from stringr import ena_unit_metadata
from stringr import ena_update_set
from stringr import ena_writeup
from stringr import ENAdata
from stringr import ENAplot
from stringr import ENArotation_set
from stringr import ENAset
from stringr import lws_positions_sq
from stringr import namesToAdjacencyKey
from stringr import ocpu_optimization
from stringr import RcppExports
from stringr import rENA
from stringr import RS_data
from stringr import utils
from stringr import utils_classes
from stringr import utils_matrix
from stringr import utils_plot
from stringr import zzz

from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP, SignatureTranslatedAnonymousPackage


string = str(cohens_d.cohens_d)+\
         str(connection_matrix.connection_matrix)+\
         str(ena+ena_accumulate_data.ena_accumulate_data)+\
         str(ena_accumulate_data_file.ena_accumulate_data_file)+\
         str(ena_conversations.ena_conversations)+\
         str(ena_correlations.ena_correlations)+\
         str(ena_generate.ena_generate)+\
         str(ena_get_stats.ena_get_stats)+\
         str(ena_group.ena_group)+\
         str(ena_make_set.ena_make_set)+\
         str(ena_optimize_set.ena_optimize_set)+\
         str(ena_plot.ena_plot)+\
         str(ena_plot_group.ena_plot_group)+\
         str(ena_plot_network.ena_plot_network)+\
         str(ena_plot_points.ena_plot_points)+\
         str(ena_plot_subtraction.ena_plot_subtraction)+\
         str(ena_plot_trajectory.ena_plot_trajectory)+\
         str(ena_plotter.ena_plotter)+\
         str(ena_rotate_by_mean.ena_rotate_by_mean)+\
         str(ena_set.ena_set)+\
         str(ena_set_creator.ena_set_creator)+\
         str(ena_svd.ena_svd)+\
         str(ena_unit_group.ena_unit_group)+\
         str(ena_unit_metadata.ena_unit_metadata)+\
         str(ena_update_set.ena_update_set)+\
         str(ena_writeup.ena_writeup)+\
         str(lws_positions_sq.lws_positions_sq)+\
         str(namesToAdjacencyKey.namesToAdjacencyKey)+\
         str(ocpu_optimization.ocpu_optimization)+\
         str(zzz.zzz)


#ENAdata = STAP(ENAdata.ENAdata, "ENAdata")
#ENAplot = STAP(ENAplot.ENAplot, "ENAplot")
#ENArotation_set = STAP(ENArotation_set.ENArotation_set, "ENArotation_set")
#ENAset = STAP(ENAset.ENAset, "ENAset")
#rENA = STAP(rENA.rENA, "rENA")
#utils = STAP(utils.utils, "utils")
#accumulate_data = STAP(accumulate_data.accumulate_data, "accumulate_data")


print("Keys")
stringr_c = STAP(string, "stringr_c")
print(stringr_c._rpy2r.keys())


#string_c.ena_plot(enaset, title = "ENA Plot", dimension.labels = c("", "")
#  font.size = 10, font.color = "#000000", font.family = c("Arial",
#  "Courier New", "Times New Roman"), scale.to = "network", ...)