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

from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP



#accumulate_data = STAP(accumulate_data.accumulate_data, "accumulate_data")

cohens_d = STAP(cohens_d.cohens_d, "cohens_d")
connection_matrix = STAP(connection_matrix.connection_matrix, "connection_matrix")
ena = STAP(ena.ena, "ena")
ena_accumulate_data = STAP(ena_accumulate_data.ena_accumulate_data, "ena_accumulate_data")
ena_accumulate_data_file = STAP(ena_accumulate_data_file.ena_accumulate_data_file, "ena_accumulate_data_file")
ena_conversations = STAP(ena_conversations.ena_conversations, "ena_conversations")
ena_correlations = STAP(ena_correlations.ena_correlations, "ena_correlations")
ena_generate = STAP(ena_generate.ena_generate, "ena_generate")
ena_get_stats = STAP(ena_get_stats.ena_get_stats, "ena_get_stats")
ena_group = STAP(ena_group.ena_group, "ena_group")
ena_make_set = STAP(ena_make_set.ena_make_set, "ena_make_set")
ena_optimize_set = STAP(ena_optimize_set.ena_optimize_set, "ena_optimize_set")
ena_plot = STAP(ena_plot.ena_plot, "ena_plot")
ena_plot_group = STAP(ena_plot_group.ena_plot_group, "ena_plot_group")
ena_plot_network = STAP(ena_plot_network.ena_plot_network, "ena_plot_network")
ena_plot_points = STAP(ena_plot_points.ena_plot_points, "ena_plot_points")
ena_plot_subtraction = STAP(ena_plot_subtraction.ena_plot_subtraction, "ena_plot_subtraction")
ena_plot_trajectory = STAP(ena_plot_trajectory.ena_plot_trajectory, "ena_plot_trajectory")
ena_plotter = STAP(ena_plotter.ena_plotter, "ena_plotter")
ena_rotate_by_mean = STAP(ena_rotate_by_mean.ena_rotate_by_mean, "ena_rotate_by_mean")
ena_set = STAP(ena_set.ena_set, "ena_set")
ena_set_creator = STAP(ena_set_creator.ena_set_creator, "ena_set_creator")
ena_svd = STAP(ena_svd.ena_svd, "ena_svd")
ena_unit_group = STAP(ena_unit_group.ena_unit_group, "ena_unit_group")
ena_unit_metadata = STAP(ena_unit_metadata.ena_unit_metadata, "ena_unit_metadata")
ena_update_set = STAP(ena_update_set.ena_update_set, "ena_update_set")
ena_writeup = STAP(ena_writeup.ena_writeup, "ena_writeup")
ENAdata = STAP(ENAdata.ENAdata, "ENAdata")
ENAplot = STAP(ENAplot.ENAplot, "ENAplot")
ENArotation_set = STAP(ENArotation_set.ENArotation_set, "ENArotation_set")
ENAset = STAP(ENAset.ENAset, "ENAset")
lws_positions_sq = STAP(lws_positions_sq.lws_positions_sq, "lws_positions_sq")
namesToAdjacencyKey = STAP(namesToAdjacencyKey.namesToAdjacencyKey, "namesToAdjacencyKey")
ocpu_optimization = STAP(ocpu_optimization.ocpu_optimization, "ocpu_optimization")
RcppExports = STAP(RcppExports.RcppExports, "RcppExports")
rENA = STAP(rENA.rENA, "rENA")
RS_data = STAP(RS_data.RS_data, "RS_data")
utils = STAP(utils.utils, "utils")
utils_classes = STAP(utils_classes.utils_classes, "utils_classes")
utils_matrix = STAP(utils_matrix.utils_matrix, "utils_matrix")
utils_plot = STAP(utils_plot.utils_plot, "utils_plot")
zzz = STAP(zzz.zzz, "zzz")


print("Keys")
print(ena_set._rpy2r.keys())

