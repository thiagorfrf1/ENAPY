# -*- coding: UTF-8 -*-
import time

import pandas as pd
from rpy2.robjects import r, pandas2ri, Formula
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
rENA = importr('rENA')

rsdata = pd.read_csv("rsdata.csv")

pandas2ri.activate()  # makes some conversions automatic

print("Keys")
pr = ro.r("print")
colMeans = ro.r("colMeans")
rescale = ro.r("scales::rescale")
rmax = ro.r("max")

print(rENA._rpy2r.keys())

# load your file

units = rsdata[['Condition','UserName']]
conversation = rsdata[['Condition','GroupName']]
codes = rsdata[['Data','Technical.Constraints','Performance.Parameters','Client.and.Consultant.Requests','Design.Reasoning','Collaboration']]
meta = rsdata[["CONFIDENCE.Change","CONFIDENCE.Pre","CONFIDENCE.Post","C.Change"]]

accum = rENA.ena_accumulate_data(units, conversation, codes, meta)

set = rENA.ena_make_set(
  enadata=accum
)

points = (set.rx2('points'))

df_points1 = points.loc[points['Condition'] == "FirstGame"]
df_points2 = points.loc[points['Condition'] == "SecondGame"]
### Subset rotated points for the first condition
first_game_points = df_points1.drop(columns=['ENA_UNIT','Condition','UserName', 'CONFIDENCE.Change', 'CONFIDENCE.Pre', 'CONFIDENCE.Post', "C.Change"])
### Subset rotated points for the second condition
second_game_points = df_points2.drop(columns=['ENA_UNIT','Condition','UserName',  'CONFIDENCE.Change', 'CONFIDENCE.Pre', 'CONFIDENCE.Post', "C.Change"])




plot = rENA.ena_plot(set, scale_to ="network", title ="Groups of Units")
plot = rENA.ena_plot_points(plot, points = first_game_points, confidence_interval ="box", colors = ("blue"))
plot = rENA.ena_plot_points(plot, points = second_game_points, confidence_interval ="box", colors = ("red"))
pr(plot)



plot = rENA.ena_plot(set, scale_to = [-1, 0, 1], title ="Groups and Means")
plot = rENA.ena_plot_points(plot, points = first_game_points, confidence_interval ="box", colors = ("blue"))
plot = rENA.ena_plot_points(plot, points = second_game_points, confidence_interval ="box", colors = ("red"))
plot = rENA.ena_plot_group(plot, first_game_points, colors = ("red"), confidence_interval ="box")
plot = rENA.ena_plot_group(plot, second_game_points, colors =("blue"), confidence_interval ="box")
pr(plot)

line_weights = (set.rx2('line.weights'))

### Subset lineweights for SecondGame and Calculate the colMeans
first_game_lineweights = line_weights.loc[points['Condition'] == "FirstGame"]
second_game_lineweights = line_weights.loc[points['Condition'] == "SecondGame"]

first_game_lineweights = first_game_lineweights.drop(columns=['ENA_UNIT','Condition','UserName', 'CONFIDENCE.Change', 'CONFIDENCE.Pre', 'CONFIDENCE.Post', "C.Change"])
### Subset rotated points for the second condition
second_game_lineweights = second_game_lineweights.drop(columns=['ENA_UNIT','Condition','UserName',  'CONFIDENCE.Change', 'CONFIDENCE.Pre', 'CONFIDENCE.Post', "C.Change"])



first_game_mean = colMeans(first_game_lineweights)
second_game_mean = colMeans(second_game_lineweights)

### Subtract the two sets of means, resulting in a vector with negative values
### indicatinag a stronger connection with the SecondGame, and positive values
### a stronger FirstGame connection
subtracted_mean = first_game_mean - second_game_mean

#> [1]  0.0769062510 -0.0409819566  0.0204737523  0.0008708646  0.0373163091
#Plot subtracted network only
plot_first = rENA.ena_plot(set, title ="FirstGame")
plot_first = rENA.ena_plot_network(plot_first, network = first_game_mean)
pr(plot_first)


plot_second = rENA.ena_plot(set, title ="SecondGame")
plot_second = rENA.ena_plot_network(plot_second, network = second_game_mean, colors = ("blue"))
pr(plot_second)




plot_sub = rENA.ena_plot(set, title ="Subtracted")
plot_sub = rENA.ena_plot_network(plot_sub, network = subtracted_mean)
pr(plot_sub)


rotation = (set.rx2('rotation'))
nodes = rotation.rx2('nodes')

# Scale the nodes to match that of the network, for better viewing
first_game_points_max = rmax(first_game_points)
second_game_points_max = rmax(second_game_points)
if(first_game_points_max> second_game_points_max):
  point_max = first_game_points_max
else:
  point_max = second_game_points_max

nodes = nodes.drop(columns=['code'])
print("nodes")
print(nodes)
max_nodes = rmax(nodes)
print("max_nodes")
print(max_nodes)
print("point_max")
print(point_max)

with localconverter(ro.default_converter + pandas2ri.converter):
  first_game_points = ro.conversion.py2rpy(first_game_points)
print(type(first_game_points))

first_game_scaled = first_game_points
second_game_scaled = second_game_points

plot = rENA.ena_plot(set, title ="Plot with Units and Network", font_family ="Times")
plot = rENA.ena_plot_points(plot, points = first_game_scaled, colors = ("red"))
plot = rENA.ena_plot_points(plot, points = second_game_scaled, colors = ("blue"))
plot = rENA.ena_plot_group(plot, point = first_game_scaled, colors =("red"), confidence_interval ="box")
plot = rENA.ena_plot_group(plot, point = second_game_scaled, colors =("blue"), confidence_interval ="box")
plot = rENA.ena_plot_network(plot, network = subtracted_mean)
pr(plot)

time.sleep(10)