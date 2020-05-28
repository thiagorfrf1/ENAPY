import pandas as pd
from rpy2.robjects import r, pandas2ri, Formula
from rpy2.robjects.packages import importr
import rpy2.robjects as ro


rsdata = pd.read_csv("example.csv")




rENA = importr("rENA")


pandas2ri.activate()  # makes some conversions automatic

units = rsdata[['Condition','UserName']]
print(units.head())
conversation = rsdata[['Condition','GroupName']]
print(conversation.head())
codes = rsdata[['Data','Technical.Constraints','Performance.Parameters','Client.and.Consultant.Requests','Design.Reasoning','Collaboration']]
print(codes.head())
meta = rsdata[["CONFIDENCE.Change","CONFIDENCE.Pre","CONFIDENCE.Post","C.Change"]]
print("meta")
print(meta.head())

accum = rENA.ena_accumulate_data(units, conversation, codes, meta)

set = rENA.ena_make_set(
  enadata=accum
)

points = (set.rx2('points'))

print("points")
print(points)

df_points1 = points.loc[points['Condition'] == "FirstGame"]
df_points2 = points.loc[points['Condition'] == "SecondGame"]

### Subset rotated points for the first condition
first_game_points = df_points1.drop(columns=['ENA_UNIT','Condition','UserName', 'CONFIDENCE.Change', 'CONFIDENCE.Pre', 'CONFIDENCE.Post', "C.Change"])
### Subset rotated points for the second condition
second_game_points = df_points2.drop(columns=['ENA_UNIT','Condition','UserName',  'CONFIDENCE.Change', 'CONFIDENCE.Pre', 'CONFIDENCE.Post', "C.Change"])


plot = rENA.ena_plot(set, scale_to = "network", title = "Groups of Units")

plot = rENA.ena_plot_points(plot, points = first_game_points, confidence_interval = "box", colors = ("blue"))
plot = rENA.ena_plot_points(plot, points = second_game_points, confidence_interval = "box", colors = ("red"))

ro.r.plot("plot")