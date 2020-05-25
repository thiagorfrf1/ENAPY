# -*- coding: UTF-8 -*-
import rpy2.robjects as robjects

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
from rpy2.robjects import IntVector, Formula
import pandas as pd

pandas2ri.activate()

teste = """
library(rENA)
data(RS.data)


units = RS.data[,c("Condition","UserName")]
head(units)
#>   Condition    UserName
#> 1 FirstGame    steven z
#> 2 FirstGame     akash v
#> 3 FirstGame alexander b
#> 4 FirstGame   brandon l
#> 5 FirstGame   brandon l
#> 6 FirstGame christian x

conversation = RS.data[,c("Condition","GroupName")]
head(conversation)
#>   Condition GroupName
#> 1 FirstGame  Electric
#> 2 FirstGame  Electric
#> 3 FirstGame  Electric
#> 4 FirstGame  Electric
#> 5 FirstGame  Electric
#> 6 FirstGame  Electric

codeCols = c(
  'Data','Technical.Constraints','Performance.Parameters',
  'Client.and.Consultant.Requests','Design.Reasoning','Collaboration'
)
codes = RS.data[,codeCols]
head(codes)
#>   Data Technical.Constraints Performance.Parameters
#> 1    0                     0                      0
#> 2    0                     0                      0
#> 3    0                     0                      0
#> 4    0                     0                      0
#> 5    0                     0                      0
#> 6    0                     0                      0
#>   Client.and.Consultant.Requests Design.Reasoning Collaboration
#> 1                              0                0             0
#> 2                              0                0             0
#> 3                              0                0             0
#> 4                              0                0             0
#> 5                              0                0             0
#> 6                              0                0             0

# optional
meta = RS.data[,c("CONFIDENCE.Change",
                  "CONFIDENCE.Pre","CONFIDENCE.Post","C.Change")]
head(meta)
#>   CONFIDENCE.Change CONFIDENCE.Pre CONFIDENCE.Post   C.Change
#> 1                 1              7               8 Pos.Change
#> 2                 2              6               8 Pos.Change
#> 3                 1              5               7 Pos.Change
#> 4                 1              5               6 Pos.Change
#> 5                 1              5               6 Pos.Change
#> 6                 0              4               4 Neg.Change





accum = ena.accumulate.data(
  units = units,
  conversation = conversation,
  codes = codes,
  metadata = meta,
  window.size.back = 4
)

### adjacency.vectors: Each Unit's Co-Occurrence Accumulation
#head(accum$adjacency.vectors)

### adjacency.matrix: Columns representing co-occurred 
### codes in the adjacency.vector
#head(accum$adjacency.matrix)










set = ena.make.set(
  enadata = accum
)

### The location in space for each unit.  Units are rows, columns are each 
### dimension in the high-dimensional space.
#head(set$points.rotated)

### The positiona of each code in the high-dimensional space
#head(set$node.positions)

### The weight of each connection. Units are rows, columns the co-occurrence
#head(set$line.weights)










### Subset rotated points for the first condition
first.game.points = as.matrix(set$points$Condition$FirstGame)

### Subset rotated points for the second condition
second.game.points = as.matrix(set$points$Condition$SecondGame)

plot = ena.plot(set, scale.to = "network", title = "Groups of Units")
plot = ena.plot.points(plot, points = first.game.points, confidence.interval = "box", colors = c("red"))
plot = ena.plot.points(plot, points = second.game.points, confidence.interval = "box", colors = c("blue"))
#plot$plot







### Using the same plot object above, we will be able to plot the means 
### alongside their corresponding units.
plot = ena.plot(set, scale.to = list(x=-1:1, y=-1:1), title = "Groups and Means")
plot = ena.plot.points(plot, points = first.game.points, 
                       confidence.interval = "box", colors = c("red"))
plot = ena.plot.points(plot, points = second.game.points, 
                       confidence.interval = "box", colors = c("blue"))
plot = ena.plot.group(plot, point = first.game.points, 
                      colors =c("red"), confidence.interval = "box")
plot = ena.plot.group(plot, point = second.game.points, 
                      colors =c("blue"), confidence.interval = "box")
#plot$plot



### Subset lineweights for FirstGame and Calculate the colMeans
first.game.lineweights = as.matrix(set$line.weights$Condition$FirstGame)

### Subset lineweights for SecondGame and Calculate the colMeans
second.game.lineweights = as.matrix(set$line.weights$Condition$SecondGame)
first.game.mean = as.vector(colMeans(first.game.lineweights))
second.game.mean = as.vector(colMeans(second.game.lineweights))

### Subtract the two sets of means, resulting in a vector with negative values
### indicatinag a stronger connection with the SecondGame, and positive values
### a stronger FirstGame connection
subtracted.mean = first.game.mean - second.game.mean

# View the first 5 elements to see the substraction
head(first.game.mean, 5)
#> [1] 0.40328303 0.27885635 0.34790907 0.08888453 0.10970647
head(second.game.mean, 5)
#> [1] 0.32637678 0.31983831 0.32743531 0.08801366 0.07239016
head(subtracted.mean, 5)
#> [1]  0.0769062510 -0.0409819566  0.0204737523  0.0008708646  0.0373163091






### Subset rotated points for the first condition
first.game.points = as.matrix(set$points$Condition$FirstGame)

### Subset rotated points for the second condition
second.game.points = as.matrix(set$points$Condition$SecondGame)

plot = ena.plot(set, scale.to = "network", title = "Groups of Units")
plot = ena.plot.points(plot, points = first.game.points, confidence.interval = "box", colors = c("red"))
plot = ena.plot.points(plot, points = second.game.points, confidence.interval = "box", colors = c("blue"))
plot$plot

"""


stringfinal = teste

stringr_c = STAP(stringfinal, "stringr_c")

print("Keys")

print(stringr_c._rpy2r.keys())
print(stringr_c.plot.find("plot" ,False))
print(type(stringr_c.plot))

print(stringr_c.codeCols)

stringr_c.plot.find("plot",False)
