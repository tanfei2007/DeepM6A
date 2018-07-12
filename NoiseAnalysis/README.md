




noise level 

fdr = fp/(tp + fp)

when smart-seq has different fdr (error rate ), we do simulation experiments to demonstate the 
the evolution of performance of the proposed method and the patterns discovered. 


## Patterns: TGCA + XOR

fdr = 0

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 1.0|1.0 |1.0 |1.0 | |


fdr = 0.1

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.999999625|0.998377636341 |0.996760528283 |1.0 | |


fdr = 0.2

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|0.99999775 | 0.998377636341 | 0.996760528283 |1.0 | |


fdr = 0.5

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.999986125| 0.99787950605 | 0.995767986059|1.0 | |



fdr = 0.8

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|0.9997775 | 0.980272025487 |  0.961307378034 |  1.0 | |



fdr = 0.9

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|0.99958375 | 0.945626477541  |0.896860986547   |1.0   | |


## Patterns: TGCA

fdr = 0

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 1.0|1.0 |1.0 |1.0 | |


fdr = 0.1

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.999998625 |0.998377636341  | 0.996760528283 |1.0 | |

fdr = 0.2

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|0.999998375  |0.998377636341 | 0.996760528283 |1.0 | |


fdr = 0.5

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.999969125  | 0.997008973081 | 0.994035785288  |1.0 | |

fdr = 0.8

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.999443625  |  0.980872976949 | 0.962463907603   |1.0 | |


fdr = 0.9

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.99956975  |  0.948091964921 |  0.901306894998  |1.0 | |


## Patterns: XOR

fdr = 0.1

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|  | |  | | |

fdr = 0.2

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|  | |  | | |
