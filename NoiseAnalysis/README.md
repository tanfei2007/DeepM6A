




noise level 

fdr = fp/(tp + fp)

when smart-seq has different fdr (error rate ), we do simulation experiments to demonstate the 
the evolution of performance of the proposed method and the patterns discovered. 


## Patterns: TGCA + XOR

fdr = 0

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 1.0|1.0 |1.0 |1.0 | |



fdr = 0.5

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.9999845| 0.997116710543|1.0 | 0.99425| |



fdr = 0.8

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.999902375| 0.965071151358 | 1.0  | 0.9325 | |



fdr = 0.9

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
| 0.99952825375| 0.943475964078  |  1.0 | 0.893  | |


## Patterns: TGCA

fdr = 0.1

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|  | |  | | |

fdr = 0.2

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|  | |  | | |



## Patterns: XOR

fdr = 0.1

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|  | |  | | |

fdr = 0.2

| auc@roc| f1 score| precision |  recall | auc@pr|
|--|--| -- | --| -- |
|  | |  | | |
