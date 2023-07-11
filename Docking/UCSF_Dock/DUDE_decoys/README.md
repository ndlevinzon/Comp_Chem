# What are decoys?
Or, more importantly, why do we care about decoys and keep hearing about decoys all the time?

"Decoys" is essentially a codename for a way of evaluating how well your docking program has done on a target (or a set of targets). Decoys refers to a set of molecules that (probably) won't bind to your target. 
## Terms:
* Ligands: A set of known ligands that bind to your protein target. Often taken from papers or a database like ChEMBL
* Known Decoys/ Known non-binders: A set of molecules that have been tested against your protein target and found not to bind. ChEMBL or the literature is also a source of these.
* Property-Matched Decoys: A set of molecules, typically from ZINC, that look like your ligands in chemical and physical property, but are not similar to the ligands by Tanimoto of a 2D fingerprint. Making these is explained here Automated_Database_Preparation#Automatic_Decoy_Generation and more information can be found in the Huang et al [1] or Verdonk et al [2] papers.
    * Random Decoys: A set of random molecules, usually chosen from ZINC, most of which won't be binders simply by chance.

Now, once you have these various sets, you can examine the enrichment of various sets, usually by looking at a ROC curve, log ROC curve, or the LogAUC of one set over the others. The most common usage is ligands over property-matched decoys. If your target does well at this, the general attitude is that you will do well at a prospective virtual screen. If you have many known decoys or known non-binders you can examine the enrichment of those over ligands, which also tends to indicate how well you are doing. Sometimes it is also illustrative to test the enrichment of ligands over random decoys (usually something like the leadlike ZINC subset, trimmed at 60% Tanimoto overlap [1]. A final thing to examine if you have known decoys is to examine their enrichment over random decoys, expecting random performance.

References:
``
[1] Huang N, Shoichet BK, Irwin JJ. Benchmarking sets for molecular docking. J Med Chem. 2006 Nov 16; 49(23):6789-801.
[2] Marcel L. Verdonk*, Valerio Berdini, Michael J. Hartshorn, Wijnand T. M. Mooij, Christopher W. Murray, Richard D. Taylor, and Paul Watson. J. Chem. Inf. Comput. Sci., 2004, 44 (3), pp 793–806.
``
# LogAUC
LogAUC is a metric to evaluate virtual screening performance that has many of the same advantages as area under the curve (AUC), but is based on a plot where the x-axis is semilog in order to focus on early enrichment.
## Motivation
When we look at virtual screening performance, we plot an ROC curve (or enrichment curve) with a base 10 semilog x-axis, because this has the advantage of focusing the graph on "early enrichment", where molecules are most likely to be selected for further testing. If we had instead plotted the curve with the usual linear x-axis, then the area under the curve (AUC) is a well-regarded metric to summarize the overall performance of a virtual screening campaign as a single number1. While AUC can be formulated alternate ways2,3, it can be mechanically constructed by simply integrating under the curve, and interpreted as the fraction of the area under the curve over the area under the best possible ROC curve. It just happens that in a linear ROC plot, the AUC of the best possible curve is the entire unit square, with an area of 1. By analogy, in our typical semilog plots, we can construct the same fraction of the area under the log curve, over the area under the perfect log curve, and define that fraction as the logAUC. The lone nuisance is that the area under the log curve is infinite in general. However, if we are practical and limit our focus to a region of log space that we can actually measure, say above a certain threshold <math>\lambda</math>, then the perfect log area is finite.
## Definition
Formally, we define <math>logAUC_\lambda</math>, where the log area computations run from <math>\lambda</math> to 1.0, and we typically refer to <math>logAUC_{0.001}</math> as simply <math>logAUC</math>, where the area is integrated from 0.1 percent (0.001) to 100 percent (1.0) of decoys found. For integrating the area under the curve, we use the trapezoidal rule as follows:

<math>LogAUC_\lambda=\frac{\displaystyle \sum_{i}^{where~x_i\ge\lambda} (\log_{10} x_{i+1} - \log_{10} x_i)(\frac{y_{i+1}+y_i}{2})}{\log_{10}\frac{1}{\lambda}}</math>
## Discussion
From similar reasoning based on semilog ROC plots, Clark and Webster-Clark construct the pROC AUC metric2, which is similar to the numerator of logAUC except that the integration is done over horizontal bars instead of vertical trapezoids. The advantage of constructing logAUC as a fraction over the ideal area is that the choice of base for the logarithm is irrelevant, because changing base simply results in a constant that cancels between numerator and denominator. Also, by explicitly defining the area of interest using λ and integrating vertically, we are able to avoid the singularity at <math>x_i=0</math> encountered in pROC. More importantly, the fixed integration area means we can more directly compare <math>logAUC_\lambda</math> values across databases of different sizes and across targets with different ratios of actives to inactives. The final advantage of logAUC is that if you are used to looking at semilog ROC plots plotted from λ to 1, and understand that logAUC is just the percentage of the total area below the curve, then you can at some point gain the same intuitive feel as AUC has for linear ROC plots. In a semilog ROC plot the random line occupies only a sliver of the total area, and indeed its logAUC is just 14.462%. In order to more easily compare a given logAUC to this random value, we instead report the “adjusted logAUC” as the calculated value minus 14.462%, so that positive values mean overall enrichments better than random.

<math>Adjusted~LogAUC=LogAUC_{0.001}-0.14462</math>
References:
```
[1] Nicholls, A., What do we know and when do we know it? J Comput Aided Mol Des 2008, 22, (3-4), 239-55.
[2] Clark, R. D.; Webster-Clark, D. J., Managing bias in ROC curves. J Comput Aided Mol Des 2008, 22, (3-4), 141-6.
[3] Truchon, J. F.; Bayly, C. I., Evaluating virtual screening methods: good and bad metrics for the "early recognition" problem. J Chem Inf Model 2007, 47, (2), 488-508.
```
