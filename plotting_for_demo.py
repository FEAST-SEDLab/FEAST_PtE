import matplotlib.pyplot as plt
import feast
feast.ResultsProcessing.plotting_functions.abatement_cost_plotter('ExampleRunScriptResults')

ax = plt.gca()
ax.set_xticklabels(["No dispatch\nthreshold", "10%\naccuracy", "50%\naccuracy", "200%\naccuracy"])
plt.show()