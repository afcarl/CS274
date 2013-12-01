import knn
import matplotlib.pyplot as plt
from pylab import *

k_values = [1,5,10,15,30,50]
p = 0.5
for k in k_values:
	knn.main("ALL.dat", "AML.dat", k, p)
	print "\n"


# plot(k_values, [0.93, 0.94, 0.83, 0.79, 0.61, 0.61])
# xlabel('k')
# ylabel('accuracy')
# savefig('question1.jpg')

print "\n\n---------------------------------\n\n"

p_values = [0, 0.05, 0.20, 0.50, 0.75, 0.95, 1.00]
results = []
k = 30
for p in p_values:
	#knn.main("ALL.dat", "AML.dat", k, p)
	print "\n"

fig, ax = plt.subplots()
sensitivity = [1.0, 1.0, 1.0, 1.0, 0.77, 0.0, 0.0]
specificity = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
to_plot_spec = [1 - x for x in specificity]
ax.plot(sensitivity, to_plot_spec)

sens_old = spec_old = -1.0
for i in range(7):
	sens_val = sensitivity[i]
	spec_val = to_plot_spec[i]
	if sens_val == sens_old and spec_val == spec_old:
		sens_old = sens_val
		spec_old = spec_val
		continue
	ax.annotate("p=" + str(p_values[i]), xy=(sensitivity[i],to_plot_spec[i]), xytext=(-5, 5),
                textcoords='offset points')
	sens_old = sens_val
	spec_old = spec_val

xlabel('sensitivity')
ylabel('1 - specificity')
# savefig('question2.jpg')
# plt.show()
