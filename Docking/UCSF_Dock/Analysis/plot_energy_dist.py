import sys, os
import matplotlib  # must import first
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np

def main():
	pwd = os.getcwd()+'/'
	infile = pwd + sys.argv[1]
	plot_in = open(infile, 'r')
	read_plt = plot_in.readlines()
	plot_in.close()

	e_stat = []
	vander = []
	polard = []
	totals = []	

	for line in read_plt:
		splitline = line.strip().split()
		e_stat_NRG =float(splitline[-10])
		vander_NRG =float(splitline[-8])
		polard_NRG =float(splitline[-7])
		totals_NRG =float(splitline[-1])
		e_stat.append(e_stat_NRG)
		vander.append(vander_NRG)
		polard.append(polard_NRG)
		totals.append(totals_NRG)
		
	#print('The electrostatic energies are:'+str(sorted(e_stat)))
        #print('The van der waals energies are:'+str(sorted(vander)))
        #print('The polar desolvation energies are:'+str(sorted(polard)))
        #print('The total energies are:'+str(totals))
	es = sorted(e_stat)
	vdw = sorted(vander)
	lig_des = sorted(polard)
	tot = sorted(totals)
	
	number_of_bins = 30
	labels = ['Electrostatics', 'Van Der Waals', 'Polar Desolvation']
	data_sets = [e_stat, vander, polard]
	hist_range = (np.min(data_sets), np.max(data_sets))
	binned_data_sets = [
		np.histogram(d, range=hist_range, bins=number_of_bins)[0]
		for d in data_sets
	]	
	binned_maximums = np.max(binned_data_sets, axis=1)
	x_locations = np.arange(-50, sum(binned_maximums), np.max(binned_maximums))
	bin_edges = np.linspace(hist_range[0], hist_range[1], number_of_bins + 1)
	centers = 0.5 * (bin_edges + np.roll(bin_edges, 1))[:-1]
	heights = np.diff(bin_edges)
	
	fig, ax = plt.subplots()
	for x_loc, binned_data in zip(x_locations, binned_data_sets):
		lefts = x_loc - 0.5 * binned_data
		ax.barh(centers, binned_data, height=heights, left=lefts)

	ax.set_xticks(x_locations)
	ax.set_xticklabels(labels)

	ax.set_ylabel("Energy Value")
	ax.set_xlabel("Energy Term")

	
	plt.savefig('energy_distributions.png', dpi=300)














	#num_bins = 25
	#fig, ax = plot.subplots()

	#n, bins, patches = plt.hist(es, num_bins, facecolor='blue', alpha=0.5)
	
	#ax.plot(bins)
	#ax.set_xlabel('Energies')
	#ax.set_ylabel('Frequency')
	
	#histogram = plt.figure()
	#bins =np.linspace(-80,30, 100)
	#plt.hist(totals, facecolor='blue')
	#plt.hist(vdw, facecolor='purple')
	#plt.hist(lig_des, facecolor='orange')
	#plt.plot(es, facecolor='green')
	#plt.hist(e_stat, density = 1, bins = 10)
	#plt.axis([-80,30, 0, 1])
	#plt.title("Energies")
	#plt.xlabel("Value")
	#plt.ylabel("Frequency")

	#fig = plt.gcf()
	#data = es	
	#plt.hist(es, bins = [-80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40])
	#plt.hist(data, bins=range(min(data), max(data) + binwidth, binwidth))
	
main()
