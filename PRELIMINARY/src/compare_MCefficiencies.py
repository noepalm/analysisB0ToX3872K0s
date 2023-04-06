import csv
import ROOT
import pandas as pd
import argparse
import matplotlib.pyplot as plt

usage = ' usage: $python ./src/compare_MCefficiencies.py -f [file.csv] --out [outdir]'
parser = argparse.ArgumentParser(usage = usage)

parser.add_argument('-p', '--path', default = './data/efficiency_X3872.csv')
parser.add_argument('-o', '--outdir', default = '/eos/user/c/cbasile/B0toX3872K0s/')


args = parser.parse_args()

data = pd.read_csv(args.path, index_col = 0, decimal=',')
data = data.fillna(0)
print(data.dtypes)

#df = ROOT.RDF.MakeCsvDataFrame(args.path)
#dataset_names = df.GetColumnNames()
#h_prof = df.Profile1D(("eff_16", "", 5, 0, 5.), "2017", "2016")
#c = ROOT.TCanvas("c", "", 800, 600)
#h_prof.Draw()
#c.SaveAs("/eos/user/c/cbasile/www/B0toX3872K0s/RECO_LEVEL/MCefficiencies_df.png");

years = ['2016preVFP', '2016', '2017', '2018']
colors_per_year = {
    '2016preVFP' : 'green', 
    '2016' : 'gold', 
    '2017' : 'blue', 
    '2018' : 'red'} 

fig, ax = plt.subplots(figsize=(10,8))
[data[y].plot(drawstyle="steps-mid",  yerr = data[y+'_err'], linewidth = 2, barsabove = True, color = colors_per_year[y]) for y in years]
plt.xlim(-0.5,4.5)
ax.set_xticklabels(['','Gen filter', 'Preselection', 'HLT-emulation', 'Signal Region\nJpsiPiPi K0s', 'MVA'])
plt.margins(x=1)
plt.grid(axis='x')
plt.grid(axis='y')
plt.legend(frameon=False, fontsize = 20)
plt.xlabel('selection-step', fontweight='bold', fontsize='16', labelpad = 10, horizontalalignment='center')
plt.xticks(fontweight='bold', fontsize='12', horizontalalignment='center')
plt.ylabel('efficiency', fontweight='bold', fontsize='16', labelpad = 10, horizontalalignment='center')
plt.yticks(fontweight='bold', fontsize='12') 
plt.yscale('log')
plt.savefig(args.outdir + ".png")
plt.savefig(args.outdir + ".pdf")
#plt.show()
