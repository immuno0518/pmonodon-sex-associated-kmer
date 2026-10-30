from release_config import *
from reviewer_genomics_common import *
import os
os.environ['MPLCONFIGDIR']=str(ROOT/'work/matplotlib_cache')
import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
LD=OUT/'local_LD'
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(8.2,8.5));fig.subplots_adjust(left=.09,right=.91,bottom=.20,top=.90,hspace=.36,wspace=.40)
cmap=copy.copy(plt.cm.viridis);cmap.set_bad('#E9ECEF')
summary=json.loads((LD/'summary.json').read_text())
for ax,group,title,letter in zip(axes.flat[:3],('all200','female100','male100'),('All individuals (n = 200)','Females (n = 100)','Males (n = 100)'),('A','B','C')):
    rows=list(csv.DictReader((LD/(group+'_DP5_heatmap_10kb.tsv')).open(),delimiter='\t'));m=np.full((60,60),np.nan)
    for r in rows:
        i=(int(r['bin1_start'])-2500001)//10000;j=(int(r['bin2_start'])-2500001)//10000
        if int(r['pair_count'])>=10:m[i,j]=m[j,i]=float(r['mean_r2'])
    im=ax.imshow(m,origin='lower',extent=(2.5,3.1,2.5,3.1),vmin=0,vmax=1,cmap=cmap,interpolation='nearest',aspect='equal')
    ax.add_patch(Rectangle((2.71,2.71),.18,.18,fill=False,edgecolor='#E8A52D',linewidth=1.3))
    ax.set_xlabel('Chr19 position (Mb)');ax.set_ylabel('Chr19 position (Mb)');ax.set_title(title+'\n'+str(summary['results']['5']['groups'][group]['informative_sites'])+' informative SNPs',fontsize=9)
    ax.text(-.20,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=13)
cax=fig.add_axes([.935,.57,.015,.25]);fig.colorbar(im,cax=cax,label='Mean unphased genotype r²')
ax=axes.flat[3]
for group,label,color in [('all200','All 200','#495A68'),('female100','Female 100','#B53C7D'),('male100','Male 100','#2375AD')]:
    rr=list(csv.DictReader((LD/(group+'_DP5_distance_10kb.tsv')).open(),delimiter='\t'))
    xx=[(int(r['distance_start_bp'])+5000)/1000 for r in rr];yy=[float(r['mean_r2']) if r['mean_r2']!='NA' else np.nan for r in rr]
    ax.plot(xx,yy,label=label,color=color,lw=1.3)
ax.set(xlabel='SNP separation (kb)',ylabel='Mean unphased genotype r²',ylim=(0,1),xlim=(0,600));ax.legend(frameon=False,fontsize=8);ax.set_title('Distance profile (10-kb bins)',fontsize=9);ax.text(-.20,1.08,'D',transform=ax.transAxes,fontweight='bold',fontsize=13)
fig.suptitle('Local genotype correlations in the full-sib family',fontsize=13,y=.972)
fig.text(.09,.055,'DP ≥ 5, GQ ≥ 20; per-sex call rate ≥ 90%; MAF ≥ 0.05. Golden box: Chr19:2.71–2.89 Mb.\nHeatmaps summarize eligible pairs in 10-kb coordinate bins; gray: fewer than 10 valid pairs.\nEach group uses its informative sites; matched-pair comparisons are supplied separately.\nDescriptive family-level correlations do not define population LD or a recombination boundary.',fontsize=8,linespacing=1.6)
fig.savefig(LD/'local_LD_family_stratified.pdf');fig.savefig(LD/'local_LD_family_stratified.svg');fig.savefig(LD/'local_LD_family_stratified.png',dpi=300);plt.close(fig)
print('LD figure written')
