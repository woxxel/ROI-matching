from Sheintuch_matching_prob import *

SH = Sheintuch_matching('/media/wollex/Analyze_AS1/linstop','246',footprints_file='results_OnACID.mat',sessions=range(1,64),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
SH.run_matching()


SH = Sheintuch_matching('/media/wollex/Analyze_AS1/linstop','231',footprints_file='results_OnACID.mat',sessions=range(1,106),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
SH.run_matching()


SH = Sheintuch_matching('/media/wollex/Analyze_AS1/linstop','232',footprints_file='results_OnACID.mat',sessions=range(13,98),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
SH.run_matching()


SH = Sheintuch_matching('/media/wollex/Analyze_AS1/linstop','236',footprints_file='results_OnACID.mat',sessions=range(1,29),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
SH.run_matching()


SH = Sheintuch_matching('/media/wollex/Analyze_AS1/Shank','918shKO',footprints_file='results_OnACID.mat',sessions=range(1,29),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
SH.run_matching()


SH = Sheintuch_matching('/media/wollex/Analyze_AS1/Shank','931wt',footprints_file='results_OnACID.mat',sessions=range(1,29),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
SH.run_matching()


SH = Sheintuch_matching('/media/wollex/Analyze_AS1/Shank','943shKO',footprints_file='results_OnACID.mat',sessions=range(1,29),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
SH.run_matching()


#SH = Sheintuch_matching('/media/wollex/Analyze_AS1/others','756',footprints_file='results_OnACID.mat',sessions=range(1,31),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
#SH.run_matching()

#SH = Sheintuch_matching('/media/wollex/Analyze_AS1/others','757',footprints_file='results_OnACID.mat',sessions=range(1,31),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
#SH.run_matching()

#SH = Sheintuch_matching('/media/wollex/Analyze_AS1/others','758',footprints_file='results_OnACID.mat',sessions=range(1,31),d_thr=12,nbins=100,use_kde=True,SNR_thr=1,r_thr=0,model='new')
#SH.run_matching()
