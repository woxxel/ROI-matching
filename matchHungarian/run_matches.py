

def run_matches(w):
  
  #results = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(1,1),thr_cost=0.7,w=w,pl=False);
  #results = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(2,2),thr_cost=0.7,w=w,pl=False);
  #results = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(3,3),thr_cost=0.7,w=w,pl=False);
  
  #_ = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=None,thr_cost=0.5,w=w,pl=False);
  _ = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(1,1),thr_cost=0.5,w=w,pl=False);
  _ = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(2,2),thr_cost=0.5,w=w,pl=False);
  _ = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(3,3),thr_cost=0.5,w=w,pl=False);
  
  _ = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=None,thr_cost=0.8,w=w,pl=False);
  _ = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(1,1),thr_cost=0.8,w=w,pl=False);
  _ = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(2,2),thr_cost=0.8,w=w,pl=False);
  _ = match_ROIs_test("/media/wollex/Analyze_AS3/Data/879/",sessions=(1,15),std=(3,3),thr_cost=0.8,w=w,pl=False);
  