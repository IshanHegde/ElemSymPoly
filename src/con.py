import pyrasch
import pandas as pd
import numpy as np
from scipy.stats import norm
import random
import time
df = pd.read_csv("../G0.txt",sep="\t",index_col="Student")
arr =df.to_numpy(dtype=float)

#print(arr)

NUMBER_PERSONS=1000
NUMBER_ITEMS=5


question_diff= np.linspace(-3.0,3.0,NUMBER_ITEMS)

abilities= np.random.normal(0,1,NUMBER_PERSONS)

rng = np.random.default_rng()
max_scores=rng.integers(1,6,(NUMBER_ITEMS,))

def prob_func(i,x):

    means=pcc_means[i]
    raw_probs=[norm.pdf(x,loc=mean,scale=1) for mean in means]
    normalized_probs=raw_probs/np.sum(raw_probs)
    return normalized_probs*100



pcc_means= []
for i in range(NUMBER_ITEMS):
    max= max_scores[i]
    
    mean_points=np.random.normal(0,3,max+1)

    ordered_mean_points= np.sort(mean_points)
    
    pcc_means.append(ordered_mean_points)



data=[]
for n in range(NUMBER_PERSONS):

    person_n=[]
    for i in range(NUMBER_ITEMS):

        probs=prob_func(i,abilities[n])
        person_n.append(random.choices(range(0,max_scores[i]+1),probs)[0])
    
    data.append(person_n)

#print(data)

start= time.time()

a = pyrasch.Rasch(np.array(data))


a.JMLE(4)

print(np.array(data))
end= time.time()

print(end-start)

#print(a.difficulty)
#print(question_diff)
#print(a.ability)
#print(abilities)
#print(a.ability)
#print(a.RA_Thresholds)




