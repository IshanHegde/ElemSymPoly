import numpy as np
from utils.functions import sigmoid

def PROX(arr,dichotomous=True,PROX_MAX=100):

    ability=np.zeros(arr.shape[0])
    difficulty=np.zeros(arr.shape[1])

    if dichotomous:
        MAX_PERSON_SCORE=arr.shape[1]
        MAX_ITEM_SCORE=arr.shape[0]

    condition=1
    iter_count=0
    while condition:

        for person in range(len(ability)):

            person_raw_score=np.sum(arr[person,:])

            if person_raw_score==0:
                person_raw_score+=1
            elif person_raw_score==MAX_PERSON_SCORE:
                person_raw_score-=1

            ability[person] = np.mean(difficulty)+ np.sqrt(1+np.var(difficulty)/2.9)*np.log(person_raw_score/(MAX_PERSON_SCORE-person_raw_score))

        for item in range(len(difficulty)):
                
            item_raw_score=np.sum(arr[:,item])

            difficulty[item] = np.mean(ability)- np.sqrt(1+np.var(ability)/2.9)*np.log(item_raw_score/(MAX_ITEM_SCORE-item_raw_score))

        iter_count+=1
        if iter_count>PROX_MAX:
            condition=0

    return ability,difficulty


if __name__ =='__main__':

    synthetic_questions = np.arange(-3, 3, 0.3)
    synthetic_students = np.random.normal(loc=0,scale=1.5,size=(100,))
    synthetic_logits = synthetic_students.reshape(-1,1) - synthetic_questions.reshape(1,-1)
    synthetic_probs = sigmoid(synthetic_logits)
    synthetic_data = (synthetic_probs > np.random.rand(synthetic_probs.shape[0],synthetic_probs.shape[1])).astype('int')

    sample_arr = synthetic_data

    print(sample_arr)

    print(PROX(sample_arr)[1])
