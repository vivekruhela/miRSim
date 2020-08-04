#!/usr/bin/python

from scipy.stats import gamma
from scipy.stats import uniform
from scipy.stats import expon
from scipy.stats import poisson

def expression_split(number, number_of_subsections, distribution, seed, min_random_number_desired = 2):
    split_number_list = []
    cumulative_sum_of_random_numbers = 0
    current_subsection = 1
    max_random_number = int(number/number_of_subsections)
    if min_random_number_desired < number:
        if min_random_number_desired > max_random_number:
#             print("WARNING: Cannot have min number as {} and split {} in {} subsections".format(min_random_number_desired, number, number_of_subsections))
            number_of_subsections = int(np.floor(number/min_random_number_desired))
            return expression_split(number,number_of_subsections,distribution, seed, min_random_number_desired)
                    
        elif distribution == 'uniform':
            split_num1 = uniform.rvs(size=number_of_subsections, loc = 10, scale=20, random_state= seed)
        elif distribution == 'gamma':
            split_num1 = gamma.rvs(a=5, size=number_of_subsections, random_state= seed)
        elif distribution == 'exponential':
            split_num1 = expon.rvs(scale=1,loc=0,size=number_of_subsections, random_state= seed)
        elif distribution == 'poisson':
            split_num1 = poisson.rvs(mu=3, size=number_of_subsections, random_state= seed)
            
        split_num1 = [int(number*v/sum(split_num1)) for v in split_num1]
        num_test = [1 for v in split_num1 if v <= min_random_number_desired]
        if sum(num_test) > int(0.25*len(split_num1)):
            return expression_split(number,int(number_of_subsections*0.75),distribution, seed, min_random_number_desired)
        for i in range(len(split_num1)):
            if split_num1[i] < min_random_number_desired:
                split_num1[i] = min_random_number_desired                    
        split_num1[-1] = number-sum(split_num1[:-1])
        if sum([1 for v in split_num1 if v<=0]) >= 1:
            return expression_split(number,int(number_of_subsections*0.75),distribution, seed, min_random_number_desired)
#             print('ERROR: Please lower the minimum Depth.')
#             sys.exit(1)
        return split_num1
    else:
        print('ERROR : minimum depth is greater than provided number and can not be splitted.')
        sys.exit(1)
