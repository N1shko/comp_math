from base_1 import *
from zhen import *

#import advanced_1
import random
import matplotlib.pyplot as plt
import numpy as np

def main():
    #basic part
    part = input()

    if(part=='base'):
        uniform()
        chebishev()
   # elif(part == 'a_tsk_1'):
     #   advanced_part_tsk_1()
    #elif (part == 'a_tsk_2'):

    elif (part == 'tsk1'):

        tsk1_2(True)
    elif(part == 'tsk3'):
        tsk3(tsk1_2(False))
    elif(part == "all"):
        tsk3(tsk1_2(True))
    #advanced part



if __name__ == '__main__':
    #np.random.seed(1000)

    main()
