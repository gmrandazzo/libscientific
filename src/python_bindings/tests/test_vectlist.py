from libscientific.vectlist import *
import random

def test_vectlist():
    
    a = [[random.random() for j in range(10)] for i in range(3)]
    d = DVectorList(a)
    #d.debug()
    #print("get value")
    #vect.print_dvector(d[1])
    #print("set list")
    d[1] = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
    dlst = d.tolist()
    #print("print the list converted")
    #for item in dlst:
    #    print(item)

    #print("Append a list - high level")
    lst = [9, 8, 7, 6]
    d.append(lst)
    #print("Reappend ppend a list - low level")
    dvector_list_append(d.dvl, lst)
    #d.debug()
